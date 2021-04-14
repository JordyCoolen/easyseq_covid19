#!/usr/bin/env nextflow

params.threads = 4
params.outDir = "./output"
params.reads = "$baseDir/test/test_OUT01_R{1,2}.fastq.gz"

//filters
call_threshold = 0.40
qual_threshold = 20
min_depth = 10
max_ambig = 0.09999

// Parsing the input parameters
sampleName       = "$params.sampleName"
outDir           = "$params.outDir"
threads          = "$params.threads"
reference  = "$baseDir/db/data/NC_045512.2/NC_045512.2.fasta"
primerfile = "$baseDir/db/primers_bedpe_V3.bed"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .set  { reads_ch1 }

// Tools paths and command prefixes
fastq_reporter   = "$baseDir/scripts/fastq_report.py"
reporter         = "$baseDir/scripts/final_report.py"
vcf2table        = "$baseDir/scripts/vcf2table.py"
parse_stats      = "$baseDir/scripts/parse_stats.py"

log.info """

NEXTFLOW EasySeq RC-PCR SARS-CoV-2/COVID-19
Variant pipeline V0.7.0
================================
sample     : $params.sampleName
reads      : $params.reads
outDir     : $params.outDir
codeBase   : $baseDir
threads    : $params.threads

~~~~~~~~~~variant filter~~~~~~~~~~~~
call_threshold: $call_threshold
qual_threshold: $qual_threshold
min_depth: $min_depth

~~~~~~~~~~lineage filter~~~~~~~~~~~~
max_ambig: $max_ambig

~~~~~~~~~~~Databases~~~~~~~~~~~
reference  : $reference
primerfile : $primerfile

~~~~~~~~~~~Authors~~~~~~~~~~~~~~
        J.P.M. Coolen
        R.A. Lammerts
          J.T. Vonk
================================
"""

// Clean reads (adapter and read length filter)
process '1A_clean_reads' {
    tag '1A'
    conda 'bioconda::fastp=0.20.1 bioconda::pyfastx=0.6.12 conda-forge::simplejson=3.17.0'
    publishDir outDir + '/QC', mode: 'copy', pattern: "*.fastp.json"
    input:
        set pairID, file(reads) from reads_ch1
    output:
        set file("${reads[0].baseName}_fastp.fastq.gz"), file("${reads[1].baseName}_fastp.fastq.gz") into fastp_2A
        file "${sampleName}.fastp.json"
        file ".command.*"
    script:
        """
        fastp -i ${reads[0]} -I ${reads[1]} -o ${reads[0].baseName}_fastp.fastq.gz -O ${reads[1].baseName}_fastp.fastq.gz \
        --trim_poly_x --length_required 100 --json ${sampleName}.fastp.json --html ${sampleName}.fastp.html --thread ${threads}
        """
}

// Process 2A: Map the reads to the reference genome
process '2A_map_paired_reads' {
    tag '2A'
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.11 bioconda::bamclipper=1.0.0'
    publishDir outDir + '/mapping', mode: 'copy'
    input:
        set file(trimmed1), file(trimmed2) from fastp_2A
    output:
        file("${sampleName}.bam") into bam_2B
        file("${sampleName}.bam.bai") into bamindex_2B
        file ".command.*"
    script:
        """
        bwa mem -M -k 10 -t ${threads} ${reference} ${trimmed1} ${trimmed2} | samtools sort -@${threads} -o ${sampleName}.bam -
        samtools index ${sampleName}.bam > ${sampleName}.bam.bai
        """
}

// Process 2B: Filter the bam file of primer sequences
process '2B_bam_clipper' {
    tag '2B'
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.11 bioconda::bamclipper=1.0.0'
    publishDir outDir + '/mapping', mode: 'copy'
    input:
        file bam from bam_2B
        file bamindex from bamindex_2B
    output:
        file("${sampleName}.final.bam") into (bam_2C, bam_3A, bam_3D)
        file("${sampleName}.final.bam.bai") into bamindex_2C
        file ".command.*"
    script:
        """
        bamclipper.sh -b ${bam} -n ${threads} -p $primerfile

        samtools sort ${sampleName}.primerclipped.bam -o ${sampleName}.final.bam
        samtools index ${sampleName}.final.bam
        """
}

// Process 2C: Genome depth
process '2C_depth' {
    tag '2C'
    conda 'bioconda::mosdepth=0.3.1'
    publishDir outDir + '/QC', mode: 'copy'
    input:
        file bam from bam_2C
        file bamindex from bamindex_2C
    output:
        file "*"
        file ".command.*"
    script:
        """
        mosdepth --threads $threads ${sampleName} ${bam}
        """
}

// Process 3A: Variant calling
process '3A_variant_caller' {
    tag '3A'
    conda "${baseDir}/conda/env-lofreq/"
    publishDir outDir + '/rawvcf', mode: 'copy'
    input:
        file bam from bam_3A
    output:
        file("${sampleName}.raw.vcf") into (rawvcf_3B, rawvcf_3C, rawvcf_3F)
        file ".command.*"
    script:
        """
        lofreq indelqual --dindel --ref ${reference} -o ${sampleName}.lofreq.bam ${bam}
        samtools index ${sampleName}.lofreq.bam

        lofreq call-parallel --call-indels --no-default-filter --pp-threads ${threads} -f ${reference} --min-cov ${min_depth} \
        -o ${sampleName}.raw.vcf ${sampleName}.lofreq.bam
        """
}

// Process 3B: filter_variants
// Filter out variants that do not adhere to the criteria, outputs plain vcf
// criteria used: too low overall DP, too low quality, percentage of alts over used depth
process '3B_filter_variants' {
    tag '3B'
    conda "${baseDir}/conda/env-lofreq/"
    publishDir outDir + '/vcf', mode: 'copy'
    input:
        file vcf from rawvcf_3B
    output:
        file("${sampleName}.variants.vcf") into (vcf_3E, vcf_3D, vcf_5A)
        file ".command.*"
    script:
        """
        lofreq filter --no-defaults --af-min ${call_threshold} -i ${vcf} -o ${sampleName}.variants.vcf
        """
}

// Process 3C: ubiquitous_variant
// Find the variants where we are uncertain that the variation is correct
// so that they can be replaced/masked in the consensus seq by N
process '3C_ubiquitous_variant' {
    tag '3C'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/uncovered', mode: 'copy'
    input:
        file vcf from rawvcf_3C
    output:
        file "${sampleName}_ubiq.bed"
        file ".command.*"
    script:
        """
		bcftools query -f "%CHROM\\t%POS\\t%END\\n" -i \
		"(%QUAL<${qual_threshold} || \
		%MAX(FORMAT/AD[0:1])/%MAX(FORMAT/DP)<${call_threshold}) && \
		INFO/DP>=${min_depth}" ${vcf} | awk '{{print(\$1 \"\\t\" \$2 \"\\t\" \$3 \"\\tubiquitous_variant\")}}' > ${sampleName}_ubiq.bed
        """
}

// Process 3D: non_covered_regions
// Find the regions where coverage was too low
// so that they can be replaced/maksed in the consensus seq by N
// uses bedtools subtract to make sure deletions are not masked
process '3D_non_covered_regions' {
    tag '3D'
    conda 'bioconda::bedtools'
    publishDir outDir + '/uncovered', mode: 'copy'
    input:
        file bam from bam_3D
        file vcf from vcf_3D
    output:
        file("${sampleName}_noncov.bed") into (noncov_4A)
        file ".command.*"
    script:
        """
		bedtools genomecov -ibam ${bam} -bga | \
        awk '\$4 < ${min_depth}' | \
        awk '{{print(\$1 \"\\t\" \$2 + 1 \"\\t\" \$3 \"\\tlow_coverage\")}}' |\
        bedtools subtract -a - -b ${vcf} > ${sampleName}_noncov.bed
        """
}

// Process 3E: index_called_vcf
// compresses and indexes the vcf file for making the consensus sequence
process '3E_index_called_vcf' {
    tag '3E'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/vcf', mode: 'copy'
    input:
        file vcf from vcf_3E
    output:
        file("${sampleName}.final.vcf")
        file("${sampleName}.final.vcf.gz") into (compressed_vcf_4A)
        file("${sampleName}.final.vcf.gz.csi") into (compressed_vcf_index_4A)
        file ".command.*"
    script:
        """
        sed '5i ##contig=<ID=NC_045512.2,length=29903>' ${vcf} > ${sampleName}.final.vcf
        bcftools view ${sampleName}.final.vcf -O b > ${sampleName}.final.vcf.gz
        bcftools index --threads ${threads} ${sampleName}.final.vcf.gz -o ${sampleName}.final.vcf.gz.csi
        """
}


// Process 3F: filter_variants notpassed
// Same as BCF filter, but in reverse. this is to report all variants that didn't pass the filters
// outputs in a tabular format
// WARNING: if filters are ever removed or added, they should change here as well
process '3F_filter_variants_notpassed' {
    tag '3F'
    conda "${baseDir}/conda/env-lofreq/"
    publishDir outDir + '/vcf/notpassed', mode: 'copy'
    input:
        file vcf from rawvcf_3F
    output:
        file("${sampleName}.notpassed.vcf")
        file ".command.*"
    script:
        """
        lofreq filter --no-defaults --af-max ${call_threshold} -i ${vcf} -o ${sampleName}.notpassed.vcf
        """
}

// 4A: create the consensus fasta
// variants that are marked as uncertain or low covered regions are masked by N
process '4A_create_consensus' {
    tag '4A'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/report', mode: 'copy'
    input:
        file vcf from compressed_vcf_4A
        file vcf_index from compressed_vcf_index_4A
        file noncov from noncov_4A
    output:
        file("${sampleName}.fasta") into (consensus_6, consensus_7B)
        file ".command.*"
  script:
        """
		cat ${reference} | bcftools consensus ${vcf} \
		-m <(cat ${noncov}) \
		${vcf} > ${sampleName}.fasta
        sed -i 's/^>NC_045512.2/>${sampleName}/g' ${sampleName}.fasta
        """
}

// 5A: annotation of the genome/consensus fasta
process '5A_annotation' {
    tag '5A'
    conda 'bioconda::snpeff=5.0'
    publishDir outDir + '/annotation', mode: 'copy'
    publishDir outDir + '/QC', mode: 'copy', pattern: "${sampleName}_snpEff.csv"
    input:
        file vcf from vcf_5A
    output:
        file("${sampleName}_snpEff.csv") into snpEffStats_7
        file("${sampleName}_annot_table.txt") into annotation_8
        file("${sampleName}_snpEff.genes.txt")
        file("snpEff_summary.html")
        file ".command.*"
  script:
        """
        snpEff eff -v NC_045512.2 ${vcf} -ud 0 -strict -hgvs1LetterAa \
        -c ${baseDir}/db/snpEff.config -csvStats ${sampleName}_snpEff.csv \
        > ${vcf.baseName}_annot.vcf
        ${baseDir}/conda/env-variantcalling/bin/python $vcf2table ${vcf.baseName}_annot.vcf --sample ${sampleName} \
        -ad -e -o ${sampleName}_annot_table.txt
        """
}

// 6: pangolin to determine lineage
process '6_lineage' {
    tag '6'
    conda "${baseDir}/conda/env-pangolin"
    publishDir outDir + '/lineage', mode: 'copy'
    input:
        file consensus from consensus_6
    output:reporter
        file("lineage_report.csv") into lineage_8
        file ".command.*"
  script:
        """
        pangolin ${consensus} --max-ambig ${max_ambig}
        """
}

// 7A: multiqc
process '7A_QC' {
    tag '7A'
    conda "bioconda::multiqc"
    publishDir outDir + '/QC', mode: 'copy'
    input:
        file snpEffStats from snpEffStats_7
    output:reporter
        file "*"
        file ".command.*"
  script:
        """
        multiqc ${outDir}/QC
        """
}

// 7B: genome completeness/stats calculation
process '7B_genome_stats' {
    tag '7B'
    conda "bioconda::ucsc-facount"
    publishDir outDir + '/QC', mode: 'copy'
    input:
        file consensus from consensus_7B
    output:
        file "stats.txt" into stats_8B
        file ".command.*"
  script:
        """
        facount ${consensus} > stats.txt
        python $parse_stats --stats stats.txt
        """
}

// Process 8A: obtain run parameters
process '8A_parameters' {
    tag '8A'
    publishDir outDir + '/report', mode: 'copy'
    input:
    output:
        file "parameters.txt" into params_8B
    script:
        """
        touch parameters.txt
        echo "Parameter\tValue" >> parameters.txt
        echo "mutation frequency:\t>=${call_threshold}" >> parameters.txt
        echo "QUAL:\t>=${qual_threshold}" >> parameters.txt
        echo "Minimum sequence depth:\t${min_depth}" >> parameters.txt
        echo "Maximum %N:\t${max_ambig}" >> parameters.txt
        """
}

// Process 8B: generate a report for interpretation by the clinician (or for research purposes)
process '8B_report' {
    tag '8B'
    conda "${baseDir}/conda/env-025066a104bf8ce5621e328d8009733a"
    publishDir outDir + '/report', mode: 'copy'
    input:
        file lineage from lineage_8
        file annotation from annotation_8
        file params from params_8B
        file stats from stats_8B
    output:
        file "${sampleName}.html"
        file "${sampleName}.pdf"
        file ".command.*"
    script:
        """
        $reporter --sampleName ${sampleName} --lineage ${lineage} \
        --annotation ${annotation} --params ${params} --stats ${stats}
        """
}
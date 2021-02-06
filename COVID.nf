#!/usr/bin/env nextflow

params.threads = 4
params.outDir = "./output"
params.beta = false
params.subsampling = true
params.max_depth = 80
params.min_depth = 30
params.max_bad_cov = 0.15
params.reads = "$baseDir/test/test_OUT01_R{1,2}.fastq.gz"

//bcf_filter
call_treshold = 0.75
qual_treshold = 20
min_depth = 5

// Parsing the input parameters
sampleName       = "$params.sampleName"
outDir           = "$params.outDir"
subsampling      = "$params.subsampling"
threads          = "$params.threads"
//reference        = "params.reference"
reference  = "$baseDir/db/data/NC_045512.2/NC_045512.2.fasta"
primerfile = "$baseDir/db/primers_bedpe.bed"
kma_index  = "$baseDir/db/kma/HV69-70"
vcf_db     = "$baseDir/db/kma/vcf/"

// special channel for the fastQ reads
Channel
      .fromFilePairs( params.reads )
      .ifEmpty { "cannot find read pairs in path"}
      .set  { reads_ch1 }

// Tools paths and command prefixes
//picard           = "java -jar $baseDir/bin/picard.jar"
//snpit            = "/miniconda/envs/mycoprofiler/bin/snpit" // https://github.com/philipwfowler/snpit
//snpit_reporter   = "$baseDir/bin/MycoprofilerUtils/snpit_wrapper.py"
//species_reporter = "$baseDir/bin/MycoprofilerUtils/species_report.py"
fastq_reporter   = "$baseDir/scripts/fastq_report.py"
reporter         = "$baseDir/scripts/final_report.py"
merge_json       = "$baseDir/scripts/merge_json.py"
vcf2table        = "$baseDir/scripts/vcf2table.py"
HV69_70          = "$baseDir/scripts/HV69-70.py"

log.info """

NEXTFLOW EasySeq RC-PCR SARS-CoV-2/COVID-19
Variant pipeline V0.3
================================
sample     : $params.sampleName
reads      : $params.reads
outDir     : $params.outDir
codeBase   : $baseDir
threads    : $params.threads

~~~~~~~~~~bcf_filter~~~~~~~~~~~~
call_treshold: $call_treshold
qual_treshold: $qual_treshold
min_depth: $min_depth

~~~~~~~~~~~Databases~~~~~~~~~~~
reference  : $reference
primerfile : $primerfile
kma_index  : $kma_index

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
        set file("${reads[0].baseName}_fastp.fastq.gz"), file("${reads[1].baseName}_fastp.fastq.gz") into (fastp_2A, fastp_2D)
        file "${sampleName}.fastp.json"
        //file "${sampleName}_*.json" into meta_json_fastq
        file ".command.*"
    script:
        """
        #$fastq_reporter -i ${reads[0]} --sample ${sampleName}

        fastp -i ${reads[0]} -I ${reads[1]} -o ${reads[0].baseName}_fastp.fastq.gz -O ${reads[1].baseName}_fastp.fastq.gz \
        --trim_poly_x --length_required 100 --json ${sampleName}.fastp.json --html ${sampleName}.fastp.html --thread ${threads}

        # merge json results to meta json
        #$merge_json --json1 ${sampleName}_*.json --json2 ${sampleName}.fastp_json --key fastp
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

// Process 2B: Map the reads to the reference genome
process '2B_bam_clipper' {
    tag '2B'
    conda 'bioconda::bwa=0.7.17 bioconda::samtools=1.11 bioconda::bamclipper=1.0.0'
    publishDir outDir + '/mapping', mode: 'copy'
    input:
        file bam from bam_2B
        file bamindex from bamindex_2B
    output:
        file("${sampleName}.primerclipped.bam") into (bam_2C, bam_3A, bam_3E)
        file("${sampleName}.primerclipped.bam.bai") into bamindex_2C
        file ".command.*"
    script:
        """
        bamclipper.sh -b ${bam} -n ${threads} -p $primerfile
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

// annotation of genome
process '2D_HV69-70' {
    tag '2D'
    conda 'bioconda::kma=1.3.9'
    publishDir outDir + '/HV69-70', mode: 'copy'
    input:
        set file(trimmed1), file(trimmed2) from fastp_2D
    output:
        file "${sampleName}_HVdel.vcf" into HVdel_3C
        file("*")
        file ".command.*"
  script:
        """
        kma -ipe ${trimmed1} ${trimmed2} -t_db ${kma_index} -o ./${sampleName}

        # will give a stdout
        python $HV69_70 -i ./${sampleName}.res -p ${vcf_db} > ${sampleName}_HVdel.vcf
        """
}

// Process 3A: Variant calling
process '3A_variant_caller' {
    tag '3A'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/rawvcf', mode: 'copy'
    input:
        file bam from bam_3A
    output:
        file("${sampleName}.vcf.gz") into (rawvcf_3B, rawvcf_3D, rawvcf_3G) //, rawvcf_4A
        //file("${sampleName}.vcf.gz.csi") into (rawvcfindex_4A)
        file ".command.*"
    script:
        """
        bcftools mpileup -L 999999 -Q 0 -q 0 -A -B -d 1000000 \
        --threads ${threads} -a \'FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR\' -f ${reference} ${bam} | \

        bcftools call --threads ${threads} --no-version -c -v --ploidy 1 -Ob -o ${sampleName}.vcf.gz
        """
}

// Process 3B: filter_variants
// Filter out variants that do not adhere to the criteria, outputs plain vcf
// criteria used: too low overall DP, too low quality, percentage of alts over used depth
process '3B_filter_variants' {
    tag '3B'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/vcf', mode: 'copy'
    input:
        file vcf from rawvcf_3B
    output:
        file("${sampleName}.vcf") into (vcf_3C)
        file ".command.*"
    script:
        """
        bcftools view -Ov -i '%QUAL>=${qual_treshold} && \
        %MAX(FORMAT/AD[0:1])/%MAX(FORMAT/DP)>=${call_treshold} && INFO/DP>=${min_depth}' ${vcf} > ${sampleName}.vcf
        """
}

// Process 3C: concat_vcf
// This process merges the HV 69-70 vcf detection with the bcftools call vcf file
// This results in a correct vcf file which will be used for further processing.
process '3C_concat_vcf' {
    tag '3C'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/vcf', mode: 'copy'
    input:
        file vcf from vcf_3C
        file HVdel from HVdel_3C
    output:
        file("${sampleName}_final.vcf") into (vcf_3E, vcf_3F, vcf_5A)
        file("${sampleName}_table.txt")
        file ".command.*"
    script:
        """

        # concat VCF files
        bcftools view ${vcf} -O b > ${vcf}.gz
        bcftools view ${HVdel} -O b > ${HVdel}.gz
        bcftools index ${vcf}.gz -o ${vcf}.gz.csi
        bcftools index ${HVdel}.gz -o ${HVdel}.gz.csi

        bcftools reheader -s ${sampleName} ${vcf}.gz -h ${vcf} -o ${sampleName}_final.vcf.gz
        bcftools reheader -s ${sampleName} ${HVdel}.gz -h ${vcf} -o ${sampleName}_HVdelfinal.vcf.gz

        bcftools index ${sampleName}_final.vcf.gz -o ${sampleName}_final.vcf.gz.csi
        bcftools index ${sampleName}_HVdelfinal.vcf.gz -o ${sampleName}_HVdelfinal.vcf.gz.csi

        bcftools concat ${sampleName}_HVdelfinal.vcf.gz ${sampleName}_final.vcf.gz -a -d all -o ${sampleName}_final.vcf

        cat ${sampleName}_final.vcf | "${baseDir}/conda/env-variantcalling/bin/python" "$vcf2table" - \
        --sample ${sampleName} > ${sampleName}_table.txt
        """
}

// Process 3D: ubiquitous_variant
// Find the variants where we are uncertain that the variation is correct
// so that they can be replaced/masked in the consensus seq by N
process '3D_ubiquitous_variant' {
    tag '3D'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/uncovered', mode: 'copy'
    input:
        file vcf from rawvcf_3D
    output:
        file("${sampleName}_ubiq.bed") into (ubiq_4A)
        file ".command.*"
    script:
        """
		bcftools query -f "%CHROM\\t%POS\\t%END\\n" -i \
		"(%QUAL<${qual_treshold} || \
		%MAX(FORMAT/AD[0:1])/%MAX(FORMAT/DP)<${call_treshold}) && \
		INFO/DP>=${min_depth}" ${vcf} | awk '{{print(\$1 \"\\t\" \$2 \"\\t\" \$3 \"\\tubiquitous_variant\")}}' > ${sampleName}_ubiq.bed
        """
}

// Process 3E: non_covered_regions
// Find the regions where coverage was too low
// so that they can be replaced/maksed in the consensus seq by N
// uses bedtools subtract to make sure deletions are not masked
process '3E_non_covered_regions' {
    tag '3E'
    conda 'bioconda::bedtools'
    publishDir outDir + '/uncovered', mode: 'copy'
    input:
        file bam from bam_3E
        file vcf from vcf_3E
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

// Process 3F: index_called_vcf
// compresses and indexes the vcf file for making the consensus sequence
process '3F_index_called_vcf' {
    tag '3F'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/vcf_index', mode: 'copy'
    input:
        file vcf from vcf_3F
    output:
        file("${sampleName}_concat.vcf.gz") into (compressed_vcf_4A)
        file("${sampleName}_concat.vcf.gz.csi") into (compressed_vcf_index_4A)
        file ".command.*"
    script:
        """
        bcftools view ${vcf} -O b > ${sampleName}_concat.vcf.gz
        bcftools index --threads ${threads} ${sampleName}_concat.vcf.gz -o ${sampleName}_concat.vcf.gz.csi
        """
}


// Process 3G: filter_variants
// Same as BCF filter, but in reverse. this is to report all variants that didn't pass the filters
// outputs in a tabular format
// WARNING: if filters are ever removed or added, they should change here as well
process '3G_filter_variants_notpassed' {
    tag '3G'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/vcf/notpassed', mode: 'copy'
    input:
        file vcf from rawvcf_3G
    output:
        file("${sampleName}_3G.txt")
        file ".command.*"
    script:
        """
        # add new python to path
        #export PATH="${baseDir}/conda/env-variantcalling/bin/:$PATH"

        bcftools view -Ov -i '%QUAL<${qual_treshold} || \
        %MAX(FORMAT/AD[0:1])/%MAX(FORMAT/DP)<${call_treshold} || INFO/DP<${min_depth}' ${vcf} \
        | "${baseDir}/conda/env-variantcalling/bin/python" "$vcf2table" - \
        --sample ${sampleName} > ${sampleName}_3G.txt
        """
}


// creates consensus fasta's
// variants that are marked as uncertain or low covered regions are masked by N
process '4A_create_consensus' {
    tag '4A'
    conda "${baseDir}/conda/env-variantcalling/"
    publishDir outDir + '/report', mode: 'copy'
    input:
        //file vcf from rawvcf_4A
        file vcf from compressed_vcf_4A
        file vcf_index from compressed_vcf_index_4A
        file ubiq from ubiq_4A
        file noncov from noncov_4A
    output:
        file("${sampleName}.fasta") into consensus_6
        file ".command.*"
  script:
        """
		cat ${reference} | bcftools consensus ${vcf} \
		-m <(cat ${noncov} ${ubiq}) \
		${vcf} > ${sampleName}.fasta
        sed -i 's/^>NC_045512.2/>${sampleName}/g' ${sampleName}.fasta
        """
}

// annotation of genome
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
        snpEff ann -v NC_045512.2 ${vcf} -ud 0 -strict \
        -c ${baseDir}/db/snpEff.config -csvStats ${sampleName}_snpEff.csv \
        > ${vcf.baseName}_annot.vcf
        ${baseDir}/conda/env-variantcalling/bin/python $vcf2table ${vcf.baseName}_annot.vcf --sample ${sampleName} \
        -ad -e -o ${sampleName}_annot_table.txt
        """
}

// pangolin
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
        pangolin ${consensus}
        """
}

// multiqc
process '7_QC' {
    tag '7'
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

// Process 8: generate a report for interpretation by the clinician (or for research purposes)
process '8_report' {
    tag '8'
    conda 'conda-forge::jinja2=2.11.1 conda-forge::weasyprint=51 \
            conda-forge::simplejson=3.17.0 conda-forge::matplotlib=3.1.2 \
            conda-forge::pandas=1.0.1 conda-forge::open-fonts=0.7.0'
    publishDir outDir + '/report', mode: 'copy'
    input:
        file lineage from lineage_8
        file annotation from annotation_8
    output:
        file "${sampleName}.html"
        file "${sampleName}.pdf"
        file ".command.*"
    script:
        """
        $reporter --sampleName ${sampleName} --lineage ${lineage} \
        --annotation ${annotation}
        """
}
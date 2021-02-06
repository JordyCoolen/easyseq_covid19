# EasySeq RC-PCR SARS-CoV-2/COVID-19
# Variant pipeline V0.3

## Table of contents
* [GENERAL-INFO](#GENERAL-INFO)
* [INSTALL](#INSTALL)
* [RUN](#RUN)
* [OUTPUT](#OUTPUT)
* [FLOW-DIAGRAM](#FLOW-DIAGRAM)
* [TOOLS](#TOOLS)
* [DOCKER](#DOCKER)
* [CONTRIBUTORS](#CONTRIBUTORS)
* [REMARKS](#REMARKS)
* [REFERENCE](#REFERENCE)
* [DISCLAIMER](#DISCLAIMER)
* [LICENSE](*LICENSE)

## General info
This github repository contains an automated pipeline
dedicated to properly analyse the EasySeq SARS-CoV-2 (COVID-19) sequence
sequencing data. All validation are done using 149 bp or 151 bp paired-end reads.

In short:
* Automated pipeline to analyse Illumina EasySeq COVID-19 samples to a variant report
* The pipeline cleans the Illumina sequencing data
* Uses the SARS-CoV-2 reference genome (NC_045512.2)
* Mutations and deletions are measured
* Fasta consensus of the sample is created
* Lineage is determined
* Output is available in a structure way
* Full QC reports are created
* PDF and HTML report as output

## INSTALL
1. install docker on your OS
2. docker pull jonovox/easyseq_covid19:latest
3. download the newest release of the pipeline via 
   https://github.com/JordyCoolen/easyseq_covid19/releases
4. extract the source code
5. go into the extracted/project folder
6. download conda environments via:
   https://surfdrive.surf.nl/files/index.php/s/Nu67ByV5P6w3wqI
7. extract conda.tar.gz into the project folder created at step 5
8. Proceed to RUN examples

## RUN
### RUN_option1
- now you have to perform the test to set everything in place
- first time running the variant pipeline will deploy more conda environments needed 
  to successfully install the pipeline. This can take a while.
- open docker runtime container from image with write rights
```bash
sh docker/run.sh covid easyseq_covid19:latest
```

- run the test sample inside the container
```bash
nextflow run COVID.nf --sampleName test -resume --outDir /workflow/output/test --reads "/workflow/input/test_OUT01_R{1,2}.fastq.gz"
```
### RUN_option2
- you can also execute multiple samples in non-parallel way
```bash
bash docker/run.sh <path to folders containing the fastq.gz file> easyseq_covid19:latest
```

## OUTPUT
```bash
/workflow/output/test
|-- HV69-70
|   |-- test.aln
|   |-- test.frag.gz
|   |-- test.fsa
|   |-- test.res
|   `-- test_HVdel.vcf
|-- QC
|   |-- multiqc_data
|   |   |-- multiqc.log
|   |   |-- multiqc_data.json
|   |   |-- multiqc_fastp.txt
|   |   |-- multiqc_general_stats.txt
|   |   |-- multiqc_snpeff.txt
|   |   `-- multiqc_sources.txt
|   |-- multiqc_report.html
|   |-- test.fastp.json
|   |-- test.mosdepth.global.dist.txt
|   |-- test.mosdepth.summary.txt
|   |-- test.per-base.bed.gz
|   |-- test.per-base.bed.gz.csi
|   `-- test_snpEff.csv
|-- annotation
|   |-- snpEff_summary.html
|   |-- test_annot_table.txt
|   |-- test_snpEff.csv
|   `-- test_snpEff.genes.txt
|-- lineage
|   `-- lineage_report.csv
|-- mapping
|   |-- test.bam
|   |-- test.bam.bai
|   |-- test.primerclipped.bam
|   `-- test.primerclipped.bam.bai
|-- rawvcf
|   `-- test.vcf.gz
|-- report
|   |-- test.fasta
|   |-- test.html
|   `-- test.pdf
|-- uncovered
|   |-- test_noncov.bed
|   `-- test_ubiq.bed
|-- vcf
|   |-- notpassed
|   |   `-- test_3G.txt
|   |-- test.vcf
|   |-- test_final.vcf
|   `-- test_table.txt
`-- vcf_index
    |-- test_concat.vcf.gz
    `-- test_concat.vcf.gz.csi
```

## FLOW-DIAGRAM
![Alt text](flowchart.png?raw=true "Flowdiagram")

## TOOLS
* nextflow
* python
* conda/bioconda
* fastp
* BWA MEM
* samtools
* bcftools
* mosdepth
* bedtools
* snpEff
* KMA
* multiQC
* pangolin v2.1.11

## DOCKER
### build your own docker image
```bash
cd easyseq_covid19
docker build --rm -t <image name> ./
```

## CONTRIBUTORS
Department of Medical Microbiology and Radboudumc Center for Infectious Diseases, Radboud university medical center, Nijmegen, The Netherlands
* J.P.M. Coolen
  (jordy.coolen@radboudumc.nl)

NimaGen B.V., Nijmegen, The Netherlands
* R.A. Lammerts (NimaGen B.V., Nijmegen, The Netherlands)
* J.T. Vonk (Student HAN Bioinformatics, Nijmegen, The Netherlands)

## REMARKS
```bash
spike S
21765-21770HV 69-70 deletion

The current EasySeq design is not completly overlapping the region 21765-21770 / HV 69-70.

---->         To solve this a template based strategy using KMA is applied.
      <---    This method measures which template matches best. Either Wildtype (NC_045512.2) or
              a variant containing the 21765-21770 / HV 69-70 deletion. The result of this strategy
              is projected in the VCF to ensure correct output. This works perfect for now because no other deletions are
              known on this exact location.
```

## REFERENCE
For citing this work please cite:

**Novel SARS-CoV-2 Whole-genome sequencing technique using Reverse Complement PCR enables easy, fast and accurate outbreak analysis in hospital and community settings**
*Femke Wolters, Jordy P.M. Coolen, Alma Tostmann, Lenneke F.J. van Groningen, Chantal P. Bleeker-Rovers, Edward C.T.H. Tan, Nannet van der Geest-Blankert, Jeannine L.A. Hautvast, Joost Hopman, Heiman F.L. Wertheim, Janette C. Rahamat-Langendoen, Marko Storch, Willem J.G. Melchers
bioRxiv 2020.10.29.360578; doi: https://doi.org/10.1101/2020.10.29.360578.*

Also cite the other programs used, see list of used tools 

The work is currently under revision.

## DISCLAIMER
This is for Research Only.
The code and pipeline is continuously under development.
We cannot guarantee a full error free result.
Especially with the fast developments in SARS-CoV-2/COVID-19 sequencing and
the continuously mutating nature of the virus.

## LICENSE

[![CC0](https://licensebuttons.net/p/zero/1.0/88x31.png)](https://creativecommons.org/publicdomain/zero/1.0/)
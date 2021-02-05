# EasySeq RC-PCR SARS-CoV-2 (COVID-19) 
# Variant pipeline V0.3

## Table of contents
* [General info](#general-info)
* [INSTALL](#INSTALL)
* [Docker](#DOCKER)
* [Test](#Test)
* [Output](#output)
* [Flow diagram](#Flow diagram)
* [CONTRIBUTORS](#CONTRIBUTORS)
* [REFERENCE](#REFERENCE)
* [DISCLAIMER](#DISCLAIMER)
* [LICENSE](*LICENSE)

## General info
This github repository contains an automated pipeline
dedicated to properly analyse the EasySeq SARS-CoV-2 (COVID-19) sequence
sequencing data. In short: 
* The pipeline cleans the Illumina sequencing data
* Uses the SARS-CoV-2 reference genome (NC_045512.2)
* Mutations and deletions are measured
* Fasta consensus of the sample is created
* Lineage is determined
* Output is available in a structure way
* PDF and HTML rapport as output

## INSTALL

install docker on your system
download docker image
pull this repo into your system
cd to repo
open the image 
```bash
sh docker/run.sh covid covid:0.3
```

## DOCKER

### build docker image
```bash
docker build --rm -t covid:0.3 ./
```

### open docker runtime container from image with write rights
```bash
sh docker/run.sh covid covid:0.3
```

## Test

### test run inside the container
```bash
nextflow run COVID.nf --sampleName test -resume --outDir /workflow/output/test --reads "/workflow/input/test_OUT01_R{1,2}.fastq.gz"
```

### directly execute pipeline from outside the container
```bash
docker run -it --rm --mount type=bind,source=${PWD},target=/workflow covid:0.3 \
nextflow run COVID.nf --sampleName test -resume --outDir /workflow/output/test \
--reads "/workflow/input/test_OUT01_R{1,2}.fastq.gz"
```

## Output
```bash
/workflow/output/test/
|-- HV69-70
|   |-- test.aln
|   |-- test.frag.gz
|   |-- test.fsa
|   `-- test.res
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
|   |-- test.vcf.gz
|   `-- test.vcf.gz.csi
|-- report
|   |-- test.fasta
|   |-- test.html
|   `-- test.pdf
|-- uncovered
|   |-- test_noncov.bed
|   `-- test_ubiq.bed
`-- vcf
    |-- notpassed
    |-- test.vcf
    `-- test_3B.txt
```

## Flow diagram
![Alt text](flowchart.png?raw=true "Flowdiagram")

## List of used tools
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
* pangolin

## CONTRIBUTORS
Radboudumc
Department of Medical Microbiology and Radboudumc Center for Infectious Diseases, Radboud university medical center, Nijmegen, The Netherlands
* J.P.M. Coolen
jordy.coolen@radboudumc.nl

NimaGen B.V., Nijmegen, The Netherlands
* R. A. Lammerts (NimaGen, Nijmegen, The Netherlands)
* J.T. Vonk (Student HAN Bioinformatics, Nijmegen, The Netherlands)

## REFERENCE
For citing this work please cite:

Novel SARS-CoV-2 Whole-genome sequencing technique using Reverse Complement PCR enables easy, fast and accurate outbreak analysis in hospital and community settings
Femke Wolters, Jordy P.M. Coolen, Alma Tostmann, Lenneke F.J. van Groningen, Chantal P. Bleeker-Rovers, Edward C.T.H. Tan, Nannet van der Geest-Blankert, Jeannine L.A. Hautvast, Joost Hopman, Heiman F.L. Wertheim, Janette C. Rahamat-Langendoen, Marko Storch, Willem J.G. Melchers
bioRxiv 2020.10.29.360578; doi: https://doi.org/10.1101/2020.10.29.360578.

The work is currently under revision.

## DISCLAIMER
The code and pipeline is continuously under development.
We cannot garantee a full error free result.
Especially with the fast developments in SARS-CoV-2/COVID-19 sequencing and
the continuously mutating nature of the virus.

## LICENSE

[![CC0](https://licensebuttons.net/p/zero/1.0/88x31.png)](https://creativecommons.org/publicdomain/zero/1.0/)
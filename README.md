# EasySeq RC-PCR SARS-CoV-2/COVID-19 WGS kit
# Variant pipeline V0.7.0 
## Use with V3 of the WGS kit, else see general info

## Table of contents
* [GENERAL-INFO](#GENERAL-INFO)
* [INSTALL](#INSTALL)
* [RUN](#RUN)
* [OUTPUT](#OUTPUT)
* [FLOW-DIAGRAM](#FLOW-DIAGRAM)
* [TOOLS](#TOOLS)
* [PANGOLIN](#PANGOLIN)
* [DOCKER](#DOCKER)
* [SINGULARITY](#SINGULARITY)
* [CONTRIBUTORS](#CONTRIBUTORS)
* [REMARKS](#REMARKS)
* [REFERENCE](#REFERENCE)
* [DISCLAIMER](#DISCLAIMER)
* [LICENSE](*LICENSE)

## General info
This github repository contains an automated pipeline
dedicated to properly analyse the EasySeq SARS-CoV-2 (COVID-19) sequence
sequencing data. Validated with 150/151 bp paired-end reads.

Advice is to redownload the conda.tar.gz after each update to
be sure that all conda environments are set in place.

version 1 or 2 of the EasySeq RC-PCR SARS-CoV-2 WGS kit
* Use version v0.5.2 of the github
https://github.com/JordyCoolen/easyseq_covid19/releases/tag/v0.5.2

version 3 of the EasySeq RC-PCR SARS-CoV-2 WGS kit
* Use code version v0.7.0 and newer
* Implemented lofreq for variant calling which gives much more
accurate calls in the report. Consensus output is mostly unaffected.

In short:
* Automated pipeline to analyse Illumina EasySeq COVID-19 samples to a variant report
* The pipeline cleans the Illumina sequencing data
* Uses the SARS-CoV-2 reference genome (NC_045512.2)
* Custom EasySeq Primer filtering and correction
* Mutations and deletions are measured
* Fasta consensus of the sample is created
* Lineage is determined
* Output is available in a structured way
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
   https://surfdrive.surf.nl/files/index.php/s/ggoLXzMoa5iSZYa
7. extract conda.tar.gz into the project folder created at step 5
8. Proceed to RUN examples

## RUN
### RUN_option1
- now you have to perform the test to set everything in place
- first time running the variant pipeline will deploy more conda environments needed 
  to successfully install the pipeline. This can take a while.
- open docker runtime container from image with write rights
```bash
sh docker/run.sh covid jonovox/easyseq_covid19:latest
```

- run the test sample inside the container
```bash
nextflow run COVID.nf --sampleName test -resume --outDir /workflow/output/test --reads "/workflow/input/test_OUT01_R{1,2}.fastq.gz"
```
### RUN_option2
- you can also execute multiple samples in non-parallel way
```bash
bash scripts/run_batch.sh <path to folders containing the fastq.gz file> <extension of files> <threads> jonovox/easyseq_covid19:latest
```

## OUTPUT
```bash
/workflow/output/test/
├── QC
│   ├── multiqc_data
│   │   ├── multiqc.log
│   │   ├── multiqc_data.json
│   │   ├── multiqc_fastp.txt
│   │   ├── multiqc_general_stats.txt
│   │   ├── multiqc_snpeff.txt
│   │   └── multiqc_sources.txt
│   ├── multiqc_report.html
│   ├── stats.txt
│   ├── test.fastp.json
│   ├── test.mosdepth.global.dist.txt
│   ├── test.mosdepth.summary.txt
│   ├── test.per-base.bed.gz
│   ├── test.per-base.bed.gz.csi
│   └── test_snpEff.csv
├── annotation
│   ├── snpEff_summary.html
│   ├── test_annot_table.txt
│   ├── test_snpEff.csv
│   └── test_snpEff.genes.txt
├── lineage
│   └── lineage_report.csv
├── mapping
│   ├── test.bam
│   ├── test.bam.bai
│   ├── test.final.bam
│   └── test.final.bam.bai
├── rawvcf
│   └── test.raw.vcf
├── report
│   ├── parameters.txt
│   ├── test.fasta
│   ├── test.html
│   └── test.pdf
├── uncovered
│   ├── test_noncov.bed
│   └── test_ubiq.bed
└── vcf
    ├── notpassed
    │   └── test.notpassed.vcf
    ├── test.final.vcf
    ├── test.final.vcf.gz
    ├── test.final.vcf.gz.csi
    └── test.variants.vcf
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
* lofreq
* mosdepth
* bedtools
* snpEff
* KMA
* multiQC
* pangolin v2.3.8 (pangoLEARN 2021-04-01) (default in conda.tar.gz)

## PANGOLIN
### to update the pangolin tool and database perform following commands
```bash
sh docker/run.sh covid jonovox/easyseq_covid19:latest
conda activate /workflow/conda/env-pangolin
pangolin --update
```

## DOCKER
### build your own docker image
```bash
cd easyseq_covid19
docker build --rm -t <image name> ./
```

## SINGULARITY
## build SINGULARITY IMAGE from dockerhub
```bash
singularity build <imagename>.simg docker://jonovox/easyseq_covid19:latest
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
21765-21770 HV 69-70 deletion

Version 1 and 2 of the EasySeq RC-PCR SARS-CoV-2 WGS kit are not completly overlapping the region 21765-21770 / HV 69-70.
If you use these versions of the WGS kit please use:

variant pipeline v0.5.2
https://github.com/JordyCoolen/easyseq_covid19/releases/tag/v0.5.2

---->         This version solves the not overlapping region of 21765-21770 by using a template based strategy using KMA.
      <---    This method measures which template matches best. Either Wildtype (NC_045512.2) or
              a variant containing the 21765-21770 / HV 69-70 deletion. The result of this strategy
              is projected in the VCF to ensure correct output. This works perfect for now because no other deletions are
              known on this exact location.

variant pipeline v0.7.0

  ---->     In Version 3 of the EasySeq RC-PCR SARS-CoV-2 WGS kit the region 21765-21770 / HV 69-70 region is           
    <----   complety overlapping by having a new primer design. This version of the variant pipeline handles the 
            data obtained using version 3 correctly.
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
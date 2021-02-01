# NEXTFLOW EasySeq RC-PCR COVID19 V0.1

## build docker image
```bash
docker build --rm -t covid:0.2 ./
```

## open docker runtime container from image with write rights
```bash
sh docker/run.sh covid covid:0.2
```

## directly execute pipeline
```bash
docker run -it --rm --mount type=bind,source=${PWD},target=/workflow covid:0.2 \
nextflow run COVID.nf --sampleName test -resume --outDir output/test \
--reads "test/test_OUT01_R{1,2}.fastq.gz"
```

## test run
```bash
nextflow run COVID.nf --sampleName test -resume --outDir output/test --reads "test/test_OUT01_R{1,2}.fastq.gz"
```

## output of test
```bash
`-- test
    |-- annotation
    |   |-- snpeff_genes.txt
    |   |-- snpeff_summary.html
    |   `-- test_annot.vcf
    |-- fastp
    |   |-- test_OUT01_R1.fastq_fastp.fastq.gz
    |   `-- test_OUT01_R2.fastq_fastp.fastq.gz
    |-- json
    |-- lineage
    |   `-- lineage_report.csv
    |-- mapping
    |   `-- test.bam
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

## Flow diagram of current pipeline:

![Alt text](flowchart.png?raw=true "Flowdiagram")
#!/bin/bash

# USAGE
# bash run_batch.sh <inputpath>

for fname in ${1}/*_R1.fastq.gz
do
    base=${fname##*/}
    base=${base%_R1*}
    echo "${base}_R1.fastq.gz"
    echo "${base}_R2.fastq.gz"

    docker run -it --rm --mount type=bind,source=${PWD},target=/workflow \
    --mount type=bind,source=${1},target=/workflow/test \
    covid:0.2 nextflow run COVID.nf --sampleName ${base} \
    --outDir /workflow/output/${base} --reads "/workflow/test/${base}_R{1,2}.fastq.gz" -resume

done
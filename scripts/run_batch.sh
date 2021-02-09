#!/bin/bash

# USAGE
# bash run_batch.sh <inputpath> <file_extension> <image>
# bash run_batch.sh /workflow/input/ _001.fastq.gz jonovox/easyseq_covid19:latest
# <file_extension> most common _001.fastq.gz

for fname in ${1}/*_R1${2}
do
    base=${fname##*/}
    base=${base%_R1*}
    echo "${base}_R1${2}"
    echo "${base}_R2${2}"

    docker run -it --rm --mount type=bind,source=${PWD},target=/workflow \
    --mount type=bind,source=${1},target=/workflow/input \
    ${3} nextflow run COVID.nf --sampleName ${base} \
    --outDir /workflow/output/${base} --reads "/workflow/input/${base}_R{1,2}${2}" -resume

done
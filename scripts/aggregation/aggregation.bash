# bash aggregation.bash <inputpath> <outputpath>

# create folders
mkdir ${2}/summary
mkdir ${2}/summary/pdf
mkdir ${2}/summary/fasta
mkdir ${2}/summary/annotation
mkdir ${2}/summary/bam

# obtain aggregated output files
cat ${1}/*/lineage/*.csv > ${2}/summary/lineages.txt
cat ${1}/*/QC/stats.txt > ${2}/summary/stats.txt
cp ${1}/*/report/*.fasta ${2}/summary/fasta/
cat ${2}/summary/fasta/*.fasta > ${2}/summary/all.fasta
cp ${1}/*/report/*.pdf ${2}/summary/pdf/
cp ${1}/*/mapping/*.final.bam* ${2}/summary/bam/
cp ${1}/*/annotation/*_annot_table.txt ${2}/summary/annotation
python ${PWD}/aggregate_mutations.py -i ${2}/summary/annotation -o ${2}/summary/annotation

# run script to merge all files to one output file
#!/bin/bash

fastq_dir=""
genome_path=""
csv_file=""
output_dir=""
localcores=16
localmem=100

# read csv, extract second column

while IFS= read -r line
do
    sample=$(echo $line | cut -d, -f2)
    # if starts with # or if named Barcode, skip
    if [[ $sample == \#* ]] || [[ $sample == "Barcode" ]]; then
        continue
    fi
    echo $sample

    ~/cellranger-9.0.0/bin/cellranger count --id=$sample \
        --transcriptome=${genome_path}/refdata-gex-GRCm39-2024-A \
        --fastqs=$fastq_dir \
        --sample=$sample \
        --localcores=$localcores \
        --localmem=$localmem \
        --create-bam true \
        --include-introns=true

    mv $sample $output_dir
done < $csv_file

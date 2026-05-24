#!/bin/bash
#../../cellranger_analysis/step0_merge_data/A*_R1_001.fastq.gz
pre='/Acropora_sc_analysis/cellranger_analysis/step0_merge_data/'
lat_a='_R1_001.fastq.gz'
lat_b='_R2_001.fastq.gz'
targets=("Acropora1" "Acropora2" "Acropora3")

for target in ${targets[@]}; do
    target_a=$pre$target$lat_a
    target_b=$pre$target$lat_b
    sbatch kallisto_index_generator.slurm $target $target_a $target_b
done 

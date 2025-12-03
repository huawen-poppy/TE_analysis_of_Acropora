#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:00:00
#SBATCH --mem=5000
#SBATCH -o fastqc.%J.out
#SBATCH -e fastqc.%J.err

# fastQC
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# below is for doing the fastqc (quality checking) for the raw data

set -vex

module load fastqc

outDir='output'
mkdir -p $outDir

threads=2

for f in ../step0_merge_data/*R2_001.fastq.gz; do
        echo "fastQC'ing file"${f}
        fastqc $f -o ${outDir} -t ${threads}
done

# specify R2 is due to R2 means cDNA sequences
# R1 means barcodes
# I1 means illumina lane info

#!/bin/bash --login
#SBATCH -N 10
#SBATCH --partition=batch
#SBATCH -J kb_index
#SBATCH -o /TE_analysis/acropora/logs/kb_bus.%J.out
#SBATCH -e /TE_analysis/acropora/logs/kb_bus.%J.err
#SBATCH --time=15:30:00
#SBATCH --mem=64G

#run the application:
conda activate bustools

cd /TE_analysis/acropora

kallisto bus -i acropora_cds_all_TEs.index -t 30 -o $1 -x 10xv3 $2 $3 



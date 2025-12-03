#!/bin/bash
#SBATCH -N 3
#SBATCH --cpus-per-task=32
#SBATCH -J cellranger_reanalyze
#SBATCH -o reanalyze.%J.out
#SBATCH -e reanalyze.%J.err
#SBATCH --time=168:00:00

# this step is also optional, it is use the cellranger's filter ceriteria to do the counting analysis on the acropora

#Load the required module
module load cellranger/5.0.1


#run the application:
cellranger reanalyze --id=Acropora_reanalyze --matrix=./Acropora_sc_analysis/cellranger_analysis/step4_aggregate/Acropora_aggr/outs/count/filtered_feature_bc_matrix.h5 --params=reanalyze_Acropora.csv


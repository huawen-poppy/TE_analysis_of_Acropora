#!/bin/bash
#SBATCH -N 3
#SBATCH --cpus-per-task=32
#SBATCH -J aggregate
#SBATCH -o aggr.%J.out
#SBATCH -e aggr.%J.err
#SBATCH --time=168:00:00


# this step is optional, just to aggregate all the samples together
#Load the required module
module load cellranger/5.0.1

#first generate the csv file from the .h5 files
#nano acropora_aggr.csv
#put below infor inside
#library_id,molecule_h5
#Acropora1,./Acropora_sc_analysis/cellranger_analysis/step3_count/run_count_Acropora1/outs/molecule_info.h5
#Acropora2,./Acropora_sc_analysis/cellranger_analysis/step3_count/run_count_Acropora2/outs/molecule_info.h5
#Acropora3,./Acropora_sc_analysis/cellranger_analysis/step3_count/run_count_Acropora3/outs/molecule_info.h5

#run the application:
cellranger aggr --id=Acropora_aggr --csv=acropora_aggr.csv

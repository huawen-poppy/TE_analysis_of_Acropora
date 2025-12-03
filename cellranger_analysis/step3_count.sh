#!/bin/bash
#SBATCH -N 3
#SBATCH --cpus-per-task=32
#SBATCH -J cellranger_count
#SBATCH -o count.%J.out
#SBATCH -e count.%J.err
#SBATCH --time=168:00:00


# we are counting the expression of each gene/transcript in this step

#in this step, the cellranger count require the sample named as Sample_S1_L00X_R1/R2/I1_001.fastq
#and also require two lanes data(more than two lanes is still ok)
#also, this step can only deal with one sample(experiment in a well) at one time
#if the data comes from the same well with multisamples, then we should use the cellranger multi 

#Load the required module
module load cellranger/5.0.1


#run the application:
sample=(Acropora1 Acropora2 Acropora3)

for element in ${sample[@]};
do
if [ "${element:0:1}" == "A" ];then
cellranger count --id="run_count_BMCGen2018_algae_"${element} \
--fastqs=./Acropora_sc_analysis/cellranger_analysis/step3_count/input_data \
--sample=${element} \
--transcriptome=./Acropora_sc_analysis/cellranger_analysis/step2_make_reference/acropora_genome

fi
done
 

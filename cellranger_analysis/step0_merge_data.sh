#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=5
#SBATCH -J merge_lane_data
#SBATCH -o merge.%J.out
#SBATCH -e merge.%J.err
#SBATCH --mail-type=ALL
#SBATCH --time=01:30:00
#SBATCH --mem=50G

#this job is to merge the raw data from lane1 and lane2 into one

# if I do it manuelly, then I have to repeat the following command several times by changing the group name(S1,S2,S3,A1,A2,A3) and index name(I1,R1,R2)
# for example:
# cat ../AipSC/Lane1/version_01/M-19-4777_S1_SI-GA-A1_S1_L001_I1_001.fastq.gz ../AipSC/Lane2/version_01/M-19-4777_S1_SI-GA-A1_S1_L002_I1_001.fastq.gz > ./AipSC_analysis/cellranger_analysis/step0_merge_data/merged_S1_I1.fastq.gz


set -vex
for file_l1 in ./Acropora_sc_analysis/raw_data/Lane1/version_01/*
do
#echo ${file_l1} >> filename.txt
file_l1_path=${file_l1}
#echo ${file_l1_path} >> file_l1_path.txt
file_l2_path1=${file_l1_path/L001/L002}  #replace L001 with L002
file_l2_path=${file_l2_path1/Lane1/Lane2}  #replace Lane1 with Lane2
#echo ${file_l2_path} >> file_l2_path.txt
if [ -f ${file_l2_path} ];then
        merged_file_name_left1=${file_l1_path##*/}  #remove from file_l1_path everything before the last "/", means get the text behind the last "/"
	merged_file_name_left=${merged_file_name_left1:10:9}   #get the text start from position 11, text length be 2, here we get the species label
	#echo ${merged_file_name_left} >> yes.txt
        merged_file_name_right=${file_l1_path#*L001_}  #remove the text infront of L001_, including L001_
	#echo ${merged_file_name_right} >> yes.txt
        merge_file_name=${merged_file_name_left}"_"${merged_file_name_right}
	echo ${merge_file_name} >> yes.txt
        cat ${file_l1_path} ${file_l2_path} > ${merge_file_name}
fi
done
 

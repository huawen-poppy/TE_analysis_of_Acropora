#!/bin/bash --login
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH -J mkref
#SBATCH -o mkref.%J.out
#SBATCH -e mkref.%J.err
#SBATCH --time=24:00:00
#SBATCH --mem=100G

# in this step is to prepare the reference for the sequence mapping

#Load the required module
module load cellranger
cd ./Acropora_sc_analysis/cellranger_analysis/step2_make_reference
#convert the gff3 file to gtf file using the tool cufflink
#eg:
#gffread 4.321.CC7v2_braker.pasaupdated.tidy.gff3 -T -o CC7.gtf

#run the application:
#the step here is to filter the GTF,to remove the non-polyA transcripts(multi-mapped genes)

#cellranger mkgtf \
#Ahem_transcript_rename_genename2.gtf \
#Ahem_transcript_rename_genename2_filtered.gtf \
#--attribute=gene_biotype:protein_coding

cellranger mkref \
--genome=acropora_genome \
--fasta=Acropora_genome.fa \
--genes=Acropora_genome_rename2.gtf 


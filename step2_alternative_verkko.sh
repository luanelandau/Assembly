#!/bin/bash
#SBATCH --account=omergokc
#SBATCH --partition=general-compute
#SBATCH --qos=nih
#SBATCH --time=40:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=900G
#SBATCH --job-name="verkko-1946"
#SBATCH --output=verkko_aug14.out
#SBATCH --error=verkko_aug14.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue

#Verkko is one of the best assemblers for human genomes. The problem is that it requires a lot of read depth, hifi & ONT reads. With the new corrector herro
#it is possible to use verkko only with ONT long reads after correction. This script was adapted from the paper: https://doi.org/10.1101/2024.05.18.594796

#this script is still being tested!! (Aug 15, 2024)

module load gcc/11.2.0  openmpi/4.1.1 snakemake/7.18.2 r/4.2.0
module load samtools seqtk/1.3 
eval "$(/projects/academic/omergokc/charikleia/anaconda3/bin/conda shell.bash hook)"
conda activate herro

#For this version of the assembler, I am using trios. Therefore, I had to download the illumina short reads for each of the parents and for the child in order to phase the assembly. The next steps
#are just me concatenating all the illumina files into one for each parent and one for the child.
#Only do this step if you do not have it already done.
#cat /vscratch/grp-omergokc/Kendra/1k_Genomes/HG01944/*.fastq.gz > HG01944.fastq.gz #You have to concatenate all the files into one file for each parent and one file for child.
#cat /vscratch/grp-omergokc/Kendra/1k_Genomes/HG01945/*.fastq.gz > HG01945.fastq.gz
#cat /vscratch/grp-omergokc/Kendra/1k_Genomes/HG01946/*.fastq.gz > HG01946_SR.fastq.gz #SR as in short reads

#One has to use the fastq or fasta file for the analysis, the zipped file wont work. 

#gzip /projects/academic/omergokc/Luane/HG01946/HG01946_allreads_Apr16.fastq
#cat /projects/academic/omergokc/Luane/HG01946/HG01944/*.fastq.gz > HG01944.fastq.gz
#cat /projects/academic/omergokc/Luane/HG01946/HG01945/*.fastq.gz > HG01945.fastq.gz

##---------HAPLOID---------------
#verkko -d HG01946_verkko --hifi HG01946_10kbp_HERRO.fasta --nano HG01946_allreads_Apr16.fastq --local-memory 900 --local-cpus 16 --cleanup --haploid
#Dont forget to match the local memory with what you asked for in the script. This has given me problems before.

##----------DIPLID----------------
#Run the following for the assembly sample and the parents’ short reads, this step was done with Meryl 1.4.1 for all samples except for I002C which used Meryl 1.4
#child

meryl count compress k=30 threads=32 memory=900 /projects/academic/omergokc/Luane/HG01946/HG01946_allreads_Apr16.fastq output child_compressed.k30.meryl

#AUG 12, I am trying to use short reads for this step. I am not sure this is the way to do it.
meryl count compress k=30 threads=32 memory=600 HG01946_SR.fastq.gz output child_compressed.k30.meryl

#mother
meryl count compress k=30 threads=32 memory=900 HG01944.fastq.gz output maternal_compressed.k30.meryl

#father
meryl count compress k=30 threads=32 memory=900 HG01945.fastq.gz output paternal_compressed.k30.meryl

#Create hapmers for verkko assembly, Merqury commit 4b4846a78e8eaa06e42e06ec9f40c7f329e664c4 (https://github.com/marbl/merqury/tree/4b4846a78e8eaa06e42e06ec9f40c7f329e664c4) with Meryl 1.4.1 was used for all samples for this step

/projects/academic/omergokc/Luane/softwares/merqury/trio/hapmers.sh ../maternal_compressed.k30.meryl ../paternal_compressed.k30.meryl ../child_compressed.k30.meryl

#c.	For diploid inputs (with parental data):

#verkko -d HG01946_verkko --hifi ../herro/HG01946_10kbp_HERRO.fasta --nano ../HG01946_allreads_Apr16.fastq --local-memory 900 --local-cpus 32 --snakeopts "--unlock" --cleanup --hap-kmers maternal_compressed.k30.hapmer.meryl paternal_compressed.k30.hapmer.meryl trio

verkko -d HG01946_verkko --hifi /projects/academic/omergokc/Luane/HG01946/herro/HG01946_10kbp_HERRO.fasta --nano ../../HG01946_allreads_Apr16.fastq --local-memory 900 --local-cpus 32 --cleanup --hap-kmers maternal_compressed.k30.hapmer.meryl paternal_compressed.k30.hapmer.meryl trio

conda deactivate

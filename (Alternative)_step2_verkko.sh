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

#Verkko (https://github.com/marbl/verkko) is one of the best assemblers for human genomes. The problem is that it requires a lot of read depth, hifi & ONT reads. 
#With the new corrector herro it is possible to use verkko only with ONT long reads after correction. 
#This script was adapted from the paper: https://doi.org/10.1101/2024.05.18.594796

#this script is still being tested!! (Sept 10, 2024)

module load gcc/11.2.0  openmpi/4.1.1 snakemake/7.18.2 r/4.2.0
module load samtools seqtk/1.3 
eval "$(/projects/academic/omergokc/charikleia/anaconda3/bin/conda shell.bash hook)"
conda activate herro #I unfortunately have meryl and verkko in two different environments. But you can have all of them in one. This is why I am activating herro and further on I activate verkko.

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

#Defining inputs. Note that this is a generic script and I am not using the same individuals here and before. But the idea is that you define the files for the two parents and for the child. These are short reads.
par1="/projects/academic/omergokc/Luane/HG01972/verkko/HG01970.fastq.gz"
par2="/projects/academic/omergokc/Luane/HG01972/verkko/HG01971.fastq.gz"
child="/projects/academic/omergokc/Luane/HG01972/verkko/HG01972_SR.fastq.gz"

#The reads are the HERRO corrected reads, that will be inputed as pac bio hifi data. 
#The raw are the raw nanopore ultra long reads.
reads="/projects/academic/omergokc/Luane/HG01972/herro/HG01972_10kbp_HERRO.fasta"
raw="/projects/academic/omergokc/Luane/HG01972/herro/HG01972_allreads.fastq"


##---------HAPLOID---------------
#Run the following if you dont have trios and cant phase your genome. 
#verkko -d HG01946_verkko --hifi $reads --nano $raw --local-memory 900 --local-cpus 32 --cleanup --haploid
#Dont forget to match the local memory with what you asked for in the script. This has given me problems before.

##----------DIPLID----------------
#Run the following step-by-step if you have the trios. This will phase the genome and generate a phased assembly based on the kmers and short reads of the parents. 
#Run the following for the assembly sample and the parents’ short reads, this step was done with Meryl 1.4.1 for all samples except for I002C which used Meryl 1.4
#child

#child - make sure you're using the short reads here
meryl count compress k=30 threads=32 memory=900 $child output child_compressed.k30.meryl #adjust memory according to script

#mother
meryl count compress k=30 threads=32 memory=900 $par1 output maternal_compressed.k30.meryl #adjust memory according to script

#father
meryl count compress k=30 threads=32 memory=900 $par2 output paternal_compressed.k30.meryl #adjust memory according to script

#Create hapmers for verkko assembly, Merqury commit 4b4846a78e8eaa06e42e06ec9f40c7f329e664c4 (https://github.com/marbl/merqury/tree/4b4846a78e8eaa06e42e06ec9f40c7f329e664c4) with Meryl 1.4.1 was used for all samples for this step
/projects/academic/omergokc/Luane/softwares/merqury/trio/hapmers.sh ../maternal_compressed.k30.meryl ../paternal_compressed.k30.meryl ../child_compressed.k30.meryl

#define the parents variables again
par1_hap="/projects/academic/omergokc/Luane/HG01972/verkko/diploid/maternal_compressed.k30.hapmer.meryl"
par2_hap="/projects/academic/omergokc/Luane/HG01972/verkko/diploid/paternal_compressed.k30.hapmer.meryl"

#c.	For diploid inputs (with parental data):
conda activate verkko #make sure to have verkko 2.2 installed. I was using an old version that gave me problems.

#you only run this if you need to resume a run that didnt finish. 
#verkko -d HG01972_verkko --hifi ${reads} --nano ${raw} --local-memory 900 --local-cpus 32 --snakeopts "--unlock" --hap-kmers $par1_hap $par2_hap trio 
#This step will unlock the folder using snakemake and then the next step will be your actual run. I do not know why, but this was the only way I found to resume runs.

verkko -d HG01972_verkko --hifi ${reads} --nano ${raw} --local-memory 900 --local-cpus 32 --hap-kmers $par1_hap $par2_hap trio #adjust memory according to script

conda deactivate


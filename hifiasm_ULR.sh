#!/bin/bash
#SBATCH --qos=nih
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=600G
#SBATCH --job-name=1922_hifiasm
#SBATCH --output=1922_hifiasm.out
#SBATCH --error=1922_hifiasm.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue

module load gcccore/11.2.0
module load hifiasm/0.19.6
module load seqkit/2.7.0

cor_reads="/projects/academic/omergokc/Luane/HG01922/herro/HG01922_10kbp_HERRO.fasta"
#par1="/projects/academic/omergokc/Luane/HG01922/verkko/diploid/HG01920.fastq.gz"
#par2="/projects/academic/omergokc/Luane/HG01922/verkko/diploid/HG01921.fastq.gz"
raw="/projects/academic/omergokc/Luane/HG01922/HG01922_allreads.fastq.gz"

out_pre="HG01922"

#If coverage is higher than 38x, then the reads need to be downsampled
#Corrected reads were first chopped to a length of 30,000 bp, and chunks shorter than 10,000 bp were filtered away using SeqKit v2.5.1:
#seqkit sliding -s 30000 -W 30000 -g ${cor_reads} > chopped.fasta
#seqkit seq -m 10000 chopped.fasta > processed.fasta


#If coverage is higher than 38x, then the reads need to be downsampled
#Corrected reads were first chopped to a length of 30,000 bp, and chunks shorter than 10,000 bp were filtered away using SeqKit v2.5.1:
#seqkit sliding -s 30000 -W 30000 -g ${cor_reads} > chopped.fasta
#seqkit seq -m 10000 chopped.fasta > processed.fasta

#October 25
seqkit seq -m 50000 $raw > raw_reads_HG01922_50kb.fastq.gz

seqkit seq -m 100000 $raw > raw_reads_HG01922_100kb.fastq.gz

par1="/projects/academic/omergokc/Luane/HG01922/hifiasm/diploid/parent1.yak"
par2="/projects/academic/omergokc/Luane/HG01922/hifiasm/diploid/parent2.yak"


/projects/academic/omergokc/Luane/softwares/hifiasm/hifiasm -t 32 -o HG01946_diploid_50kb --ul raw_reads_HG01946_50kb.fastq --ul-cut 50000 --dual-scaf -1 $par1 -2 $par2 $cor_reads 

/projects/academic/omergokc/Luane/softwares/hifiasm/hifiasm -t 32 -o HG01946_diploid_100kb --ul raw_reads_HG01946_100kb.fastq --ul-cut 100000 --dual-scaf -1 $par1 -2 $par2 $cor_reads   

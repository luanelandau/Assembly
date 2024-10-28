#!/bin/bash
#SBATCH --qos=nih
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=600G
#SBATCH --job-name=1972_hifiasm
#SBATCH --output=1972_hifiasm.out
#SBATCH --error=1972_hifiasm.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue

module load gcccore/11.2.0
module load hifiasm/0.19.6
module load seqkit/2.7.0

cor_reads="/projects/academic/omergokc/Luane/HG01972/herro/HG01972_10kbp_HERRO.fasta"
#par1="/projects/academic/omergokc/Luane/HG01972/verkko/HG01971.fastq.gz"
#par2="/projects/academic/omergokc/Luane/HG01972/verkko/HG01970.fastq.gz"

out_pre="HG01972"

raw="/projects/academic/omergokc/Luane/HG01972/HG01972_allreads.fastq.gz"


#October 25
seqkit seq -m 50000 $raw > raw_reads_HG01972_50kb.fastq.gz

seqkit seq -m 100000 $raw > raw_reads_HG01972_100kb.fastq.gz

par1="/projects/academic/omergokc/Luane/HG01972/hifiasm/phased/parent1.yak"
par2="/projects/academic/omergokc/Luane/HG01972/hifiasm/phased/parent2.yak"


/projects/academic/omergokc/Luane/softwares/hifiasm/hifiasm -t 32 -o HG01946_diploid_50kb --ul raw_reads_HG01972_50kb.fastq.gz --ul-cut 50000 --dual-scaf -1 $par1 -2 $par2 $cor_reads 

/projects/academic/omergokc/Luane/softwares/hifiasm/hifiasm -t 32 -o HG01946_diploid_100kb --ul raw_reads_HG01972_100kb.fastq.gz --ul-cut 100000 --dual-scaf -1 $par1 -2 $par2 $cor_reads   

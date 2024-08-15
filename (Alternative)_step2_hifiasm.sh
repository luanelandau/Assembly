#!/bin/bash
#SBATCH --qos=nih
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=18:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=600G
#SBATCH --job-name=1946_hifiasm
#SBATCH --output=1946_hifiasm.out
#SBATCH --error=1946_hifiasm.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue

#hifiasm: https://github.com/chhylp123/hifiasm

module load gcccore/11.2.0
module load hifiasm/0.19.6
module load seqkit/2.7.0

cor_reads="/projects/academic/omergokc/Luane/HG01946/herro/HG01946_10kbp_HERRO.fasta" #Reads corrected using herro
par1="/projects/academic/omergokc/Luane/HG01946/verkko/HG01944.fastq.gz" #parent 1 short reads from illumina sequencing
par2="/projects/academic/omergokc/Luane/HG01946/verkko/HG01945.fastq.gz" #parent 2 short reads from illumina sequencing

#If coverage is higher than 38x, then the reads need to be downsampled
#Corrected reads were first chopped to a length of 30,000 bp, and chunks shorter than 10,000 bp were filtered away using SeqKit v2.5.1:
seqkit sliding -s 30000 -W 30000 -g ${cor_reads} > chopped.fasta
seqkit seq -m 10000 chopped.fasta > processed.fasta


#please see to install yak: https://github.com/lh3/yak
/projects/academic/omergokc/Luane/softwares/yak/yak count -b37 -t 32 -o parent1.yak <(zcat ${par1})  
/projects/academic/omergokc/Luane/softwares/yak/yak count -b37 -t 32 -o parent2.yak <(zcat ${par2})   


/projects/academic/omergokc/Luane/softwares/hifiasm/hifiasm -t 32 -o output_prefix --ul ${cor_reads} --ul-cut 50000 --dual-scaf -1 parent1.yak -2 parent2.yak processed.fasta   

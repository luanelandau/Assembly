#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=500G
#SBATCH --job-name=2252_hifiasm_haploid
#SBATCH --output=2252_hifiasm_haploid.out
#SBATCH --error=2252_hifiasm_haploid.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
##SBATCH --requeue

#hifiasm: https://github.com/chhylp123/hifiasm
#Following the pipeline published here: https://doi.org/10.1101/2024.05.18.594796

module load gcccore/11.2.0
module load hifiasm/0.19.6
module load seqkit/2.7.0

cor_reads="/projects/academic/omergokc/Luane/HG02252/herro/HG02252_10kb_HERRO.fasta" #reads corrected with fasta

#For Haploid samples
/projects/academic/omergokc/Luane/softwares/hifiasm/hifiasm -t 32 -o HG02006_hifiasm_haploid -l0 ${cor_reads}

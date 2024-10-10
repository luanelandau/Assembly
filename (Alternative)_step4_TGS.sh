#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=700G
#SBATCH --job-name="TGS2006"
#SBATCH --output=2006_TGS_hifiasm.out
#SBATCH --error=2006_TGS_hifiasm.err
#SBATCH --mail-user=luanejan@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --export=NONE

#TGS-GapCloser (https://github.com/BGI-Qingdao/TGS-GapCloser) is a program that uses the reads to fill the 
#gaps of your assemblies and create more contiguous assemblies. This is the last step of my assemblies.

#conda envirnoment I have it installed
eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate Assembly

#Inputs 
hap1="/projects/academic/omergokc/Luane/HG02006/ragtag/ragtag_output_hap1/ragtag.scaffold.fasta"
hap2="/projects/academic/omergokc/Luane/HG02006/ragtag/ragtag_output_hap2/ragtag.scaffold.fasta"
reads="/projects/academic/omergokc/Luane/HG01922/herro/HG01922_10kbp_HERRO.fasta" #your corrected reads
threads=32

#Haplotype 1
/projects/academic/omergokc/Luane/softwares/TGS-GapCloser/tgsgapcloser2  \
        --scaff  $hap1 \
        --reads  $reads \
        --output HG2006 \
        --ne \
        --thread $threads \

mkdir -p alternative
cd alternative

#Haplotype 2
/projects/academic/omergokc/Luane/softwares/TGS-GapCloser/tgsgapcloser2  \
        --scaff  $hap2 \
        --reads  $reads \
        --output HG2006_alt \
        --ne \
        --thread $threads \

conda deactivate

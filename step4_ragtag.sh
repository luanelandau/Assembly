#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=54
#SBATCH --mem=900G
#SBATCH --job-name=ragtag_HG01946
#SBATCH --output=ragtag_HG01946.out
#SBATCH --error=ragtag_HG01946.err
#SBATCH --constraint=AVX512
#SBATCH --export=NONE

#ragtag is a software (https://github.com/malonge/RagTag) for correcting and scaffolding a draft assembly. 

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate ragtag #this is my conda environment in which I have ragtag installed.

#input
assembly="/projects/academic/omergokc/Luane/HG01946/pilon/pilon4_out/HG01946_10kbp_flye_Medaka_4pilon.fasta" #path to draft genome
THRDS=40 #number threads
sample="HG01946" #sample name
Ref="/projects/academic/omergokc/Kendra/Ancients/CONGA/T2T.fa" #path to reference genome

# correct a query assembly
ragtag.py correct ${Ref} ${assembly}

# scaffold a query assembly
ragtag.py scaffold ${Ref} ragtag_output/ragtag.correct.fasta

# make joins and fill gaps in target.fa using sequences from query.fa
#ragtag.py patch ${Ref} ragtag_output/ragtag.scaffold.fasta #I generally do not trust this step, but it is doable. Please read the manual for deciding wether to use it or not.

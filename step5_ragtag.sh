#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=300G
#SBATCH --job-name=ragtag1946
#SBATCH --output=ragtag_HG01946.out
#SBATCH --error=ragtag_HG01946.err
#SBATCH --export=NONE

#ragtag is a software (https://github.com/malonge/RagTag) for correcting and scaffolding a draft assembly. 

#I have it installed in my conda environment
eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate ragtag

#input
phased_assembly="/projects/academic/omergokc/Luane/HG01946/verkko/diploid/HG01946_verkko_set5/assembly.haplotype1.fasta" #path to draft genome
phased_assembly_2="/projects/academic/omergokc/Luane/HG01946/verkko/diploid/HG01946_verkko_set5/assembly.haplotype2.fasta" #path to draft genome
THRDS=32 #number threads
sample="HG01946" #sample name
Ref="/projects/academic/omergokc/hg38.fa"

#RagTag is only compatible with fasta files. 
#Therefore, make sure your output from hifiasm is converted using these commands:
#gfa tools is installed in the herro environment
#gfatools gfa2fa $phased_assembly > HG01972.dip.hap1.p_ctg.fasta

# scaffold a query assembly
ragtag.py scaffold ${Ref} ${phased_assembly}

#ragtag doesnt have an option to define the output name, so I need to rename the output folder 
#before running it for the next haplotype.
mv ragtag_output/ ragtag_output_hap1

ragtag.py scaffold ${Ref} ${phased_assembly_2}

mv ragtag_output/ ragtag_output_hap2

conda deactivate

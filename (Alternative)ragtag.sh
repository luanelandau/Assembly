#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=300G
#SBATCH --job-name=ragtag2006
#SBATCH --output=ragtag_HG02006.out
#SBATCH --error=ragtag_HG02006.err
#SBATCH --export=NONE

#ragtag is a software (https://github.com/malonge/RagTag) for correcting and scaffolding a draft assembly. 
#I am not correcting the assembly here, because the herro program is already very accurate.

#I have it installed in my conda environment
eval "$(/projects/academic/omergokc/Luane/softwares/anaconda_new/bin/conda shell.bash hook)"
conda activate ragtag

#input
#in the case of an unphased assembly, hifiasm generates an "alternative" file (which here im leaving as phased_assembly_2) that contains unresolved contigs. I want to take a look at that.
phased_assembly="/projects/academic/omergokc/Luane/HG02006/hifiasm/with_readcorrection/HG02006_hifiasm_haploid.bp.p_ctg.fasta" #path to draft genome haplotype 1
phased_assembly_2="/projects/academic/omergokc/Luane/HG02006/hifiasm/with_readcorrection/HG02006_hifiasm_haploid.bp.a_ctg.fasta" #path to draft genome haplotype 2
THRDS=32 #number threads
sample="HG02006" #sample name
Ref="/projects/academic/omergokc/hg38.fa"

#RagTag is only compatible with fasta files. 
#Therefore, make sure your output from hifiasm is converted using these commands:
#gfa tools is installed in the herro environment
#gfatools gfa2fa HG02006_hifiasm_haploid.bp.p_ctg.gfa > HG02006_hifiasm_haploid.bp.p_ctg.fasta

# scaffold a query assembly
ragtag.py scaffold ${Ref} ${phased_assembly}

#ragtag doesnt have an option to define the output name, so I need to rename the output folder 
#before running it for the next haplotype.
mv ragtag_output/ ragtag_output_hap1

ragtag.py scaffold ${Ref} ${phased_assembly_2}

mv ragtag_output/ ragtag_output_hap2

conda deactivate

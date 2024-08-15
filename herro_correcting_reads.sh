#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --mem=400G
#SBATCH --gpus-per-node=1
#SBATCH --job-name="DC_HG01946"
#SBATCH --output=Dorado_correct.out
#SBATCH --error=Dorado_correct.err
#SBATCH --export=NONE

##########THIS SCRIPT WAS COPIED AND ADAPTED FROM DAN MACGUIGAN. ##########################
#IMPORTANT NOTES: herro corrector is an AI developped software that is supposed to do wonders with the correction step. 
#After correcting with this method, one can use the ONT ULR as if they were Pac-Bio HiFi reads, and use verkko assembler and hifiasm assembler. 
#This is great for using better assemblers. 

# ABOUT ################################################
# author: Dan MacGuigan
# designed to work on UB CCR cluster

# script to run Nanopore read correction using HERRO (via Dorado)
# requires basecalled reads fastq file
# outputs a fasta file of corrected reads

# example usage:
# first modify variables in the INPUTS section
# then run "sbatch AISO_Dorado_readCorrection.sh"

#### notes from Dorado Github: ############
# Dorado supports single-read error correction with the 
# integration of the HERRO algorithm. HERRO uses all-vs-all
# alignment followed by haplotype-aware correction using a 
# deep learning model to achieve higher single-read accuracies. 
# The corrected reads are primarily useful for generating de 
# novo assemblies of diploid organisms.

# Dorado correct only supports FASTX(.gz) as the input and 
# generates a FASTA file as output. An index file is generated 
# for the input FASTX file in the same folder unless one is already 
# present. Please ensure that the folder with the input file is 
# writeable by the dorado process and has sufficient disk space 
# (no more than 10GB should be necessary for a whole genome dataset).

# The error correction tool is both compute and memory intensive. 
# As a result, it is best run on a system with multiple high 
# performance CPU cores ( > 64 cores), large system memory ( > 256GB) 
# and a modern GPU with a large VRAM ( > 32GB).

# INPUTS ################################################
WD="/projects/academic/omergokc/Luane/HG01946/herro" # working directory
READS="HG01946_allreads_Apr16.fastq" # uncorrected reads, must be located in WD
MIN_READ_LEN=10000 # min read length (bp), HERRO documentation suggests >10 kb
SP="HG01946" # species or sample prefix
CORR_READS="HG01946_10kbp_HERRO.fasta" # name for corrected reads, must be .fasta
THREADS=56


# SCRIPT ################################################

module load gcc/11.2.0 
module load fastp/0.23.2

# filter by min read length they say 10kbp min length
fastp -i ${READS} -o ${SP}_L${MIN_READ_LEN}.fastq -h ${SP}_L${MIN_READ_LEN}.html -j ${SP}_L${MIN_READ_LEN}.json -l ${MIN_READ_LEN} --disable_quality_filtering --disable_adapter_trimming --disable_trim_poly_g -V -w 8

# convert to bgzip
module load htslib
#gunzip ${READS}.fastq.gz
#bgzip ${READS}.fastq

gunzip HG01946_L10000.fastq.gz
bgzip HG01946_L10000.fastq

# run read correction
module load cuda/11.8.0
/projects/academic/tkrabben/software/Dorado/dorado-0.7.1-linux-x64/bin/dorado correct -v -t ${THREADS} ${SP}_L${MIN_READ_LEN}.fastq.gz > ${CORR_READS}

# OUTPUT ################################################
echo "read correction complete"

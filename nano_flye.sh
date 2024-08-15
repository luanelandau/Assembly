#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --mem=600G
#SBATCH --job-name=10kbnanoq_HG01946
#SBATCH --output=10kbULRnanoq_HG01946.out
#SBATCH --error=10kbULRnanoq_HG01946.err
#SBATCH --constraint=AVX512
#SBATCH --export=NONE

#These parameters are for using in the University at Buffalo CCR cluster.

####################### NANOQ ##################################
#First we need to correct the reads. nanoq (https://github.com/esteinig/nanoq) is the program I use to correct the reads, but I highly recommend checking out the herro (https://github.com/lbcb-sci/herro) pipeline. I am using the herro for some of the samples, so you should find also a script for it around here.

#For nanoq I have to create the outdirectories.
outdir="./results_nanoq_10kb_HG01946" #name it as wish, depending on the previous or next steps
mkdir $outdir

#I have nanoq installed in a conda environment, so make sure you have nanoq installed to submit this job
eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

#nanoq also allows you to filter the reads by lenght. Therefore, I provide here three different options that might generate different assemblies at the end. There is no perfect choice, it depends on the quality of the reads, read depth, organism, etc. So you may play around and test which works best.
#For filtering for reads above 10000bp (this is my go-to):
nanoq -i ${reads}  -l 10000 > ./${outdir}/HG01946_10kbreads.fastq
gzip ./${outdir}/${sample}_10kbreads.fastq #I am zipping the file to use in the next step - flye

#15kb and above:
#nanoq -i ${reads}  -l 15000 > ./${outdir}/HG01922_15kbreads.fastq
#gzip ./${outdir}/HG01922_15kbreads.fastq

#20kb and above:
#nanoq -i ${reads}  -l 20000 > ./${outdir}/HG01946_20kbreads.fastq
#gzip /projects/academic/omergokc/Luane/HG01946/results_nanoq_10kb_HG01946_ULR/HG01946_10kbreads.fastq

####################### FLYE ##################################
#Flye is a good and fast assembler for ONT nanopore reads (https://github.com/mikolmogorov/Flye). However, it is not at good for humans as other assemblers, such as canu, verkko, or hifiasm. New pipelines have been developped to correct ULR using herro and, using only ONT generated ULR, run verkko and hifiasm, which I highly recommend. Please take time to look over which assembler is the best for your data.

#INPUT
FLYOVLP="1000" #minimum overlap between reads for Flye
ITERATIONS="1" #number of polishing iterations to run the default is 1
reads="/path/to/your/reads.fastq.gz"
size="3.2g" #estimated genome size in GB. This is set for human
outdir="./results_nanoq_10kb_HG01946_ULR_kh" #dont need to define this again, as it was defined for nanoq
threads=32
sample="sample_name"

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

flye --nano-hq ./${outdir}/HG01946_10kbreads.fastq -o ${outdir} -t ${threads}

conda deactivate

# if you need to resume an interrupted assembly, add the --resume option to the above command


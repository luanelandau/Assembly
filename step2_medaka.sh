#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=900G
#SBATCH --job-name=Medaka_HG01946
#SBATCH --output=Medaka_HG01946.out
#SBATCH --error=Medaka_HG01946.err
#SBATCH --constraint=AVX512
#SBATCH --export=NONE

#all these modules were load while troubleshooting for some problems with the script. They do not necessarily have to be loaded but I will leave them here, because it is working.
module load CCRconfig
module load gcc/11.2.0
module load StdEnv/2023.01
module load gentoo/2023.01 
module load gcccore/11.2.0
module load samtools/1.16.1

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)" #this is to load conda
conda activate Assembly #this is the environment in which I have medaka installed. For more information look into: https://github.com/nanoporetech/medaka
#For running this I have to have installed minimap2 and samtools!

#input
reads="path/to/your/reads.fastq.gz" #path to fastq files
assembly="assembly.fasta" #path to draft genome
sample="HG01946_10kbp_flye" #sample name
model="r1041_e82_400bps_sup_g615" #basecalling model used for fastq files
mapthreads=56 #Mapping threads
conthreads=8 #compiling threads Creators of Medaka say over 2 is not useful but use 8 as their example?!?!?!

#mkdir results

#align reads to assembly
#mini_align -i ${reads} -r ${assembly} -P -m -f -p ./results/${sample}_draft -t ${mapthreads}

# you can run lots of jobs like this at the same time by using the --region flag and providing a list of a subset of the contigs found in you draft assembly.
medaka consensus results/${sample}_draft.bam results/${sample}_new2.hdf --model ${model} --threads ${conthreads}
medaka stitch results/${sample}_new2.hdf ${assembly} ${sample}_Medaka_new2.fasta

conda deactivate

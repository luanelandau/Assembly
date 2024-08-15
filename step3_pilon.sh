#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=700G
#SBATCH --job-name=pilon_HG01946
#SBATCH --output=pilon_HG01946.out
#SBATCH --error=pilon_HG01946.err
#SBATCH --export=NONE

#Pilon (https://github.com/broadinstitute/pilon) uses illumina short reads to polish an assembly. For this, you would need to have generated the short reads for the same individual, and a forward and a reverse file. 

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)" #activating conda
source activate Assembly #the conda enviroment I have all softwares installed

#INPUT
REF="/projects/academic/omergokc/Luane/HG01946/medaka/HG01946_10kbp_flye_Medaka_new3.fasta" #Reference assembly to be polished
FORWARD="/projects/academic/omergokc/Luane/HG01946/1000KG_sequences/ERR3988956_1.fastq.gz" #file name for forward read
REVERSE="/projects/academic/omergokc/Luane/HG01946/1000KG_sequences/ERR3988956_2.fastq.gz" #file name for reverse read
prefix="HG01946_10kbp_flye_Medaka"
THRDS=40 #number threads

#STEP 1: map trimmed Illumina reads to reference using bwa. Then conver
bwa-mem2 index $REF
bwa-mem2 mem -t ${THRDS} $REF ${FORWARD} ${REVERSE} | samtools view -@ ${THRDS} -S -q 15 -b | samtools sort -@ ${THRDS} -o ${prefix}_pilon_illumina.mapped.bam

##Index bam
samtools index -@ ${THRDS} ${prefix}_pilon_illumina.mapped.bam

mkdir pilon_out
#you have to have pilon installed with java to run this next step. In this following case, I am calling pilon exactly where I installed it in my folder on the server.
java -Xmx700G -jar /projects/academic/omergokc/Luane/pilon-1.24.jar --genome ${REF} --frags  ${prefix}_pilon_illumina.mapped.bam --output ${prefix}_pilon --outdir ./pilon$_out --fix all --threads ${THRDS}

#OUTPUT
#.bam file for mapped Illumina reads
# indexed files from bwa and samtools
#Updated assembly corrected with short reads (directed to ./out within working directory)

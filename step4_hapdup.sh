#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=500G
#SBATCH --job-name=HapDup_Opossum
#SBATCH --output=HapDup_Opossum.out
#SBATCH --error=HapDup_Opossum.err
#SBATCH --export=NONE

#input
reads="/projects/academic/omergokc/Luane/Opossum1/Opossum1_allreads.fastq.gz" #path to fastq files
assembly="/projects/academic/omergokc/Luane/Opossum1/medaka/Opossum1_10kbp_flye_Medaka.fasta" #path to draft genome
sample="Opossum1_10kbp_flyeMedakapilon" #sample name
HD_DIR="/projects/academic/omergokc/Luane/softwares"


eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly


minimap2 -ax map-ont -t 30 ${assembly} ${reads} | samtools sort -@ 4 -m 4G > ${sample}_lr_mapping.bam
samtools index -@ 4 ${sample}_lr_mapping.bam

conda deactivate

export SINGULARITY_CACHEDIR="/projects/academic/omergokc/Luane/singularity3"
export SINGULARITY_LOCALCACHEDIR="${SINGULARITY_CACHEDIR}/cache"
export SINGULARITY_TMPDIR="${SINGULARITY_CACHEDIR}/tmp"
export APPTAINER_CACHEDIR="${SINGULARITY_CACHEDIR}"
export APPTAINER_TMPDIR="${SINGULARITY_TMPDIR}"
singularity exec --bind /cvmfs:/cvmfs:ro --bind /util:/util:ro --bind /vscratch:/vscratch --bind /projects:/projects --bind /scratch:/scratch docker://mkolmogo/hapdup:0.12 hapdup --assembly ${assembly} --bam ${sample}_lr_mapping.bam --out-dir hapdup -t 40 --rtype ont

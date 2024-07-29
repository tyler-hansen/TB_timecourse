#!/bin/bash
#SBATCH --job-name=STAR_aligner_run1_test
#SBATCH --output=STAR_aligner_run1_test.out
#SBATCH --error=STAR_aligner_run1_test.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=64G
#SBATCH --account=pi-lbarreiro

#paths
hg38_ref='/project/lbarreiro/USERS/tyler/general_purpose/hg38'
STAR_dir='/project/lbarreiro/USERS/tyler/general_purpose/STAR_builds'

#Run1
INPUT_DIR='/project/lbarreiro/USERS/ilia/time_course/trimmed_fastqs/trimmed_renamed'
OUTPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads/run1'

#load STAR software through my STAR conda env
module restore anaconda
source activate STAR

#align trimmed reads with STAR - Run 1
STAR --runThreadN 12 --genomeDir ${STAR_dir}/STAR_hg38_gencode.v44_100 \
    --outSAMtype BAM Unsorted \
    --readFilesIn ${INPUT_DIR}/10_EU118_T0_Mtb_MOI_5.13579.fq.gz \
    --outFileNamePrefix $OUTPUT_DIR \
    --readFilesCommand zcat

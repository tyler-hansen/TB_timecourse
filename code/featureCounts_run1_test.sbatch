#!/bin/bash
#SBATCH --job-name=feature_counts_run1_test
#SBATCH --output=feature_counts_run1_test.out
#SBATCH --error=feature_counts_run1_test.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=64G
#SBATCH --account=pi-lbarreiro

#Run1
INPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads/run1'
OUTPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/counts'

GTF='/project/lbarreiro/USERS/tyler/general_purpose/hg38/gencode.v44.primary_assembly.basic.annotation.gtf'

#load STAR software through my STAR conda env
module restore anaconda
source activate subread

#run1 is reverse stranded, so set -s 2. Use exon as feature and gene_id as metafeature. report time for estimating resources. 
time featureCounts -T 12 -s 2 -F GTF -a $GTF -t exon -g gene_id -o ${INPUT_DIR}/run1_counts_test.txt ${INPUT_DIR}/14_EU144_T4*Aligned.out.bam
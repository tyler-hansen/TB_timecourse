#!/bin/bash
#SBATCH --job-name=feature_counts
#SBATCH --output=feature_counts.out
#SBATCH --error=feature_counts.err
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64G
#SBATCH --account=pi-lbarreiro

INPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads'
OUTPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/counts'

GTF='/project/lbarreiro/USERS/tyler/general_purpose/hg38/gencode.v44.primary_assembly.basic.annotation.gtf'
GTF_PC='/project/lbarreiro/USERS/tyler/general_purpose/hg38/gencode.v44.primary_assembly.basic.annotation.protein_coding-only.gtf'

#load STAR software through my STAR conda env
module restore anaconda
source activate subread

#Make two versions, one with full GTF, other with protein coding only GTF. 
#Use exon as feature and gene_id as metafeature.
#Run1,2,3 are reverse stranded, so set -s 2. 

for i in {1..3}
do
    featureCounts -T 4 -s 2 -F GTF -a $GTF -t exon -g gene_id -o ${OUTPUT_DIR}/run${i}_counts.txt ${INPUT_DIR}/run${i}/*Aligned.out.bam
    featureCounts -T 4 -s 2 -F GTF -a $GTF_PC -t exon -g gene_id -o ${OUTPUT_DIR}/run${i}_counts_protein-coding.txt ${INPUT_DIR}/run${i}/*Aligned.out.bam
done

#Run4 is stranded, so set -s 1.
featureCounts -T 4 -s 1 -F GTF -a $GTF -t exon -g gene_id -o ${OUTPUT_DIR}/run4_counts.txt ${INPUT_DIR}/run4/*Aligned.out.bam
featureCounts -T 4 -s 1 -F GTF -a $GTF_PC -t exon -g gene_id -o ${OUTPUT_DIR}/run4_counts_protein-coding.txt ${INPUT_DIR}/run4/*Aligned.out.bam
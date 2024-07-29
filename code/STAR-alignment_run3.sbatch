#!/bin/bash
#SBATCH --job-name=STAR_aligner_run3
#SBATCH --output=STAR_aligner_run3.out
#SBATCH --error=STAR_aligner_run3.err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=64G
#SBATCH --account=pi-lbarreiro

#paths
STAR_dir='/project/lbarreiro/USERS/tyler/general_purpose/STAR_builds'

#Run1
INPUT_DIR='/project/lbarreiro/USERS/ilia/time_course/trimmed_fastq_run3'
OUTPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads/run3'

#load STAR software through my STAR conda env
module restore anaconda
source activate STAR

for i in ${INPUT_DIR}/*L001*.fq.gz; do
# Check if the file is a regular file
    if [ -f "$i" ]; then
        # Extract the basename of the input file
        name=$(basename -s _L001_R1_001_trimmed.fq.gz $i)
        #align trimmed reads with STAR - Run 3
        STAR --runThreadN 12 --genomeDir ${STAR_dir}/STAR_hg38_gencode.v44_100 \
            --outSAMtype BAM Unsorted \
            --readFilesIn ${INPUT_DIR}/${name}_L001_R1_001_trimmed.fq.gz,${INPUT_DIR}/${name}_L002_R1_001_trimmed.fq.gz,${INPUT_DIR}/${name}_L003_R1_001_trimmed.fq.gz,${INPUT_DIR}/${name}_L004_R1_001_trimmed.fq.gz \
            --outFileNamePrefix ${OUTPUT_DIR}/${name}_ \
            --readFilesCommand zcat
        echo "Aligned: ${name}"
    fi
done
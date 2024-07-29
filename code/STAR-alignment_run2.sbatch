#!/bin/bash
#SBATCH --job-name=STAR_aligner_run2
#SBATCH --output=STAR_aligner_run2.out
#SBATCH --error=STAR_aligner_run2.err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
#SBATCH --mem=64G
#SBATCH --account=pi-lbarreiro

#paths
STAR_dir='/project/lbarreiro/USERS/tyler/general_purpose/STAR_builds'

#Run2
INPUT_DIR='/project/lbarreiro/USERS/ilia/time_course/trimmed_fastq_run2'
OUTPUT_DIR='/project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads/run2'

#load STAR software through my STAR conda env
module restore anaconda
source activate STAR

#Text file containing filenames to skip
skip_file="/project/lbarreiro/USERS/tyler/TB_timecourse/code/skipfile_run2.txt"

for i in ${INPUT_DIR}/*.fq.gz; do
# Check if the file is a regular file
    if [ -f "$i" ]; then
        # Extract the basename of the input file
        name=$(basename -s .fq.gz $i)
        
        # Check if the filename is in the skip file
        if grep -Fq "$(basename $i)" "$skip_file"; then
            echo "Skipping file: ${name}"
        else
            #align trimmed reads with STAR - Run 1
            STAR --runThreadN 12 --genomeDir ${STAR_dir}/STAR_hg38_gencode.v44_100 \
                --outSAMtype BAM Unsorted \
                --readFilesIn ${INPUT_DIR}/${name}.fq.gz \
                --outFileNamePrefix ${OUTPUT_DIR}/${name}_ \
                --readFilesCommand zcat
            echo "Aligned: ${name}"
        fi
    fi
done
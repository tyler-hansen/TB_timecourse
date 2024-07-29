#!/bin/bash
#SBATCH --job-name=samtools_index_and_sort-run4-fix
#SBATCH --output=%x.out
#SBATCH --error=%x.err
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=64G
#SBATCH --account=pi-lbarreiro

# This script sorts and indexes mapped reads. 

#load software through conda
module restore anaconda
source activate samtools

#change directory
cd /project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads

#run4 only:
for n in 4
do
    for file in run${n}/*out.bam # this ensures only unsorted files are processed. 
    do    
	samtools sort -@ 24 -o run${n}/$(basename -s .bam $file).sorted.bam $file
        samtools index -@ 24 run${n}/$(basename -s .bam $file).sorted.bam run${n}/$(basename -s .bam $file).sorted.bam.bai
	if [ -e "run${n}/$(basename -s .bam $file).sorted.bam" ]; then
	        # Remove the file
		rm "$file"
	fi
     done
done

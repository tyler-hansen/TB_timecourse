#!/bin/bash
#SBATCH --job-name=infer_experiments
#SBATCH --output=infer_experiments.out
#SBATCH --error=infer_experiments.err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --account=pi-lbarreiro

# This script reports whether the mapped reads are standed and, if so, in which direction. This is necessary to know for featureCounts. 

#load software through conda
module restore anaconda
source activate rseqc

#paths
BED='/project/lbarreiro/USERS/tyler/general_purpose/hg38/gencode.v44.primary_assembly.basic.annotation.bed'
# I generated this with the following code (using bedops): awk '{ if ($0 ~ "transcript_id") print $0; else print $0" transcript_id \"\";"; }' gencode.v44.primary_assembly.basic.annotation.gtf | gtf2bed - > gencode.v44.primary_assembly.basic.annotation.bed

#change directory
cd /project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads

#loop through runs 1-4:
for n in {1..4}
do
    #remove incase this is a re-run
    rm ../../results/infer_experiments_output_run${n}.txt

    #loop through files and run infer_experiment.py from the rseqc package
    for i in run${n}/*Aligned.out.bam
        do infer_experiment.py -s 10000 -r $BED -i $i >> ../../results/infer_experiments_output_run${n}.txt
    done
done

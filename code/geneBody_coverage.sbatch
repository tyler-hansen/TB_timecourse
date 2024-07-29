#!/bin/bash
#SBATCH --job-name=geneBody_coverage
#SBATCH --output=geneBody_coverage.out
#SBATCH --error=geneBody_coverage.err
#SBATCH --time=1-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --account=pi-lbarreiro

# This script reports the gene body coverage of mapped reads. 

#load software through conda
module restore anaconda
source activate rseqc

#paths
BED='/project/lbarreiro/USERS/tyler/general_purpose/hg38/gencode.v44.knownGene.UCSC.bed'
# I generated this from UCSC genome browser. One gene per line. 

#change directory
cd /project/lbarreiro/USERS/tyler/TB_timecourse/data/mapped_reads

#loop through runs 1-4:
for n in {1..4}
do
    geneBody_coverage.py -r $BED -i run${n}/ -o ../../results/geneBody_coverage/run${n}_geneBody_coverage
done

#summarize with multiQC
multiqc ../../results/geneBody_coverage/ -n ../../results/geneBody_coverage/multiqc_report_all-runs.html

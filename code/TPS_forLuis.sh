# Use TPS to identify most optimal timepoints
# Tyler Hansen
# 02.21.2024

# I installed TPS in "Computation/programs/external_programs/TPS_code". It requires installation of python packages, which I installed via conda. Prior to running, launch TPS conda env. 

# It requires data formatted as a csv, where gene is the first column, timepoint 1 value is second, and so on. Row 1 is the column names for the timepoints. 

# Using R, adjust filtered counts. See format_cts_for_TPS.Rmd. 

# run TPS using two strageies: 
    # 1) For each individual for Mtb only
    # 2) Average across individuals for Mtb only. 


#run TPS
#instructions: 
#   Command Line Arguments:
        #-n number of time points: List of number of time point samples taken(List of integers), example: 8 10 12 (*)
        #-d data: Filename of dense time series provided format described below, example: input.data.txt (*)
        #-p plots: Path of folder to put splineplots , default: plots/  ,
        #-e error_plot: Filename of plot to produce comparing the average error for different numbers of time point samples (will only be produced if more than one value is given for n)
        #        , default: avg_error.pdf
        #-m method: Initilization method as described in supplement, options: metricA, metricB,metricC, max_dist (recommended)
        #*required arguements

OUT='../results/TPS_cts'

for file in ../results/TPS_cts/cts*.txt
do
    base=$(basename -s ".txt" $file) 
    python ../../../programs/external_programs/TPS_code/TPS.py -d $file -n {4..20} -m max_dist \
        -e ${OUT}/average_error_${base} 
done
#!/bin/bash

#directories
source_directory="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/kallisto_samples"
key_value_file="/mnt/c/Users/tylerhansen/OneDrive - The University of Chicago/5_Postdoc/Computation/projects/TB_timecourse/data/sample_conversion_all_samples_combined.txt"

# Check if the source directory exists
if [ ! -d "$source_directory" ]; then
    echo "Error: Source directory not found: $source_directory"
    exit 1
fi

# Check if the key-value file exists
if [ ! -f "$key_value_file" ]; then
    echo "Error: Key-value file not found: $key_value_file"
    exit 1
fi

# Loop through each line in the key-value file
while IFS= read -r line; do
    # Split the line into key and value
    key=$(echo "$line" | cut -f 1)
    value=$(echo "$line" | cut -f 2)

    # Check if the folder with the old name exists
    if [ -d "${source_directory}/kallisto_$key" ]; then
        # Rename the folder
        mv "${source_directory}/kallisto_$key" "${source_directory}/kallisto_$value"
        echo "Renamed kallisto_$key to kallisto_$value"
    else
        echo "Error: Folder kallisto_$key not found"
    fi
done < <(tail -n +2 "$key_value_file") #skips the first line, which is a header
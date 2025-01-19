#!/bin/bash

# Variable to store the previous job ID
previous_job_id=""

# Output file
data="Dubcova2_6"
output_file="${data}_output.txt"

# Clear the output file at the beginning
> "$output_file"

# Loop through ntasks and submit jobs sequentially
# for ntasks in "${ntasks_list[@]}"; do
for ntasks in {1..64}; do

previous_job_id=$(sbatch --ntasks=$ntasks --output="$output_file" --open-mode=append skibidi.sh "$data" | awk '{print $4}')
 
    echo "Submitted job with ntasks=$ntasks, waiting for job ID $previous_job_id to finish"
done

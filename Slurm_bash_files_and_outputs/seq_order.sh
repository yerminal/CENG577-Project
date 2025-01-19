#!/bin/bash

# Variable to store the previous job ID
previous_job_id=""

# Output file
data="Kuu_5"
output_file="${data}_output.txt"

# Clear the output file at the beginning
> "$output_file"

# Function to wait for a job to finish
wait_for_job() {
    local job_id=$1
    while true; do
        # Check if the job is still in the queue
        if squeue -j "$job_id" 2>/dev/null | grep -q "$job_id"; then
            echo "Job $job_id is still running. Waiting..."
            sleep 0.5 # Check every 5 seconds
        else
            echo "Job $job_id has completed."
            break
        fi
    done
}

# Loop through ntasks and submit jobs sequentially
for ntasks in {1..64}; do
    if [ -z "$previous_job_id" ]; then
        # Submit the first job without any dependency
        previous_job_id=$(sbatch --ntasks=$ntasks --output="$output_file" --open-mode=append skibidi.sh "$data" | awk '{print $4}')
    else
        # Wait for the previous job to finish
        wait_for_job "$previous_job_id"

        # Sleep for an additional second before submitting the next job
        sleep 1

        # Submit subsequent jobs with a dependency on the previous job
        previous_job_id=$(sbatch --ntasks=$ntasks --output="$output_file" --open-mode=append skibidi.sh "$data" | awk '{print $4}')
    fi
    echo "Submitted job with ntasks=$ntasks, waiting for job ID $previous_job_id to finish"
done

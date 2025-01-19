#!/bin/bash
# This file test.sh is a sample script to run many jobs using multiple cores.
# To have output files defined as slurm-%A_%a, run your script with sbatch -o slurm-%A_%a.out test.sh
# command where %A is the job ID and %a is {0..ntasks}
# To have one output file for all simultaneous tasks, run your script with sbatch test.sh

# Set the partition where the job will run
#SBATCH --partition=577Q

# Define the name of the submitted job
#SBATCH --job-name=slurm_test

# Default output file if the script is run with sbatch test.sh
#SBATCH --output=slurm_output-%j.txt

# Set the number of nodes and tasks per node (this will run 4 tasks simultaneously)
#SBATCH --nodes=1
#SBATCH --ntasks=64

# Send mail to this address (replace with your email address)
#SBATCH --mail-user=e244307@metu.edu.tr

# Mail alert at start, end, and abortion of execution
#SBATCH --mail-type=ALL

# Launch the command/application
mpirun ~/file_op/a.out Kuu_5
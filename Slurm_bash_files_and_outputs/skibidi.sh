#!/bin/bash

# Set the partition where the job will run
#SBATCH --partition=577Q

# Define the name of the submitted job
#SBATCH --job-name=Hani_5Dk

# Set the number of nodes
#SBATCH --nodes=1

# Launch the command/application and append its output
mpirun ~/file_op/a.out $1

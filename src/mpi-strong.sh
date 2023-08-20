#!/bin/bash

# Load the required modules
module load gcc
module load openmpi

make

# Define the maximum number of nodes
max_nodes=16

# Loop through increasing numbers of nodes
for ((nodes=1; nodes<=$max_nodes; nodes++)); do
    output_file="output_$nodes.txt"
    echo "Running with $nodes nodes..."
    mpirun -n $nodes ./main -s 0.067 > "output-mpi-strong/$output_file" 2>&1
    echo "Output saved to output-mpi-strong/$output_file"
done

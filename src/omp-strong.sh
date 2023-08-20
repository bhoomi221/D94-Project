#!/bin/bash

# Load the required modules
module load gcc
module load openmpi

make

# Define the maximum number of nodes
max_threads=16

# Loop through increasing numbers of nodes
for ((threads=1; threads<=$max_threads; threads++)); do
    output_file="output_$threads.txt"
    echo "Running with $threads threads..."
    export OMP_NUM_THREADS=$threads
    { time ./main -s 0.067 > "output-omp-strong/$output_file"; } 2> "output-omp-strong/$output_file"
    echo "Output saved to output-omp-strong/$output_file"
done
export OMP_NUM_THREADS=1
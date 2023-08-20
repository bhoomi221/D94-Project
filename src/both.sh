#!/bin/bash

# Load the required modules
module load gcc
module load openmpi

make

# Specify the range of thread and MPI process combinations to test
THREADS=(1 2 4 8)
MPI_PROCESSES=(2 4 8 16)
# Create a CSV file to store results

echo "Threads,MPI_Processes,Time" > results.csv

# Loop through thread and MPI process combinations
for t in "${THREADS[@]}"; do
    for p in "${MPI_PROCESSES[@]}"; do
        # Set OpenMP environment variables
        export OMP_NUM_THREADS=$t

        # Set the command to run your parallel executable
        RUN_COMMAND="mpirun -np $((p)) ./main -s 0.067"

        # Run the parallel program and capture execution time
        execution_time=$(TIMEFORMAT=%R; { time $RUN_COMMAND; } 2>&1)

        # Append results to CSV file
        echo "$t,$p,$execution_time" >> results.csv
    done
done

echo "Experiments complete."

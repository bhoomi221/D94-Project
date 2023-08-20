# D94-Project
SPH Project that is parallelized with OMP and MPI


Makefile\
Run make (only works in cluster) remove -lmpi_cxx to run locally\
./main to run \
use ./main -h flag to see options for parameter settings

To visualize the output\
make view

Can use the shell scripts to reproduce the results\
Results found in:
- src/output-mpi-strong
- src/output-omp-strong
- src/results.csv

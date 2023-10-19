
.PHONY: all clean 

exe: main
CC       = g++ -fopenmp
CFLAGS   = -std=c++14 -Wall -g
OPTFLAGS = -O3
LIBS     = -lm -lmpi -lmpi_cxx


# =======

main: main.o params.o math.o integration.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS)

main.o: main.cpp params.hpp particle.hpp math.hpp integration.hpp

params.o: params.cpp params.hpp
math.o: math.cpp math.hpp particle.hpp params.hpp
integration.o: integration.cpp integration.hpp particle.hpp

%.o: %.cpp
	$(CC) -c $(CFLAGS) $(OPTFLAGS) $<

# =======
view: 
	java -jar Bouncy.jar run.out
mpi:
	mpirun -n 3 ./main
tag: 
	mpirun --tag-output -n 3 ./main	
clean:
	rm -f *~ *.o

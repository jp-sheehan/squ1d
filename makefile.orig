#
#   makefile for Flux
#
OBJ = main_PIC.o solver.o mathfunctions.o mesh.o initializesolver.o writeoutput.o boundary.o solvertype.o   #..Objects to compile
CC = mpic++ #g++ #..Select the Compiler
LOADER = mpic++ #g++ #gfortran   #..Select the Loader
IFLAGS = -I"/opt/intel/mkl/include/"#..Include flags for MKL
FLIB = /opt/intel/mkl/lib/intel64#..Library location for MKL
LFLAGS = -L/opt/intel/mkl/lib/intel64/ -lmkl_intel_ilp64 -lmkl_core -lpthread -lmkl_sequential -lm -std=gnu++0x -lstdc++ -Wall -O2 -Wextra -Wformat=2 -lmpi -g -pg #..Compiler flag for MKL
CFLAGS = 

PIC.exe:    $(OBJ) $(LIBRARIES)
	$(LOADER) $(CFLAGS) -m64 -w -DMKL_ILP64 $(OBJ) $(LFLAGS) -o PIC.exe 
	export LD_LIBRARY_PATH="$(FLIB)":$(LD_LIBRARY_PATH);

.cpp.o:
	$(CC) -c -m64 -w -DMKL_ILP64 $(CFLAGS) $< $(IFLAGS) $(LIBRARIES)  $(LFLAGS) -o $@ 

clean:  
	rm *.o *.exe 

cleandat: 
	-@rm *.dat *.out *.gif *.png *.jpg
	ls | grep -P ".*\.o[0-9]+$$" | xargs -d"\n" rm

cleanjunk:
	rm *~ .*~ .*.swp .*.swo

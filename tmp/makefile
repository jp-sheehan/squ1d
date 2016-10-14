#
#   makefile for umich FLUX
#
OBJ = main_PIC.o solver.o mathfunctions.o mesh.o initializesolver.o writeoutput.o boundary.o solvertype.o   #..Objects to compile
CC = mpiCC #g++                 #..Compiler
LOADER = mpiCC #g++ #gfortran   #..Loader
IFLAGS = -I"/home/software/rhel6/intel-14.0/mkl/include" #..Include flags for MKL
FLIB = "/home/software/rhel6/intel-14.0/mkl/lib/intel64"     #..Library location for MKL
LFLAGS = -L"/home/software/rhel6/intel-14.0/mkl/lib/intel64" -L"/usr/cac/rhel6/intel-14.0/lib/intel64" -L"/home/software/rhel6/intel-2013.1.046/lib/intel64" -L"/opt/intel/mic/coi/host-linux-release/lib" -L"/opt/intel/mic/myo/lib" -L"/usr/local/torque/lib" -L"/usr/cac/rhel6/openmpi/1.6.5-intel-14.0.2/lib" -L"/usr/cac/rhel6/openmpi/1.6.5-intel-14.0.2/lib/openmpi" -L"/opt/mellanox/mxm/lib" -lpthread -lm -std=gnu++0x -lstdc++ -Wall -Wextra -Wformat=2 -O2 -lmpi -Wl,--start-group $(MKL_LIB)/libmkl_intel_ilp64.a $(MKL_LIB)/libmkl_core.a $(MKL_LIB)/libmkl_sequential.a -Wl,--end-group#..Compiler flags


PIC.exe:  $(OBJ) $(LIBRARIES)
	$(LOADER) -m64 -w -DMKL_ILP64 $(OBJ) $(LFLAGS) -o PIC.exe 
	export LD_LIBRARY_PATH="$(FLIB)":$(LD_LIBRARY_PATH);

.cpp.o:
	$(CC) -c -m64 -w -DMKL_ILP64 $< $(IFLAGS) $(LIBRARIES)  $(LFLAGS) -o $@ 

clean:  
	rm *.o *.exe 

cleandat: 
	rm *.dat *.out     

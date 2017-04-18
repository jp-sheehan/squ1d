# Notes on SQu1D

These are an assortment of notes compiled to aid new users in understanding
and executing SQu1D.  They are more freeform than the official readme and
get into the underlying operation of the code, not just how to run it.

## Installation

[Intel MKL Link Advisor]
(https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/)

MKL install location: /opt/intel/parallel_studio_xe_2016.3.067

See note: Academic Research Performance Libraries from Intel

### To install

1. update and upgrade packages
	1. [Fix screen resolution in VirtualBox](http://askubuntu.com/questions/73589/higher-screen-resolution-for-virtualbox)
	2. [If you are still having screen size problems in VirtualBox](http://askubuntu.com/questions/588943/experiencing-small-resolution-issue-in-ubuntu-14-04-2-with-virtualbox-getting-s)
2. [Download](https://software.intel.com/en-us/qualify-for-free-software/academicresearcher) and install inlet MPI and MKL libraries 
3. install MPICH2: `sudo apt-get install libcr-dev mpich2 mpich2-doc`
4. Install LAPACK 
	1. `sudo apt-get install libblas-dev liblapack-dev liblapacke-dev`
	2. Download latest tar ball
	3. unzip
	4. copy `make.inc.example` to `make.inc`
	5. run `make blaslib` [to ensure a successful build](http://icl.cs.utk.edu/lapack-forum/viewtopic.php?f=4&t=3029)
	6. run `make`
	7. Change directories to the LAPACKE folder
	8. run `make`
5. Download SQu1D files: 
6. `make`
	1. for intel compiler and libraries, just run `make`
	2. for lapack, run `make -f "makefile_LP"`
	3. on flux, run `make -f "makefile_FLUX"`
7. Add environmental variable `export LD_LIBRARY_PATH=/opt/intel/mkl/lib/intel64`
8. run `./PIC.exe`


### Graphing Software

1. [TecPlot](http://www.tecplot.com/)
2. [GnuPlot](http://www.gnuplot.info/)

## Execution

1. Execution folder must contain:
	1. PIC.exe
	2. Cross section folder
	3. SolverType.inp
	4. SolverInput.inp
2. Run `./PIC.exe`

## Output files

(****: species name; #: number)

### ****Output_cField#.dat and ****Output_pField#.dat
 
* Gives values of cell centered and edge values of parameters
* position (2D), density, velocity (3D), temperature (3D), energy (1D), electric field (3D), magnetic field (3D), heat flux (I think?) (3D)

### ****Particles#.dat
 
* Gives values of individual particles
* Position (2D), velocity (3D), energy
* Provides some parameters of the simulation

### ****phasespace#.dat
 
* Binned phase space data
* Position (1D), velocity (1D), f (?), Vi (?), fx, fy, fz

### ****totalvdf#.dat
 
* Velocity distribution function
* total velocity, total vdf, Vx, fx, Vy, fy, Vz, fz

### gmon.out
 
* binary file

### restart#.out
 
* A few values without labels

### Output_Field#.dat
 
* A variety of values
* Position (2D), number density, charge density, potential energy, kinetic energy, total energy

### pmesh.out
 
* position of edges
* 2D mesh, X and Y

### cmesh.out
 
* position of centers
* 2D mesh, X and Y

### connect.out
 
* Empty

## Viewing results

* GnuPlot 
	* Linux: install gnuplot (using version 4.6 here) and gnuplot-x11 (for graphics)
* Use script plot.sh for quick plotting to view results

## Files and their Class Definitions

1. solvertype.cpp
	1. sType
2. initializesolver.cpp
	1. initializeSolver
3. variables.h
	1. solverVars
	2. fields
	3. boundvars
	4. spclvars
	5. particles
	6. contnm
4. mesh.h
	1. Mesh
5. solver.h
	1. solver
6. boundary.h
	1. boundary
7. mathfunctions.h
	1. mathFunctions
8. writeoutput.h
	1. writeOutput
9. exprtk.hpp
	1. external package
	2. equation parser

## What Happens During Execution

1. run the program (main_PIC.cpp)
	1. Initialize the MPI environment
		1. Uses the open MPI environment through the wrapper compiler mpic++.  
	2. Set flags for the compiler type and simulation type
	3. Initialize the solver (initializeSolver.cpp)
		1.  Creates object with methods but no data
		2. Only used to set the solver flag
	4. Initialize the solver type (solverType.cpp)
	5. Set the solver flag (initializeSolver::readSolver)
		1. Opens SolverType.inp and reads first line
	6. Execute one of the solvers (solverType.cpp)
		1. Electrostatic 1D PIC
			1. This is the main process where the simulation happens
			2. Define some variables
				1. Create uninitialized vector of particles (variables.h)
				2. Create uninitialized vector of continuum parameters (variables.h)
				3. Create uninitialized MPI variables (variables.h)
			3. Initialize the code (initializesolver::readdata)
				1. Set MPI size and rank
				2. Seed random number generator
				3. Create uninitialized solver (initializeSolver.cpp)
				4. Create uninitialized solver variables (variables.h)
				5. Create uninitialized mesh (mesh.h)
				6. Create uninitialized fields (variables.h)
				7. Create uninitialized solver (solver.h)
				8. Create uninitialized boundary (boundary.h)
				9. Create math functions (mathFunctions.h)
				10. Create uninitialized boundary variables (variables.h)
				11. Create uninitialized special variables (variables.h)
				12. Read data into solver variables, mesh, boundaries, and special variables (initializeSolver::readdata)
					1. Define input file names (hard coded)
					2. Read a bunch of parameters
					3. Set the mesh flags
						1. write
						2. symmetry
						3. solver type
					4. Create mesh parameters
						1. Set length, starting point, boundaries, number of points on mesh
						2. Set cell weight, inlet area, and number of ghost cells
						3. Set particle distribution (What does this mean?)
				13. Read data into information for outputs in solver variables
				14. Read data into solver specific information into solver variables
				15. Read data into particle weight, interpolation scheme,
				16. Read data into particle mover algorithm
				17. Read data into background charge density
				18. Read solver input file
				19. Read data into initial condition variables
			4. Initialize output object (writeOutput.cpp)
			5. Create mesh (mesh::generatePICmesh)
				1. Create mesh points and centers
				2. Create array of mesh point and center neighbors (centers to points, centers to centers, points to centers, points to points)
			6. Initialize PIC (initializeSolver::initializePIC)
				1. Initialize phi and E and B
				2. Set mesh area
				3. Load restart file, if being used (initializeSolver::initializeRestart)
				4. Initialize particle locations
				5. Initialize continuum values
				6. Place particles on mesh
				7. Apply boundary conditions
			7. Seed particles (boundary::seedParticles)  
				1. Set density, velocity, and temperature (at the boundary, I think)
				2. Call the seeding algorithm
			8. Main Loop (inside solvertype::electrostatic1DPIC)
				1. Solve Poisson's equation for phi (solver::poisson1D)
					1. Depends on boundary conditions
					2. Create the vectors of the solvertridiagonal matrix elements
					3. Solve the matrix equation
					4. Assign values to phi
					5. Set E field BC, if appropriate
				2. Calculate the E field from the electric potential (solver::phitoE)/in
				3. Apply E field boundary (solver::applyEfieldBoundary)
				4. Update the particle velocities (solver::updatePartVel)
					1. Get the E field at the particle location (solver::weighEfield)
					2. Get new particle velocity
				5. Redistribute particles (solver::redistributeparticles)
				6. Output files
					1. Get continuum values (n, u, v, ...)
						1. Weigh the contributions from the individual particles
						2. Smooth the data
						3. Apply continuum boundary conditions
					2. Write PIC fields (and continuum values) for each particle type (writeOutput::writePICField)
					3. Write the individual particle information (writeOutput::writeParticles)
					4. Write the vdfs (writeOutput::findvdf)
					5. Write the phase space (writeOutput::findPhaseSpace)
					6. Write restart files (writeOutput::writeRestart)
					7. Update particle positions
			2. Electromagnetic 2D PIC
				1. Not implemented
			3. Euler 1D
				1. Not implemented
			4. Euler 2D
				1. Not implemented
		7. Close out the MPI environment

## Using Supercomputers

### On FLUX

Modules to load
1. gcc/4.8.5     
2. intel/16.0.3   
3. openmpi/1.10.2/intel/16.0.3
4. mkl/11.3.3

Interactive node: `qsub -I -V -A sheehanj_fluxod -q fluxod -l nodes=1:ppn=1,pmem=1gb,walltime=1:00:00,qos=flux`

alias: `interact`


### On Stampede

`ssh sheehanj@stampede.tacc.utexas.edu`

interactive: `idev -m 15`  (15 minutes)

load module cxx11



# Notes on SQu1D

These are an assortment of notes compiled to aid new users in understanding
and executing SQu1D.  They are more freeform than the official readme and
get into the underlying operation of the code, not just how to run it.

## Installation

[Intel MKL Link Advisor](https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/)
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




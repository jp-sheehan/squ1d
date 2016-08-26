       _____ ____        ______
      / ___// __ \__  __< / __ \
      \__ \/ / / / / / / / / / /
     ___/ / /_/ / /_/ / / /_/ /
    /____/\___\_\__,_/_/_____/
---
1. [Introduction](#introduction)
2. [Installing](#installing)
3. [Compiling](#compiling)
4. [Input Files](#input)
5. [Running](#running)
6. [References](#references)

---

<a name="introduction"></a>
## 1. Introduction

The Self-consistent Quasi-1D (SQu1D) code is an object-oriented
code written in the C++ language and is intended to be run in
parallel on multiple processors. The code is used for simulating
magnetic field guided plasmas using a quasi-one-dimensional
fully kinetic approach. The code is a variant of Particle-In-Cell (PIC) technique.

<a name="installing"></a>
## 2. Installing

SQu1D requires Message Passing Interface (MPI) and linear algebra libraries.
It is preferrable to use the [Intel Performance Libraries][intel-libs], which is
free for academic research.  These libraries give SQu1D the best perforamnce.
Download and install the MPI and MKL libraries.

Alternatively, [MPICH2][mpich2] and [LAPACK][lapack] libraries can be used.
To install MPICH2 on a Debian Linux system, such as Ubuntu, run
```sh
sudo apt-get install libcr-dev mpich2 mpich2-doc
```
To install LAPACK, download the latest tar ball and unpack it.  Copy
`make.inc.example` to `make.inc`.  Run the following commands
```sh
make blaslib
make
cd lapacke
make
```

Additionally, [python3][python3] and [gnuplot][gnuplot] are needed for 
running supporting scripts.  [Tecplot][tecplot] can also be used for
visualizing the results.

<a name="compiling"></a>
## 3. Compiling

### 3.1. Clone this project to a local folder.
```sh
git clone https://github.com/jp-sheehan/squ1d
```

### 3.2. Modify `solverdefs.h` to the desired compiler/code settingns

There are two types of compiler settings in the file `solverdefs.h` that must
be set before compiling. First the type of compiler/libraries which will be used
must be selected. If Intel/MKL are used set the COMPILER variable to 1, for
g++/LAPACK set the COMPILER variable to 0. There are also two sections
below the compiler setting which must be commented or uncommented to correspond
with this choice. Finally the type of simulation is chosen. Typically
the Full simulation model will be used (0), but a number of other options exist
which are used for verification and testing simulations. (Note: There are also
parameters in `constants.h` which must be changed for some of the verification
simulations.)

### 3.3. Make sure environment variables are correctly set. Use "export" command. (NOTE: Usually export does not work)

Before compiling and running the code make sure that the environment
variables are set correctly to link to the libraries used. In particular, make sure
to set LD LIBRARY PATH. Below is an example of a command used to set this
variable when compiling with MKL:
```sh
export LD LIBRARY PATH=/opt/intel/mkl/lib/intel64:$(LD LIBRARY PATH)
```
Note that this path depends on the machine being used as well as the type
of compiler and/or libraries used. You can check if the path is implemented
correctly by using the following command:
```sh
echo $LD LIBRARY PATH
```

### 3.4. Change the relevant makefile as needed to link to the correct libraries
 
Three makefiles are included in `makefiles\` to serve as examples for compiling the code.
The default makefile (`makefile_default`) is an example of a makefile used to compile
with Intel MKL Compilers. The makefile `makefile LP` is an example of
a makefile used to compile using g++ and LAPACK libraries. The makefile
`makefile_FLUX` is an example of a makefile used to compile on a supercomputer
such as the University of Michigan Flux Supercomputer. The makefile
`makefile_STAMPEDE` is an example of a makefile used to compile on the
Stampede supercomputer which is part of XSEDE. These makefiles should be
altered so that the libraries are correctly linked to their location on the local
computer.

### 3.5. Type the command: `make -f \name of makefile`

To run any of the makefiles other than the default makefile use the command
below:
```sh
make -f \name of makefile"
```

### 3.6. The code should compile without errors
 
After typing the make command of your choice, the code should compile
without errors. A typical error is that the libraries are linked incorrectly, check
if the makefile is using the correct locations for the libraries. Another error
could be that you are using the incorrect compiler, the compiler should notify
you if it does not recognize the compiler.

### 3.7. Recompile with changes to the code
 
If compiled correctly, the code should be ready to run. The only reason to
recompile the code after this would be if changes are made to the source code,
if you select different definitions in `solverdefs.h`, or change the constants in
`constants.h`. When recompiling, the code will generally look for files that are
changed and only recompile those, but sometimes the makefile does not register
changes in particular if changes are made to the header files. In this case use the
command `make clean` to clean the compiled objects then recompile using the
makefile. This is good practice if for some reason a compilation error can not
be resolved or if the code is giving unexpectedly strange results. The command
`make cleandat` can also be used to remove .dat and .out files in the current
directory.


<a name="input"></a>
## 4. Input Files

### 4.1. SolverType.inp

The `SolverType.inp` file specifies the type of solver used by the code. The
variables used, what they represent, and the options are summarized in the
table below. All variables use MKS units. The text describing the variable is
just a placeholder used to help arrange the input file. The parameters must be
input in this format in the order listed.



### 4.2. SolverInput.inp

All solver input files are arranged in the way described in this section. The
text are just placeholders used to help arrange the input file. The parameters
must be input in this format in the given order. Sections are used to specify
different parameters. First there is a section where each of the groups of particles
or species is declared. Next a section is included for the collisions modeled,
specifically there is a Coulomb and neutral collisions section. After this the
boundary conditions are selected. This is followed by a section for the initial
conditions. Finally a section is included for applying special regions such as
heating and particle sources. More details for each of these sections is given
below. Note that profiles can be given for many of the parameters according to
the parser. More details can be found online. For the parser the parameter x is
used as the axial direction.

<a name="running"></a>
## 5. Running the code

The code is run by typing either of the following commands, the first being for a
non-parallel version and the second for the parallel version with, for example, four processors:

1. `./PIC.exe` 
2. `mpirun -n 4 ./PIC.exe` 

The code should try to run and will be successful if the input files are included
correctly. If any of the libraries are linked incorrectly an error may appear. If
this is the case make sure that paths are declared correctly by the environment
variables.

When running the code looks for two input files `SolverType.inp` and `SolverInput.inp`
which are read to initialize the code. These input files have a very
specific layout which must be followed which is presented below. Each contains
a EOF to signify the end of the file which should catch some errors in the input,
but not all. The inputs typically include text describing the variable to
be read in which is followed by a INT, DOUBLE, or STRING read into the
code. For simulations which utilize species specific cross section data the folder
`CrossSectionData/` is also necessary with the desired cross sections.

<a name="references"></a>
## 6. References

1. F. H. Ebersohn, J. P. Sheehan, A. D. Gallimore, and J. V. Shebalin, “Quasi-one-dimensional particle-in-cell simulation of magnetic nozzles,” in International Electric Propulsion Conference, 2015, pp. IEPC–2015–357. [pdf][IEPC-2015-357]
2. F. H. Ebersohn, “Kinetic Method for Quasi-One-Dimensional Simulation of Magnetic Nozzle Plasmadynamics,” University of Michigan, 2016. [pdf][ebersohn-dissertation]


[IEPC-2015-357]:                   http://pepl.engin.umich.edu/pdf/IEPC-2015-357.pdf
[ebersohn-dissertation]:           http://pepl.engin.umich.edu/pdf/2016_Ebersohn_Thesis.pdf
[intel-libs]:                      https://software.intel.com/en-us/qualify-for-free-software/academicresearcher
[python3]:                         https://www.python.org/
[gnuplot]:                         http://www.gnuplot.info/
[tecplot]:                         http://www.tecplot.com/

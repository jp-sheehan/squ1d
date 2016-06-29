#ifndef SOLVERDEF_H

#define COMPILER 1   // Specifies type of compiler:     0 = GCC, 1 = INTEL/MKL

//  Comment below for INTEL/MKL, Uncomment fo GCC. Make sure it is consistent with COMPILER.
//#pragma GCC push_options
//#pragma GCC optimize ("-O0")
//extern "C" void dgtsv_(int *n,int *nrhs, double *dl, double *d,double *du,double *b,int *ldb,int *info);
//#pragma GCC pop_options

//  Uncomment below for INTEL/MKL, comment for GCC. Make sure it is consistent with COMPILER.
#include "mkl.h"
#include "mkl_scalapack.h"


#define SIMTYPE 0  //  Specifies simulation type:    0=Full, 1=Magnetic Mirror, 2=Collision Check, 3=Single Particle, 4=Collisional Magnetic Mirror


#endif

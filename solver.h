#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <complex.h>
#include <mpi.h>
#include "variables.h"
#include "constants.h"
#include "mesh.h"
#include "mathfunctions.h"
#include "time.h"
#include <iomanip>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <chrono>
#include "exprtk.hpp"
//#include "mkl.h"
#include "omp.h"
#include <complex.h>  //Before FFTW
//#include <fftw3.h>

#ifndef SOLVER_H
#define SOLVER_H
class solver
{
   private:

   public:
      double totalTime,deltaT;
      double coulForce;

      void updateEulerFlow(flow &,const mesh &,double, int,int);//..(CFD Project)
      void updateEulerFlow2D(flow &,const mesh &,double, int,int); //..(CFD Project)
      void timestep(const particles &, const std::vector<double> &,std::vector<double>, std::vector<int>,double,double,int, int); //..Determine timeste
      void weighDensityPIC(const std::vector<particles> &, std::vector<contnm> &,const mesh &, MPIvars &); //..Interpolate properties to cell (N,En)
      void weighContinuumPIC(const std::vector<particles> &, std::vector<contnm> &,const mesh &, MPIvars &); //..Interpolate properties to cell (N,En)
      void weighEfield(std::vector<double> &,const particles&,const fields&,const mesh &,int); //..Interpolate electric field to particles
      void weighBfield(std::vector<double> &,particles&,const fields&,const mesh &,int);  //..Interpolate magnetic field to particles
      void updatePartVel(particles &, const fields &,const mesh &,int);  //..Update particle velocities
      void updatePartPos(particles &, const mesh &, int); //..Update particle position
      void updateNeutralBackground(contnm &, const mesh &, double); //..Update Neutral Background
      void clearallParticles(particles &);  //..Clear particles
      void cleanFlow(std::vector<contnm> &,int,int); //..Set Flow properties to zero
      void clearallFlow(contnm &); //..Clear flow properties
      void poisson1D(mesh,std::vector<double> &, const std::vector<contnm> &, std::vector<boundvars> &, const std::vector<particles>&, double, double); //..Solve Poisson's equation in 1D
      void thomas1Dpoisson(mesh,std::vector<double> &, const std::vector<contnm> &);  //..Solve Poisson's equation in 1D with Thomas algorithm
      void updateAppliedEfield(std::vector<std::string>, std::vector <double> &, const mesh &); //..Update applied electric field
      void updateAppliedBfield(std::vector<std::string>, std::vector <double> &, const mesh &); //..Update applied magnetic field
      void updateAppliedPotential(std::string, std::vector <double> &, const mesh &);  //..Update applied potential
      void phitoE(std::vector<double> &, std::vector<double> &, const mesh &);  //..Calculate electric field from potential
      void checkNAN(const particles, const fields, const mesh);  //..Check particles and fields for NANS
      void particlefluxsource1D(particles &, const mesh &, std::string, std::string, std::vector<double>, double,double,double,int,int,int); //..Particle flux source
      void particlecellsource(particles &, const mesh &, std::string, std::string, std::vector<double>, double,double,int); //..Partice source at cell
      void collideParticles(std::vector<particles> &,const mesh &,contnm &,double, int);  //..Collide particles
      void coulombCollisions(std::vector<particles> &,const mesh &,std::vector<contnm> &, MPIvars &);  //..Coulomb Collisions of particles
      void seedSingleParticle(particles &,const mesh &,std::vector<double>,std::vector<double>, double, int, int); //..Seed a single particle
      void findseedproc(particles &); //..Find processor to seed to for load balancing
      void redistributeparticles(particles &,int,int); //..Find processor to seed to for load balancing

      double eqnparser(std::string, std::vector<double>);        //..Equation parser
      void randomselect(std::vector<int> &,int,int);  //..Create array of randomly selected indexes
      //void scatterParticle(particles &,double,double,double,double,double,double,int,int);
      void scatterParticle(std::vector<double> &,double,double,double);

      /*
      double errorCheck(const flow&,const mesh &); //..CFD
      void flux1D(std::vector<double> &,double,double,double,double); //..CFD
      void xflux2D(std::vector<double> &,double,double,double,double,double); //..CFD
      void yflux2D(std::vector<double> &,double,double,double,double,double); //..CFD
      void flux1DLF(std::vector<double> &, const std::vector<double> &,const std::vector<double> &,double,double,double,double,double,double,double); //..CFD
      void flux2DLF(std::vector<double> &, const std::vector<double> &,const std::vector<double> &,double,double,double,double,double,double,double,double,double); //..CFD
      double thetaMUSCL(double,double,double);  //..CFD
      double queueL(double,double,double);  //..CFD
      double queueR(double,double,double); //..CFD
      */

};
#endif

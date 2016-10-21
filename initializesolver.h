#include <iostream>
#include <algorithm>
#include <vector>
#include <map>
#include <string>
#include <regex>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <ctime>
#include <chrono>
#include <cstdio>
#include "constants.h"
#include "variables.h"
#include "mesh.h"
#include "writeoutput.h"
#include "boundary.h"
#include "exprtk.hpp"
#include "mathfunctions.h"
#include "solver.h"
#include "solverdefs.h"

#ifndef INITIALIZE_H
#define INITIALIZE_H

class initializeSolver
{
	private:

	std::string rd_inputfilename;
        std::string rd_solverfilename;
        std::string rd_solverlogfilename;
        std::string rd_meshfilename, rd_meshfileformat,rd_meshwrite;
        double loopcurrent, loopradius;  //..current loop parameters
        double loopcenter[3]; 

        int rd_meshread;      
 
	public:

        int readsolver(); //..Read solver type
	void readdata(solverVars &,mesh &, std::vector<boundvars> &, std::vector<spclvars> &); //..Read data for specific run and solver type
        void initializePIC(std::vector<particles> &,fields &, std::vector<contnm> &, mesh &,std::vector<boundvars> &,std::vector<spclvars> &,solverVars &, writeOutput &,MPIvars &);//..Intialize PIC solver
        //void initializeEuler(flow &, mesh &,boundary &,solverVars);  //..Initialize Euler solver (CFD project)
        void initializeParticles(particles &,mesh,solverVars,int); //..Initialize particles
        void initializeParticles(particles &,mesh,solverVars,int,int); //..Initialize particle from restart file
        void initializeEfield(std::vector<double> &, std::vector<double>,int,int); //..Initialize electric field
        void initializeBfield(std::vector<double> &, std::vector<double>,int,int,bool); //..Initialize magnetic field
        void initializeQ1D(std::vector<particles> &, fields &, mesh &,std::vector<boundvars> &); //..Initialize quasi-1D algorithms
        void initializephi(std::vector<double> &, const mesh &);  //..Initialized electric potential
        void initializephitoE(std::vector<double> &, std::vector<double> &, const mesh &); //..Initialized electric potential and electric field
        void initializeEperp(std::vector<double> &, const mesh &,std::vector<spclvars> &); //..Initialized electric potential and electric field
        void initializeFluid(flow &, std::vector<double>,int,int); //..Initializes fluid properties
        void initializeFlow(flow &, std::vector<double>,int,int); //..Initializes flow properties
        void initializeFlowPIC(flow &, std::vector<double>,int,int);  //..Initializes flow properties for PIC simulation
        void initializeContinuum(contnm &, const particles &, std::vector<double>,int,int);  //..Initializes continuum properties
        void initializeMPIvars(MPIvars &, std::vector<double>,int,int,int);  //..Initializes continuum properties
        void initializeNeutralBackground(contnm &, const mesh &);  //..Initializes neutral background
        void initializeCollisions(std::vector<particles> &, mesh &); //..Initialize collision parameters
        void readcrsstable(std::vector<double> &, std::vector<double> &, std::vector<double> &,std::vector<double> &, std::vector<double> &, std::string, std::string,double &,double &); //..Read cross section tables
        void initializeRestart(writeOutput &, double &, int, int); //..Initialize Restart

        double eqnparser(std::string, std::vector<double>);  //..Parses equations from text file
        std::vector<double> loop_cartesian(std::vector<double>,int,int); //..Current loop for cartesian coordinates
        std::vector<double> init_CL_Bfield(fields &, std::vector<double>); //..Initializes current loop magnetic field
        void seedSingleParticle(particles &, mesh,solverVars,int);  //..Seeds a singe particle (for debugging)
	void clearallspclv(spclvars &);  //..Clears special region variables

        std::vector<std::string> init_U;  //..Initial velocity
        std::vector<std::string> init_B;  //..Initial magnetic field
        std::vector<std::string> init_E;  //..Initial electric field
        std::vector<std::string> init_dens,init_Temp; //..Initial density and temperature
        std::string init_p;  //..Initial pressure
        std::string init_phi,init_rho; //..Initial charge density and electric potential
        std::vector<std::string> dens_dist; //..Density distribution
        std::vector<std::string> therm_dist; //..Velocity/thermal distribution
        std::vector<double> dens_pert;  //..Density perturbation
        std::string sptype,setype; //..Solver types

        std::vector<double> scharge, smass; //..Species charge and mass
        std::vector<std::string> sname;  //..Species name
        std::vector<bool> smag;  //..Species magnetization
        std::vector<int> scycIT;  //..Species subcycle iteration
   
        std::vector<std::string> scllsnname;  //..Collision Type
        std::vector<std::string> scllsntype;  //..Collision Type
        std::vector<std::string> scrsstype;  //..Cross-section type
        std::vector<double> scrsssctn;  //..Cross-section
        std::string neutdens,neuttemp;  //..Background neutral density and temperature
        std::string rhoback;  //..Background charge density

        int nsp; //..Number of species
        int nct; //..Number of collision types
        int sflag; //..Solver flag
        double areain;  //..Area at inlet
 
        unsigned long num_particles; //..Number of particles
};
#endif



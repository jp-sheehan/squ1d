#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <cstdio>
#include <numeric>
#include "constants.h"
#include "variables.h"
#include "mesh.h"
#include "solver.h"
#include "exprtk.hpp"
#include <mpi.h>

#ifndef BOUNDARY_H
#define BOUNDARY_H
class boundary
{
	private:

	public:

        void clearallbdv(boundvars &); //..Clear all boundary values	
	void seedParticles(mesh, std::vector<contnm> &, std::vector<particles> &, solverVars, std::vector<boundvars> &b_bvar); //..Seeds particles at boundary
	void cleanParticles(particles&, std::vector<double>, std::vector<boundvars> &,mesh &,int, int,int); //..Move,remove, and counts particles
	void collectParticles(particles&, std::vector<double>, int, int,int); //..
	void applyParticleBoundaryConditions(mesh, std::vector<contnm> &, std::vector<particles> &, solverVars, std::vector<boundvars> &); //..Apply particle BC
        void initializeBoundary(mesh &,std::vector<boundvars> &,std::vector<particles> &); //..Initializes boundary functions
        void initializeSpecialRegions(mesh &,std::vector<spclvars> &,std::vector<particles> &); //..Initializes special regions
        void initializeParticleCounter(std::vector<boundvars> &,std::vector<particles>); //..Initialize particle counter for those leaving domain
        void setRestartVariables(std::vector<boundvars> &,std::vector<spclvars> &, std::vector<particles>,int,int); //..Initialize particle counter for those leaving domain
        void findBoundaryCells(mesh &,std::vector<boundvars> &);  //..Find cells on boundary
        void findSpecialRegionCells(mesh &,std::vector<spclvars> &);  //..Find cells in special regions
        void applyPeriodicBoundaryConditions(mesh,flow &); //..Apply periodic BC (CFD Proj)
        void applyBoundaryConditionsEuler(mesh,flow &);  //..Apply Euler BC (CFD Proj)
        void applyContinuumBoundaryConditions(mesh,fields &,std::vector<contnm> &,std::vector<boundvars> &); //..Apply Continuum BC
        void applyBoundaryConditionsPIC(mesh,fields &,flow &); //..
        void applyEfieldBoundary(mesh &, fields &, std::vector<boundvars> &); //..Apply boundary conditions on Electric Field
	void applySpecialRegionsPIC(mesh, std::vector<contnm> &, std::vector<particles> &, solverVars, std::vector<spclvars> &, fields &); //..Apply special regions
        double beqnparser(std::string, std::vector<double>);  //..Parser for boundary condition inputs


        int nbound;
        int nsp;


};
#endif

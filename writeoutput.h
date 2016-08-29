#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <string>
#include <fstream>
#include <sstream>
#include "mesh.h"
#include "variables.h"
#include "constants.h"
#include "solver.h"


#ifndef WRITE_H
#define WRITE_H
class writeOutput 
{
	private:

	public:
        int write_val,totcont_fnum,flow_fnum,rst_fnum;  //..file numbers
        std::vector<int> part_fnum, cont_fnum, field_fnum,vdf_fnum,ps_fnum; //..File numbers

        writeOutput(int,int); //..Initialize writing
        void writepMesh(const mesh &); //..Write corner point mesh
        void writecMesh(const mesh &);  //..Write cell center mesh
        void writefMesh(std::vector<double>, std::vector<int>, int);  //..Write face mesh
        void writeParticles(particles, std::vector<int>, int,int,int,double);  //..Write particles (2D)
        void writeParticles(particles,fields,mesh,double);  //..Write particles (1D)
        void writeSingleParticle(particles,fields,mesh,double);  //..Write single particle
        void writeFields(fields,int,int);  //..Write field output
        void writeFlowEuler(mesh, flow); //..Write Euler (CFD)
        void writePICField(mesh, contnm, fields,double);  //..Write continuum and field for PIC
        void writePICField(mesh, std::vector<contnm>, fields,double);  //..Write energies for PIC
        void writeRestart(const std::vector<boundvars> &, const std::vector<spclvars> &, double,int,int,int);  //..Write restart file

        void findvdf(particles,mesh, int, double);  //..Find vdf
        void findphasespace(particles,mesh, int, double);  //..Find x-vx phasespace
        void findglobalenergy(const std::vector<particles> &, const fields &, const mesh &, double time); //..Find global energy
};
#endif

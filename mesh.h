#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <math.h>
#include "variables.h"
#include "omp.h"
#include "constants.h"

#ifndef MESH_H
#define MESH_H
class mesh
{
   private:

   public:

      // mesh(std::vector<double>,std::vector<double>,std::vector<int>,int,int,int,int);
      void generateEulermesh(double);  //..Generate Euler mesh
      void generatePICmesh(double);  //..Generate PIC mesh
      void meshPoint(std::vector<double> &,std::vector<double>,double); //..Generate points (corner points) mesh
      void meshCenter(std::vector<double> &, std::vector<double>,double); //..Generate mesh for cell centers
      void meshPointNeighbors(std::vector<int> &, std::vector<int> &); //..Generate array of point neighbors
      //void meshPointNeighbors(); //..Generate array of point neighbors
      //void meshCenterNeighbors(std::vector<int> &, std::vector<int> &); //..Generate array of neighbors for cell centers
      void meshCenterNeighbors(std::vector<int> &, std::vector<int> &); //..Generate array of neighbors for cell centers
      void meshFace(double); //..Generate mesh for cell faces
      void initconnectpartandmesh(particles &); //..Initially connect particles and mesh
      void connectpartandmesh(particles &);  //..Connect particles and mesh
      void disconnectpartandmesh(particles &);  //..Disconnect particles and mesh
      void mesharea(fields); //..Find mesh area

      int meshdims,vecdims;
      int Bloc,Eloc,philoc; // 0: cell centers; 1: cell corners
      int numghost,pdist;
      int smthcnt;
      int perflag;
      int sflag;
      int wrtflag;
      int  intscheme,mvscheme;
      int  axisflag,q1dflag;
      double deltax,hdeltax; // what are these?
      double areain;
      std::string sym;

      std::vector<double> pmesh; //..Corner points
      std::vector<double> cmesh; //..Cell centers
      std::vector<double> parea; //..Corner cross-section areas (normal=z-direction)
      std::vector<double> carea; //..Center cross-section areas (normal=z-direction)
      std::vector<double> rarea; //..Flux tube-areas (normal=r-direction)
      std::vector<double> fmesh; //..Face points
      std::vector<double> c2mesh;
      std::vector<double> p2mesh;
      std::vector<double> meshshift; //..Shifting of mesh
      std::vector<double> cellweight; //..
      std::vector<double> meshlength; //..Length of mesh
      std::vector<double> meshstart; //..Starting point of mesh
      std::vector<double> meshbd; //..
      std::vector<int> numpoints; //..Number of points
      std::vector<int> pofpneigh; //..Corner neighbors of corner 
      std::vector<int> cofpneigh; //..Center neighbors of corner
      std::vector<int> cofcneigh; //..Center neighbors of center
      std::vector<int> pofcneigh; //..Corner neighbors of center

      int nbound; //..
      int npbound; //..
      int ncbound;  //..
      std::vector<int> cboundnum;
      std::vector<int> cboundcells;
      std::vector<int> cintcells;
      std::vector<int> pboundnum;
      std::vector<int> pboundcells;
      std::vector<int> pintcells;


      int  nspcl;  //..Number of special regions
      int  nsource;  //..Number of source regions

      int mesh_size(std::vector<double>);
      std::vector<int> pneighofc(int) const; //..Point neighbors of cell center
      std::vector<int> cneighofp(int) const; //..Cell center neighbors of point
      std::vector<int> cneighofc(int) const; //..Cell center neighbors of cell center
      std::vector<int> pneighofp(int) const; //..Point neighbors of point
      std::vector<int> pneighofp(std::vector<double>,int,int,int) const; //..
      std::vector<int> fneighofc(int) const; //..Face neighbors of cell
      std::vector<int> cneighoff(int) const; //..Cell center neighbors of face
      std::vector<int> pneighoff(int) const; //..Point neighbors of face
      int incneighofc(int) const; //..Internal cell center neighbors of cell center


      std::vector<double> pt_areaweight(double,double,int) const;  //..Area weighting of point based on particle location
      std::vector<double> pt_lineweight(double,int) const; //..Line weighting of point based on particle location
      std::vector<double> pt2part_lineweight(double,std::vector<int>) const; //..Line weighting from point to particle
      std::vector<double> cell_areaweight(double,double,int) const; //..Area weighting of cell center based on particle location
      std::vector<double> cell_lineweight(double,int) const;  //..Line weighting of cell center based on particle location
      std::vector<double> cell2part_lineweight(double,std::vector<int>) const; //..Line weighting from cell center to particle
      void p2p_lineweight(double,double [],int []) const; //..Line weighting from point to particle
      void c2p_lineweight(double, double [],int []) const; //..Line weighting from cell center to particle
      void pt_lw(double [],double,int) const; //..Line weighting of point based on particle location
      void cell_lw(double [],double,int) const; //..Line weighting of cell based on particle location

      int nearp(particles const&,int) const; //..Find nearest point to particle
      int nearp1D(particles const&,int) const;  //..Find nearest point to particle (1D)
      int nearp1Dnew(particles const&,int) const;  //..Find nearest point to particle (1D)
      int nearp(std::vector<double>) const;  //..Find nearest point to data
      int nearc(particles const&,int) const; //..Find nearest cell center to particle (1D)
      int nearc1D(particles const&,int) const; //..Find nearerst cell center to particle (1D)
      int nearc1Dnew(particles const&,int) const; //..Find nearerst cell center to particle (1D)
      int nearc(std::vector<double>) const; //..Find nearest cell center to data
      int nearf(particles const&,int) const; //..Find nearset face to particle
      double dist(std::vector<double>, std::vector<double>) const; //..Distance between two points
      double dist1D(double, double) const;  //..Distance between two points (1D)
      double cellvolume(int) const;  //..Cell center volume
      double pointvolume(int) const; //..Volume around a point
      double facearea(int) const;  //..Area of Face
      double squarearea(double,double,double,double) const; //..Area of square
      std::vector<double> extravec_c2p(std::vector<double>, int); //..Extrapolate vector from center to point
      std::vector<double> extravec_p2c(std::vector<double>, int); //..Extrapolate vector from point to center
      double extrascal_c2p(std::vector<double>, int);  //..Extrapolate scalar from center to point
      double extrascal_p2c(std::vector<double>, int);  //..Extrapolate scalar from point to center

};
#endif

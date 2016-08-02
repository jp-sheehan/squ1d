#include <string>
#include <cstring>
#include <map>
#include <iostream>
#include <fstream>
#include <vector>
#include "variables.h"
#include "constants.h"
#include "mesh.h"

#ifndef MATHFUNCTIONS_H
#define MATHFUNCTIONS_H

class mathFunctions 
{
   public:
      double el_first_comp(double); //..First complete integral
      double el_second_comp(double); //..Second complete integral
      double el_cylinderical();
      std::vector<double> el_cartesian(double,double,double*,double*);
      double agm(double, double);  //..Arithmetic-Geometric Mean
      double magm(double, double); //..Modified arithmetic-geometric mean
      double mag(std::vector<double>, int); //..Magnenitude
      std::vector<double> cross(std::vector<double>,std::vector<double>, int); //..Cross product of vector 
      double cross(std::vector<double>, std::vector<double>, int, int); //..Cross product of vector (single comp)
      double erfinv(double);

      std::vector<double> pcurl2D(mesh, std::vector<double>, int);

      double ddx1Dpwp(const mesh &, const std::vector<double> &, int);   //..Derivative: 1D point with point (scalar)
      double ddx1Dpwp(const mesh &, std::vector<double>, int, int);  //..Derivative: 1D point with point (vector)
      double ddx1Dcwc(const mesh &, std::vector<double>, int); //..Derivative:  1D cell with cell (scalar) 
      double ddx1Dcwc(const mesh &, std::vector<double>, int, int); //..Derivative: 1D cell with cell (vector)
      double ddx1Dpwc(const mesh &, std::vector<double>, int); //..Derivative:  1D point with cell (scalar)
      double ddx1Dpwc(const mesh &, std::vector<double>, int,int); //..Derivative:  1D point with cell (scalar)
      double ddx1Dpwc(const mesh &, double,std::vector<double>, int); //..Derivative:  1D point with cell (scalar)
      double ddx1Dcwp(const mesh &, std::vector<double>, int); //..Derivative:  1D cell with point (point)
      double ddx1Dcwp(const mesh &, std::vector<double>, int, int); //..Derivative:  1D cell with point (vector)

      void smoothData(std::vector<double> &,int); //..Smooth Scalar Data
      void smoothData(std::vector<double> &,int,int,int); //..Smooth Vector Data

      /*double ddxpwc(mesh, std::vector<double>, int, int,int) ;
        double ddxpwc(mesh, std::vector<double>, int,int) ;
        double ddxcwp(mesh, std::vector<double>, int, int,int) ;
        double ddxcwp(mesh, std::vector<double>, int,int) ;
        double ddxp(mesh, std::vector<double>, int,int) ;
        double ddxp(mesh, std::vector<double>, int, int,int) ;
        double ddxc(mesh, std::vector<double>, int,int) ;
        double ddxc(mesh, std::vector<double>, int, int,int) ;
        */

   private:  
      double el_loopcurrent, el_loopradius, el_r,el_theta;
      double el_loopcenter[3], el_centroid[3];
};
#endif

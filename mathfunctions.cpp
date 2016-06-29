#include "mathfunctions.h"

const double el_tol = (1E-10); //..Tolerance


//..Method for current loop based on a paper by Semjon Adlaj entitled 'An Eloquent Formula for the Perimeter of an Ellipse'..//


/************************Arithmetic-Geometric Mean******************/
double mathFunctions::agm(double el_x1, double el_y1)
 {
    double el_x2, el_y2, el_tolval;
    do
     {
        el_x2 = 0.5*(el_x1+el_y1);
        el_y2 = sqrt(el_x1*el_y1);
        el_x1 = el_x2;
        el_y1 = el_y2;
        el_tolval = fabs(el_x1-el_y1);
     }
     while(el_tolval>el_tol);

     return el_x1;
 }
/*******************************************************************/

/****************Modified Arithmetic-Geometric Mean*****************/
double mathFunctions::magm(double el_x1, double el_y1)
 {
    double el_x2, el_y2, el_z2, el_z1, el_tolval;
    el_z1 = 0.0;
    do
     {
        el_x2 = 0.5*(el_x1+el_y1);
        el_y2 = el_z1 + sqrt((el_x1-el_z1)*(el_y1-el_z1));
        el_z2 = el_z1 - sqrt((el_x1-el_z1)*(el_y1-el_z1));
        el_x1 = el_x2;
        el_y1 = el_y2;
        el_z1 = el_z2;
        el_tolval = fabs(el_x1-el_y1);
     }
     while(el_tolval>el_tol);

     return el_x1;
 }
/*******************************************************************/


/******************First Complete Integral**************************/

double mathFunctions::el_first_comp(double el_k)
 {
  double el_val;
  el_val = agm(1.0,sqrt(1.0-el_k*el_k));

  return pi/(2.0*el_val);
 }

/******************Second Complete Integral************************/

double mathFunctions::el_second_comp(double el_k)
 {
  double el_val;
  el_val = magm((1.0),(1.0-el_k*el_k))/(agm(1.0,sqrt(1.0-el_k*el_k)));

  return pi/2.0*el_val;
 }

/****************************************************************/

double mathFunctions::el_cylinderical()
 {

 return 0;
 }

/***************************************************************/

//..Magnetic field due to current loop calculation..//

std::vector<double> mathFunctions::el_cartesian(double el_loopcurrent, double el_loopradius, double el_loopcenter[3],double el_centroid[3])
{

  double el_x, el_y, el_z, el_B0, el_r, el_theta;
  double el_alpha, el_beta, el_gamma, el_Q, el_kay, el_K, el_E;
  std::vector<double>  el_B(3);
  double  el_Br;	

  el_x = el_centroid[0] - el_loopcenter[0];
  el_y = el_centroid[1] - el_loopcenter[1];
  el_z = el_centroid[2] - el_loopcenter[2];

  el_B0= el_loopcurrent*mu0/(2.0*el_loopradius);
 
  el_r = sqrt(el_y*el_y + el_z*el_z);
  el_theta = atan2(el_y,el_z);

  el_alpha = el_r/el_loopradius;
  el_beta = el_x/el_loopradius;
  el_gamma = el_x/el_r;

  if (el_alpha == 1.0 && el_x == 0.0)  el_alpha = el_alpha + 1E-3;

  el_Q = (1.0+el_alpha)*(1.0+el_alpha)+el_beta*el_beta;
  el_kay = sqrt(4.0*el_alpha/el_Q);
                       
  el_K = el_first_comp(el_kay);
  el_E = el_second_comp(el_kay);

  el_B[0] = (el_B0/(pi*sqrt(el_Q)))*(el_E*(1.0-el_alpha*el_alpha-el_beta*el_beta)/(el_Q-4.0*el_alpha)+el_K) ; 
  el_Br = el_B0*el_gamma/(pi*sqrt(el_Q))*(el_E*(1.0+el_alpha*el_alpha+el_beta*el_beta)/(el_Q-4.0*el_alpha)-el_K);

  if (el_r == 0.0)  el_Br = 0.0;
	       
  el_B[1] =  el_Br*sin(el_theta);
  el_B[2] =  el_Br*cos(el_theta);
	
  return el_B;
 }

//..Magnitude of vector..//

double mathFunctions::mag(std::vector<double> vec, int s_vec_dims)
{
   int i;
   double mag = 0;

   for(i=0; i<s_vec_dims; i++) mag = mag + vec[i]*vec[i];
   return sqrt(mag);
}

//..Cross product of vector..//

std::vector<double> mathFunctions::cross(std::vector<double> vec1, std::vector<double> vec2, int s_vec_dims)
{
   int i;
   std::vector<double> c_cvec;

   for(i=0; i<s_vec_dims; i++) c_cvec.push_back(0);

   c_cvec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
   c_cvec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
   c_cvec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

   return c_cvec;
} 

//..Cross product of vector..//

double mathFunctions::cross(std::vector<double> vec1, std::vector<double> vec2, int s_vec_dims, int s_dim)
{
   int i;
   std::vector<double> c_cvec;

   for(i=0; i<s_vec_dims; i++) c_cvec.push_back(0);

   c_cvec[0] = vec1[1]*vec2[2] - vec1[2]*vec2[1];
   c_cvec[1] = vec1[2]*vec2[0] - vec1[0]*vec2[2];
   c_cvec[2] = vec1[0]*vec2[1] - vec1[1]*vec2[0];

   return c_cvec[s_dim];
}

//..Inverse error function..//

double mathFunctions::erfinv(double x)
{
   //....Approximation from "Handbook of Mathematical Function", Abramowitz and Segun...(wikiversion,check later)......//

   double a = 8*(pi-3)/(3*pi*(4-pi));
   double value,term1,term2,term3;

   //value = (0.5)*sqrt(pi)*(x + pi*pow(x,3)/12.0 + 7.0*pow(pi,2)*pow(x,5)/480.0 + 127.0*pow(pi,3)*pow(x,7)/40320.0 + 4369.0*pow(pi,4)*pow(x,9)/5806080.0 + 34807.0*pow(pi,5)*pow(x,11)/182476800.0);  
   term1 = (2.0/(pi*a)+log(1-x*x)/2.0);
   term2 = log(1-x*x)/a;
   term3 = 2.0/(pi*a)+log(1-x*x)/2.0;
   
   if(x >= 0.0)   value = sqrt( sqrt(term1*term1-term2)-term3);
   else if(x < 0.0)   value = -1.0*sqrt( sqrt(term1*term1-term2)-term3);

   return value;
}

//..Curl for a point (fix)..//

std::vector<double> mathFunctions::pcurl2D(mesh s_msh, std::vector<double> s_vector, int p_index)
{
   int i,j,k;
   std::vector<double> vec_return;
  
   double dAxdx, dAxdy;
   double dAydx, dAydy;
   double dAzdx, dAzdy;
   
   for(i=0;i<s_msh.vecdims;i++) vec_return.push_back(0.0); 

}

//..Derivative: 1D Center with Center (scalar)..//

double mathFunctions::ddx1Dcwc(const mesh &s_msh, std::vector<double> s_vector, int p_index)
{
   int i,j,k;
   double dAdx;

   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofcneigh[2*p_index+j];
   dAdx = (s_vector[s_msh.cofcneigh[2*p_index+1]]-s_vector[s_msh.cofcneigh[2*p_index]])/(s_msh.cmesh[s_msh.cofcneigh[2*p_index+1]]-s_msh.cmesh[s_msh.cofcneigh[2*p_index]]) ;

   return dAdx;
}



//..Derivative: 1D Center with Center (vector)..//

double mathFunctions::ddx1Dcwc(const mesh &s_msh, std::vector<double> s_vector,int vecelem, int p_index)
{
   int i,j,k;
   double dAdx;
  
   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofcneigh[2*p_index+j];

   dAdx = (s_vector[s_msh.vecdims*s_msh.cofcneigh[2*p_index+1]+vecelem]-s_vector[s_msh.vecdims*s_msh.cofcneigh[2*p_index]+vecelem])/(s_msh.cmesh[s_msh.cofcneigh[2*p_index+1]]-s_msh.cmesh[s_msh.cofcneigh[2*p_index]]) ;

   return dAdx;
}

//..Derivative: 1D Point with Point (scalar)..//

double mathFunctions::ddx1Dpwp(const mesh &s_msh, const std::vector<double> &s_vector, int p_index)
{
   int i,j,k;
   double dAdx;
  
   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofpneigh[2*p_index+j];

   dAdx = (s_vector[s_msh.pofpneigh[2*p_index+1]]-s_vector[s_msh.pofpneigh[2*p_index]])/(2*s_msh.deltax) ;

   return dAdx;
}


//..Derivative: 1D Point with Point (vector)..//

double mathFunctions::ddx1Dpwp(const mesh &s_msh, std::vector<double> s_vector,int vecelem, int p_index)
{
   int i,j,k;
   double dAdx;
  
   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofpneigh[2*p_index+j];
   dAdx = (s_vector[s_msh.vecdims*s_msh.pofpneigh[2*p_index+1]+vecelem]-s_vector[s_msh.vecdims*s_msh.pofpneigh[2*p_index]+vecelem])/(2*s_msh.deltax) ;

   return dAdx;
}


//..Derivative: 1D Point with Cell (scalar)..//

double mathFunctions::ddx1Dpwc(const mesh &s_msh, std::vector<double> s_vector, int p_index)
{
   int i,j,k;
   double dAdx;
  
   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*p_index+j];

   dAdx = (s_vector[s_msh.cofpneigh[2*p_index+1]]-s_vector[s_msh.cofpneigh[2*p_index]])/(s_msh.cmesh[s_msh.cofpneigh[2*p_index+1]]-s_msh.cmesh[s_msh.cofpneigh[2*p_index]]) ;

   return dAdx;
}

//..Derivative: 1D Point with Cell (scalar) efficient..//      //EFF

double mathFunctions::ddx1Dpwc(const mesh &s_msh, double dx, std::vector<double> s_vector, int p_index)
{
   return (s_vector[s_msh.cofpneigh[2*p_index+1]]-s_vector[s_msh.cofpneigh[2*p_index]])/dx ;
}


//..Derivative: 1D Point with Cell (vector)..//

double mathFunctions::ddx1Dpwc(const mesh &s_msh, std::vector<double> s_vector,int vecelem, int p_index)
{
   int i,j,k;
   std::vector<int> s_neighbors;
   double dAdx;
  
   for(i=0;i<2;i++) s_neighbors.push_back(0); 

   //s_neighbors = s_msh.cneighofp(p_index);
   for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*p_index+j];

   dAdx = (s_vector[s_msh.vecdims*s_neighbors[1]+vecelem]-s_vector[s_msh.vecdims*s_neighbors[0]+vecelem])/(s_msh.cmesh[s_neighbors[1]]-s_msh.cmesh[s_neighbors[0]]) ;

   return dAdx;
}


//..Derivative: 1D Cell with Point (scalar)..//

double mathFunctions::ddx1Dcwp(const mesh &s_msh, std::vector<double> s_vector, int p_index)
{
   int i,j,k;
   double dAdx;
  
   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*p_index+j];
   dAdx = (s_vector[s_msh.cofpneigh[2*p_index+1]]-s_vector[s_msh.cofpneigh[2*p_index]])/(s_msh.pmesh[s_msh.cofpneigh[2*p_index+1]]-s_msh.pmesh[s_msh.cofpneigh[2*p_index]]) ;

   return dAdx;
}

//..Derivative: 1D Cell with Point (vector)..//

double mathFunctions::ddx1Dcwp(const mesh &s_msh, std::vector<double> s_vector,int vecelem, int p_index)
{
   int i,j,k;
   double dAdx;

   //for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*p_index+j];
   dAdx = (s_vector[s_msh.vecdims*s_msh.pofcneigh[2*p_index+1]+vecelem]-s_vector[s_msh.vecdims*s_msh.pofcneigh[2*p_index]+vecelem])/(s_msh.pmesh[s_msh.pofcneigh[2*p_index+1]]-s_msh.pmesh[s_msh.pofcneigh[2*p_index]]) ;

   return dAdx;
}


//...Smooth scalar quantity......//

void mathFunctions::smoothData(std::vector<double> &m_scal,int m_nsmooth)
{
    int i,j,k;

    std::vector<double> m_filt;

    for(i=0;i<m_scal.size();i++) m_filt.push_back(0.0);

    m_filt[0] = m_scal[0];
    m_filt[m_scal.size()-1] = m_scal[m_scal.size()-1];

    for(k=0;k<m_nsmooth;k++)
    {
      for(i=1; i<m_scal.size()-1; i++)     m_filt[i] = 0.25*(m_scal[i-1] + 2.0*m_scal[i] + m_scal[i+1]);
      for(i=1; i<m_scal.size()-1; i++)     m_scal[i] = m_filt[i];
    }
    //std::cout << "x";
}


//...Smooth vector quantity......//

void mathFunctions::smoothData(std::vector<double> &m_vec,int m_vecdims,int m_dir,int m_nsmooth)
{
    int i,j,k;

    std::vector<double> m_filt;

    for(i=0;i<m_vec.size()/m_vecdims;i++) m_filt.push_back(0.0);

    m_filt[0] = m_vec[0+m_dir];
    m_filt[m_vec.size()-1] = m_vec[m_vecdims*(m_filt.size()-1)+m_dir];

    for(k=0;k<m_nsmooth;k++)
    {
      for(i=1; i<m_filt.size()-1; i++)     m_filt[i] = 0.25*(m_vec[m_vecdims*(i-1)+m_dir] + 2.0*m_vec[m_vecdims*i+m_dir] + m_vec[m_vecdims*(i+1)+m_dir]);
      for(i=1; i<m_filt.size()-1; i++)     m_vec[m_vecdims*i+m_dir] = m_filt[i];
    }
    //std::cout << "x";
}

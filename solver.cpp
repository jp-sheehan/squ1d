#include "solver.h"
//#include "mkl_lapacke.h"
#include "solverdefs.h"

//..Update Particle Velocities..//

void solver::updatePartVel(particles &s_part, const fields &s_EM,const mesh &s_msh,int s_mvscheme)
{
   int i,j,k;
   double coef,en_temp,fact,dt,deltaTQ1D;
   double qom = s_part.wcharge/s_part.wmass;
   double beta,eps;
   double vzp,vthp,vzm,vthm,vthg,en;
   //std::vector<double> vel,vminus,vplus,vprime,vupdate,tangent;
   //std::vector<double> int_B,int_gradB,int_E,cvec;

   long int np = s_part.pos.size()/s_msh.meshdims; // total number of particles
   long int npg  = np - s_part.gp; // number of particles that aren't ghost cells
   double sign;


   //std::cout <<std::endl <<"npg:  " << npg << std::endl;

   mathFunctions mth;

   fact = 1.0;

   dt = deltaT;  //s_part.cycIT*deltaT;
   deltaTQ1D= dt;

   //..Setup for backward half-step move scheme..//
   // I do not think this is used anymore

   if(s_mvscheme<0)
   {
      if(fabs(s_mvscheme)<3)
      {
         qom = -qom*0.5;   
      } 
      else if(fabs(s_mvscheme)==3)   
      {
         deltaTQ1D= -deltaTQ1D*0.5;
      } 
   }

   //std::cout << "\n\tAdvancing Particle Velocities...";

   //..Simple electric field..//
   if((fabs(s_mvscheme)==1) || (s_part.mag==0))  // mvscheme==1: SIMPLE
   {
      //#pragma omp parallel private(en_temp,i,j) //OMP
      //{ 
      std::vector<double> vel;
      std::vector<double> int_E;

      for(j=0; j<s_msh.vecdims; j++) 
      {
         int_E.push_back(0);
         vel.push_back(0);
      }

      //#pragma omp for //OMP
      for(i=0; i<npg; i++)
      {
         en_temp = 0.0;
         //std::cout << "\nParticle:  " << i << std::endl;
         weighEfield(int_E,s_part,s_EM,s_msh,i); // 0th or 1st order interpolation
         //std::cout << int_E[0] << "\t" << int_E[1] << "\t" << int_E[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++)  vel[j] = s_part.vel[i*s_msh.vecdims+j] + qom*dt*int_E[j];
         //std::cout << qom << "\t" << dt << "\t" << qom*dt*int_E[0] << "\n";
         for(j=0; j<s_msh.vecdims; j++) en_temp = en_temp + vel[j]*vel[j]+s_part.vel[i*s_msh.vecdims+j]*s_part.vel[i*s_msh.vecdims+j];  // Avg before and after
         //for(j=0; j<s_msh.vecdims; j++) en_temp = en_temp + vel[j]*s_part.vel[i*s_msh.vecdims+j];
         for(j=0; j<s_msh.vecdims; j++) s_part.vel[i*s_msh.vecdims+j] = vel[j];
         //std::cout << s_part.vel[i*s_msh.vecdims] << "\t" << s_part.vel[i*s_msh.vecdims+1] << "\t" << s_part.vel[i*s_msh.vecdims+2] << std::endl;
         s_part.en[i] = fabs(en_temp*(0.25)*s_part.mass);         //DC Multiply by 1/4 not 1/2 because of avging
         //s_part.en[i] = fabs(en_temp*(0.5)*s_part.mass);         
         //if(s_part.en[i]>1.6e-13) std::cout << "\n Entemp:  " << s_part.en[i] << "\t" << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t";
      } 
      //}  // OMP
   }

   //..Boris Method..//
   else if(fabs(s_mvscheme)==2) // mvscheme==2: BORIS
   { 
      //#pragma omp parallel private(i,j,en_temp) //OMP
      //{
      std::vector<double> vel,vminus,vplus,vprime,vupdate,tangent;
      std::vector<double> int_B,int_gradB,int_E,cvec;

      for(j=0; j<s_msh.vecdims; j++)  
      {
         vupdate.push_back(0);
         vminus.push_back(0);
         vplus.push_back(0);
         vprime.push_back(0);
         cvec.push_back(0);
         tangent.push_back(0);
         int_B.push_back(0);
         int_E.push_back(0);
         vel.push_back(0);
      }
      //#pragma omp for//OMP
      for(i=0; i<npg; i++)
      {
         en_temp = 0.0;
         //std::cout << "\nParticle:  " << i << std::endl;
         weighEfield(int_E,s_part,s_EM,s_msh,i);
         weighBfield(int_B,s_part,s_EM,s_msh,i);
         //std::cout << int_E[0] << "\t" << int_E[1] << "\t" << int_E[2] << std::endl;
         //std::cout << "\nintB:"  << int_B[0] << "\t" << int_B[1] << "\t" << int_B[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++)  vel[j] = s_part.vel[s_msh.vecdims*i+j];
         //std::cout << vel[0] << "\t" << vel[1] << "\t" << vel[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++)  vminus[j] = vel[j] + (0.5)*qom*dt*int_E[j];
         //std::cout << qom << "\t" << dt << "\t" << (0.5)*qom*dt*int_E[0] << "\n";
         //std::cout << vminus[0] << "\t" << vminus[1] << "\t" << vminus[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++) tangent[j] = (0.5)*int_B[j]*qom*dt;
         //std::cout << tangent[0] << "\t" << tangent[1] << "\t" << tangent[2] << std::endl;
         cvec = mth.cross(vminus, tangent, s_msh.vecdims);
         //std::cout << cvec[0] << "\t" << cvec[1] << "\t" << cvec[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++) vprime[j] = vminus[j] + cvec[j];
         //std::cout << vprime[0] << "\t" << vprime[1] << "\t" << vprime[2] << std::endl;
         coef = 2.0/(1+mth.mag(tangent,s_msh.vecdims)*mth.mag(tangent,s_msh.vecdims));
         for(j=0; j<s_msh.vecdims; j++) tangent[j] = coef*tangent[j];
         //std::cout << tangent[0] << "\t" << tangent[1] << "\t" << tangent[2] << std::endl;
         cvec = mth.cross(vprime,tangent,s_msh.vecdims);
         for(j=0; j<s_msh.vecdims; j++) vplus[j] = vminus[j] + cvec[j];
         //std::cout << vplus[0] << "\t" << vplus[1] << "\t" << vplus[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++) vel[j] = vplus[j] + (0.5)*qom*dt*int_E[j];
         //std::cout << vel[0] << "\t" << vel[1] << "\t" << vel[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++) en_temp = en_temp + vel[j]*vel[j]+s_part.vel[i*s_msh.vecdims+j]*s_part.vel[i*s_msh.vecdims+j];  // Use avg of before and after
         //for(j=0; j<s_msh.vecdims; j++) en_temp = en_temp + vel[j]*s_part.vel[i*s_msh.vecdims+j];
         for(j=0; j<s_msh.vecdims; j++) s_part.vel[i*s_msh.vecdims+j] = vel[j];
         //std::cout << s_part.vel[i*s_msh.vecdims] << "\t" << s_part.vel[i*s_msh.vecdims+1] << "\t" << s_part.vel[i*s_msh.vecdims+2] << std::endl;
         s_part.en[i] = fabs(en_temp*(0.25)*s_part.mass);       //DC  Multiply by 1/4 not 1/2 because of avging
         //s_part.en[i] = fabs(en_temp*(0.5)*s_part.mass);         
      }
      //}  // OMP
   }

   //..Quasi-1D method..//
   else if(fabs(s_mvscheme)==3) // mvscheme==3: Q1D
   {
      //#pragma omp parallel private(i,j,en_temp)  //OMP
      //{

      long double A,C,vzm,vperpm,vz,vperp;
      std::vector<double> int_B,int_gradB,int_E,vel;
      double vperp_init,vperp_m,vperp_g,vperp_p;

      //vthm = vel[2];  // theta velocity   
      //vthg = vel[2];  // Initial guess for theta velocity


      //eps = (1.E-2)*fabs(qom*deltaTQ1D*vthm*int_B[0]);  // Find order of magnitude for change in velocity and set convergence criteria
      for(j=0; j<s_msh.vecdims; j++)  
      {
         int_B.push_back(0.0);
         int_E.push_back(0.0);
         vel.push_back(0.0);
      }

      /*std::vector<double> vel,vminus,vplus,vprime,vupdate,tangent;
        std::vector<double> int_B,int_gradB,int_E,cvec;
        double vperp,vperp_init;
        double vperp_m,vperp_p,vperp_g;



        for(j=0; j<s_msh.vecdims; j++)  
        {
        vupdate.push_back(0);
        vminus.push_back(0);
        vplus.push_back(0);
        vprime.push_back(0);
        cvec.push_back(0);
        tangent.push_back(0);
        int_B.push_back(0);
        int_E.push_back(0);
        vel.push_back(0);
        }*/
      //#pragma omp for //OMP
      for(i=0; i<npg; i++)
      {

         en_temp = 0.0;
         //sign = 1.0;
         //std::cout << "\nParticle:  " << i << std::endl;
         weighEfield(int_E,s_part,s_EM,s_msh,i);
         weighBfield(int_B,s_part,s_EM,s_msh,i);
         //std::cout << int_E[0] << "\t" << int_E[1] << "\t" << int_E[2] << std::endl;
         //std::cout << "\nintB:"  << int_B[0] << "\t" << int_B[1] << "\t" << int_B[2] << std::endl;
         for(j=0; j<s_msh.vecdims; j++)  vel[j] = s_part.vel[s_msh.vecdims*i+j];

         for(j=0; j<s_msh.vecdims; j++)  vel[j] = vel[j] + (0.5)*qom*deltaTQ1D*int_E[j];

         beta = -deltaTQ1D*int_B[1]/(4*int_B[0]);

         vzm = vel[0];   // z-velocity
         vperp_init = sqrt(vel[1]*vel[1]+vel[2]*vel[2]);  // perpendicular velocity   
         vperp_m = vperp_init;
         vperp_g = vperp_init;
         vperp_p = vperp_init;
         vperp = vperp_init;

         //vthm = vel[2];  // theta velocity   
         //vthg = vel[2];  // Initial guess for theta velocity


         //eps = (1.E-2)*fabs(qom*deltaTQ1D*vthm*int_B[0]);  // Find order of magnitude for change in velocity and set convergence criteria
         eps = (1.E-2)*fabs(qom*deltaTQ1D*vperp_m*int_B[0]);  // Find order of magnitude for change in velocity and set convergence criteria

         for(j=0;j<1000;j++)  // Arbitrary choice of 1000 iteration to try to converge
         {
            //vzp = vzm-beta*(vthg+vthm)*(vthg+vthm);
            //vthp = vthm + beta*(vzp+vzm)*(vthg+vthm);

            //en_temp = vzp*vzp + vthp*vthp;

            vzp = vzm-beta*(vperp_g+vperp_m)*(vperp_g+vperp_m);
            vperp_p = vperp_m + beta*(vzp+vzm)*(vperp_g+vperp_m);

            en_temp = vzp*vzp + vperp_p*vperp_p;

            //std::cout << "vthg: " << vthg << std::endl;
            //std::cout << "vthp: "<<  vthp << std::endl;
            //std::cout << "vzp: " <<  vzp << std::endl;
            //std::cout << "j: " << j << std::endl;
            //std::cout << "en: " << en_temp << std::endl;
            //std::cout << "eps: " << eps << std::endl;


            if (fabs(vperp_g-vperp_p)/fabs(vperp_g)<eps) break;  //Check for convergence

            vthg = vthp;  // Set new guess for theta velocity
            vperp_g = vperp_p;  // Set new guess for theta velocity
         }

         vel[0] = vzp;
         //vel[2] = vthp;         
         vel[1] = vel[1]*vperp_p/vperp_init;         
         vel[2] = vel[2]*vperp_p/vperp_init;         

         /*

            vzm = vel[0];   // z-velocity
            vperpm = sqrt(vel[1]*vel[1]+vel[2]*vel[2]);  // perpendicular velocity   
            A = -int_B[1]/(int_B[0]);
            C = vzm*vzm+vperpm*vperpm;
            vz = sqrt(C)*tanh(-sqrt(C)*A*deltaTQ1D+atanh(vzm/sqrt(C)));
            vperp = sqrt(C-vz*vz); 

            vel[0] = vz;
            vel[1] = vel[1]*vperp/vperpm;         
            vel[2] = vel[2]*vperp/vperpm;      
            */

         for(j=0; j<s_msh.vecdims; j++) vel[j] = vel[j] + (0.5)*qom*deltaTQ1D*int_E[j];

         if(std::isnan(vel[0])==true || std::isnan(vel[1])==true || std::isnan(vel[2])==true) 
         {
            std::cout << "\n.....FOUND NAN......\n";
            std::cout << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t" << C << "\t";
            exit(EXIT_FAILURE);
         }

         en_temp = 0.0;

         for(j=0; j<s_msh.vecdims; j++) en_temp = en_temp + vel[j]*vel[j]+s_part.vel[i*s_msh.vecdims+j]*s_part.vel[i*s_msh.vecdims+j];  // Use avg of before and after

         for(j=0; j<s_msh.vecdims; j++) s_part.vel[i*s_msh.vecdims+j] = vel[j];
         //std::cout << s_part.vel[i*s_msh.vecdims] << "\t" << s_part.vel[i*s_msh.vecdims+1] << "\t" << s_part.vel[i*s_msh.vecdims+2] << std::endl;
         s_part.en[i] = fabs(en_temp*(0.25)*s_part.mass);       //DC  Multiply by 1/4 not 1/2 because of avging
         //if(s_part.en[i]>1.6e-15) std::cout << "\n EntempQ1D:  " << s_part.en[i] << "\t" << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t";

      }
      //}  // OMP
   }
   //std::cout << "x";

}

//..Update Particle Position..//

void solver::updatePartPos(particles &s_part, const mesh &s_msh ,int s_mvscheme)

{
   int i,j;
   double x2,y2,r2;
   double alpha,dt;

   dt = deltaT; //s_part.cycIT*deltaT;
   //std::cout << "\n\tAdvancing Particle Position...";

   //..Simple particle mover..//

   if(s_mvscheme == 1)
   {

      // #pragma omp parallel for private(j) schedule(static) //OMP
      for(i=0; i<s_part.pos.size()/s_msh.meshdims; i++)
      {
         for(j=0; j<s_msh.meshdims; j++)
         {
            s_part.pos[s_msh.meshdims*i+j] = s_part.pos[s_msh.meshdims*i+j] + s_part.vel[s_msh.vecdims*i+j]*dt;
         }
      }
   }

   //..Boris Particle mover..//
   else if(s_mvscheme == 2)    
   {
      //#pragma omp parallel for private(j) schedule(static)  //OMP
      for(i=0; i<s_part.pos.size()/s_msh.meshdims; i++)
      {
         for(j=0; j<s_msh.meshdims; j++)
         {
            //std::cout << deltaT;
            //std::cout << s_part.pos[s_mesh_dims*i+j] << "\t"  << s_part.vel[s_vec_dims*i+j] << "\t";
            s_part.pos[s_msh.meshdims*i+j] = s_part.pos[s_msh.meshdims*i+j] + s_part.vel[s_msh.vecdims*i+j]*dt;
            //std::cout << s_part.pos[s_mesh_dims*i+j] << "\n";
         }
      }
   }
   //..Quasi-1D particle mover..//
   else if(s_mvscheme == 3) 
   {
      //#pragma omp parallel for private(j) schedule(static) //OMP
      for(i=0; i<s_part.pos.size(); i++)
      {
         s_part.pos[i] = s_part.pos[i] + s_part.vel[s_msh.vecdims*i]*dt;

      }
   }
   //std::cout << "x";
}


//..Calculate time-step size..//

void solver::timestep(const particles &s_part, const std::vector<double> &s_pmesh, std::vector<double> s_meshlength, std::vector<int> s_numpoints, double s_dt, double s_cfl, int s_mesh_dims, int s_vec_dims)
{
   int i;
   double maxvel, mindx,temp;

   maxvel = *std::max_element(s_part.vel.begin(), s_part.vel.end());

   mindx = 1e10;

   for(i=0; i<s_mesh_dims; i++)
   {
      temp = s_meshlength[i]/(s_numpoints[i]-1);
      mindx = std::min(mindx,temp);
   }

   deltaT = std::min(s_dt, s_cfl*mindx/maxvel);
}

//..Clean flow, set to zero..//

void solver::cleanFlow(std::vector<contnm> &s_cont,int s_mesh_dims, int s_vec_dims)
{
   int i,j,k;

   //std::cout << "\n\tCleaning Flow...";

   for(k=0; k<s_cont[0].nsp; k++)
   {
      for(i=0; i<s_cont[k].N.size(); i++)
      {
         s_cont[k].N[i] = 0.0;
         //s_cont[k].N_all[i] = 0.0;
         //s_cont[k].Ntotal[i] = 0.0;
         s_cont[k].En[i] = 0.0;
         //s_cont[k].En_all[i] = 0.0;
         //s_cont[k].Entotal[i] = 0.0;
         s_cont[k].T[i] = 0.0;
         for(j=0; j<s_vec_dims; j++)  s_cont[k].U[s_vec_dims*i+j] = 0.0;
         //for(j=0; j<s_vec_dims; j++)  s_cont[k].U_all[s_vec_dims*i+j] = 0.0;
         //for(j=0; j<s_vec_dims; j++)  s_cont[k].Utotal[s_vec_dims*i+j] = 0.0;
         for(j=0; j<s_vec_dims; j++)  s_cont[k].Tvec[s_vec_dims*i+j] = 0.0;
         for(j=0; j<s_vec_dims; j++)  s_cont[k].Qvec[s_vec_dims*i+j] = 0.0;
         //for(j=0; j<s_vec_dims; j++)  s_cont[k].Tvec_all[s_vec_dims*i+j] = 0.0;
         //for(j=0; j<s_vec_dims; j++)  s_cont[k].Tvectotal[s_vec_dims*i+j] = 0.0;
      }
   }
   //std::cout << "x";
}

//..Clear flow properties..//

void solver::clearallFlow(contnm &s_cont)
{
   //std::cout << "\n\tClearing continuum...";

   s_cont.N.clear();
   s_cont.T.clear();
   s_cont.U.clear();
   s_cont.Tvec.clear();
   s_cont.Qvec.clear();
   s_cont.En.clear();

   //std::cout << "x";
}



//..Interpolate density to cells (N)..//

void solver::weighDensityPIC(const std::vector<particles> &s_part, std::vector<contnm> &s_cont, const  mesh &s_msh, MPIvars &s_mpiv)
{
   int i,j,k,m;
   double mok;

   cleanFlow(s_cont,s_msh.meshdims,s_msh.vecdims);

   //std::cout << "\n\tInterpolating Density...";

   if(s_msh.meshdims==1)
   {
      if(s_msh.philoc==1)
      {
         if(s_msh.intscheme == 0)  //Nearest Point Interpolation
         {
            //#pragma omp parallel private(i,j,k) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            for(k=0; k<s_part[0].nsp;k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;


               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].pt[i];
                  //#pragma omp atomic // OMP
                  s_cont[k].N[s_index] += 1;
               }
               //#pragma omp master //OMP
               //{
               if(s_msh.perflag==0)
               {
                  //std::cout << "\nN:" << s_cont[k].N[0] << "\t" << s_cont[k].N[s_msh.pmesh.size()-1] << "\n";
                  s_cont[k].N[0] += s_cont[k].N[0];
                  s_cont[k].N[s_msh.pmesh.size()-1] += s_cont[k].N[s_msh.pmesh.size()-1];
               }
               else if(s_msh.perflag==1)
               {
                  s_cont[k].N[0] += s_cont[k].N[s_msh.pmesh.size()-1];
                  s_cont[k].N[s_msh.pmesh.size()-1] = s_cont[k].N[0];
               }
               //}  //OMP
               //#pragma omp for schedule(static)   //OMP

               for(i=0; i<(s_msh.pmesh.size()); i++)
               {
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.pointvolume(i);
               }
            }
            //}  //OMP
         }
         else if(s_msh.intscheme == 1)  //Linear Interpolation
         {
            //#pragma omp parallel private(i,j,k,m) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            std::vector<int> s_neighbors;
            //std::vector<double> s_arearatio;
            std::vector<double> s_pos;
            int size = 2;
            double s_arearatio[2];

            for(j=0;j<2; j++) s_neighbors.push_back(0);
            //for(j=0;j<2; j++) s_arearatio.push_back(0.0);
            for(j=0;j<2; j++) s_pos.push_back(0.0);

            for(k=0; k<s_part[0].nsp; k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;

               //#pragma omp for schedule (static) //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].cell[i];
                  for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];
                  //s_cellvolume = s_msh.pointvolume(s_index);
                  //s_arearatio = s_msh.pt_lineweight(s_part[k].pos[i*s_msh.meshdims],s_index);
                  //s_msh.pt_lw(s_arearatio,s_part[k].pos[i*s_msh.meshdims],s_index);
                  s_msh.cell_lw(s_arearatio,s_part[k].pos[i*s_msh.meshdims],s_index);

                  for(j=0;j<size; j++) s_arearatio[j] = s_arearatio[j]/s_msh.deltax;//s_cellvolume;

                  for(j=0;j<s_neighbors.size();j++) 
                  {
                     //#pragma omp atomic // OMP
                     s_cont[k].N[s_neighbors[j]] += s_arearatio[j];
                  }
               }

               if(s_msh.perflag==0)
               {
                  //std::cout << "\nN:" << s_cont[k].N[0] << "\t" << s_cont[k].N[s_msh.pmesh.size()-1] << "\n";
                  s_cont[k].N[0] += s_cont[k].N[0];
                  s_cont[k].N[s_msh.pmesh.size()-1] += s_cont[k].N[s_msh.pmesh.size()-1];
                  //std::cout << "\nN:" << s_cont[k].N[0] << "\t" << s_cont[k].N[s_msh.pmesh.size()-1] << "\n";
               }
               else if(s_msh.perflag==1)
               {
                  s_cont[k].N[0] += s_cont[k].N[s_msh.pmesh.size()-1];
                  s_cont[k].N[s_msh.pmesh.size()-1] = s_cont[k].N[0];
               }

               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<(s_msh.pmesh.size()); i++)
               {
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.pointvolume(i);
               }
               //std::cout << "\nN:" << s_cont[k].N[0] << "\t" << s_cont[k].N[s_msh.pmesh.size()-1] << "\n";
               //#pragma omp master //OMP
               //{
               //s_cont[k].N[0] = s_part[k].pwght*s_cont[k].N[0]/s_msh.cellvolume(1);
               //s_cont[k].N[s_msh.cmesh.size()-1] = s_part[k].pwght*s_cont[k].N[s_msh.cmesh.size()-1]/s_msh.cellvolume(s_msh.cmesh.size()-2);
               //}  //OMP
            }
            //}  //OMP
         }
      }
      else
      {
         if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
         {
            //#pragma omp parallel private(i,j,k) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            for(k=0; k<s_part[0].nsp;k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;


               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].cell[i];
                  //#pragma omp atomic // OMP
                  s_cont[k].N[s_index] += 1;
               }
               //#pragma omp master //OMP
               //{
               if(s_msh.perflag==1) //CHECK
               {
                  s_cont[k].N[0] = 0.0;
                  s_cont[k].N[s_msh.cmesh.size()-1] = 0.0;
               }
               //}  //OMP
               //#pragma omp for schedule(static)   //OMP

               for(i=s_msh.numghost; i<(s_msh.cmesh.size()-s_msh.numghost); i++)
               {
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);
               }
            }
            //}  //OMP
         }
         else if(s_msh.intscheme == 1)  //Linear Interpolation
         {
            //#pragma omp parallel private(i,j,k,m) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            std::vector<int> s_neighbors;
            //std::vector<double> s_arearatio;
            std::vector<double> s_pos;
            int size = 2;
            double s_arearatio[2];

            for(j=0;j<2; j++) s_neighbors.push_back(0);
            //for(j=0;j<2; j++) s_arearatio.push_back(0.0);
            for(j=0;j<2; j++) s_pos.push_back(0.0);

            for(k=0; k<s_part[0].nsp; k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;

               //#pragma omp for schedule (static) //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].pt[i];
                  for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];
                  //s_cellvolume = s_msh.pointvolume(s_index);
                  //s_arearatio = s_msh.pt_lineweight(s_part[k].pos[i*s_msh.meshdims],s_index);
                  s_msh.pt_lw(s_arearatio,s_part[k].pos[i*s_msh.meshdims],s_index);

                  for(j=0;j<size; j++) s_arearatio[j] = s_arearatio[j]/s_msh.deltax;//s_cellvolume;

                  for(j=0;j<s_neighbors.size();j++) 
                  {
                     //#pragma omp atomic // OMP
                     s_cont[k].N[s_neighbors[j]] += s_arearatio[j];
                  }
               }

               if(s_msh.perflag==0)
               {
                  s_cont[k].N[1] += s_cont[k].N[0];
                  s_cont[k].N[s_msh.cmesh.size()-2] += s_cont[k].N[s_msh.cmesh.size()-1];
               }
               else if(s_msh.perflag==1)
               {
                  s_cont[k].N[1] += s_cont[k].N[s_msh.cmesh.size()-1];
                  s_cont[k].N[s_msh.cmesh.size()-2] += s_cont[k].N[0];
               }

               //#pragma omp for schedule(static)  //OMP
               for(i=s_msh.numghost; i<(s_msh.cmesh.size()-s_msh.numghost); i++)
               {
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);
               }
               //#pragma omp master //OMP
               //{
               s_cont[k].N[0] = s_part[k].pwght*s_cont[k].N[0]/s_msh.cellvolume(1);
               s_cont[k].N[s_msh.cmesh.size()-1] = s_part[k].pwght*s_cont[k].N[s_msh.cmesh.size()-1]/s_msh.cellvolume(s_msh.cmesh.size()-2);
               //}  //OMP
            }
            //}  //OMP
         }
      } 
   }
   else if(s_msh.meshdims==2) //FIX,GS 2D
   {
      if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
      {
         long int npg,np;
         long int s_index;
         for(k=0; k<s_part[0].nsp;k++)
         {
            np = s_part[k].pos.size()/s_msh.meshdims; 
            npg = np - s_part[k].gp;

            for(i=0; i<npg; i++)
            {
               s_index= s_part[k].cell[i];
               s_cont[k].N[s_index] = s_cont[k].N[s_index] + 1;
            }

            for(i=0; i<(s_msh.cmesh.size()/s_msh.meshdims); i++)
            {
               s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);
            }
         }
      }
      else if(s_msh.intscheme == 1)  //Linear Interpolation
      {
         long int s_index;
         double s_cellvolume;
         std::vector<int> s_neighbors;
         std::vector<double> s_arearatio;

         for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
         for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);

         for(k=0;k<s_part[0].nsp;k++)
         {
            for(i=0; i<(s_part[k].pos.size()/s_msh.meshdims); i++)
            {
               //for(j=0;j<2*s_msh.meshdims;j++) s_pos[j] = s_part[k].pos[s_msh.meshdims*i+j];
               s_index= s_part[k].pt[i];
               //s_neighbors = s_msh.cneighofp(s_index);
               for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];
               s_cellvolume = s_msh.pointvolume(s_index);
               //s_arearatio = s_msh.pt_areaweight(s_part[k].pos[i*s_msh.meshdims],s_part[k].pos[i*s_msh.meshdims+1],s_index);

               for(j=0;j<s_arearatio.size(); j++) s_arearatio[j] = s_arearatio[j]/s_cellvolume;

               for(j=0;j<s_neighbors.size();j++) 
               {
                  s_cont[k].N[s_neighbors[j]] = s_cont[k].N[s_neighbors[j]] + s_arearatio[j];
               }
            }

            for(i=0; i<(s_msh.cmesh.size()/s_msh.meshdims); i++)
            {
               s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);
            }
         }
      }
   }

   int numprocs,procid; //MPI
   //std::vector<double> total;
   //std::vector<double> N_all;

   //for(i=0;i<s_cont[0].N.size()*s_part[0].nsp;i++) N_all.push_back(0.0);

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   if(numprocs>1)
   {
      //for(i=0;i<s_cont[0].N.size()*s_part[0].nsp;i++) total.push_back(0.0);
      for(i=0;i<s_cont[0].N.size()*s_part[0].nsp;i++) s_mpiv.Ntotal[i] = 0.0;

      for(i=0;i<s_part[0].nsp;i++)  //..MPI
      {
         //for(j=0;j<s_cont[0].N.size();j++) N_all.push_back(s_cont[i].N[j]);
         for(j=0;j<s_cont[0].N.size();j++) s_mpiv.N_all[i*s_cont[0].N.size()+j] = s_cont[i].N[j];
      }

      MPI_Allreduce(&s_mpiv.N_all.front(),&s_mpiv.Ntotal.front(),s_cont[0].N.size()*s_part[0].nsp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      for(i=0;i<s_part[0].nsp;i++) //MPI
      { 
         for(j=0;j<s_cont[0].N.size();j++) s_cont[i].N[j] = s_mpiv.Ntotal[i*s_cont[0].N.size()+j];
      }
   }


   //std::cout << "x";
}


//..Interpolate continuum properties to cells (N,U,En,Tvec)..//

void solver::weighContinuumPIC(const std::vector<particles> &s_part, std::vector<contnm> &s_cont, const  mesh &s_msh, MPIvars &s_mpiv)
{
   int i,j,k,m;
   double mok;

   cleanFlow(s_cont,s_msh.meshdims,s_msh.vecdims);

   //std::cout << "\n\tInterpolating Flow...";

   /* Only 1D for now
   if(s_msh.meshdims==1)
   {
   */

      if(s_msh.philoc==1) // at corner
      {
         //std::cout << "\n CF:  " << coulForce << std::endl;
         if(s_msh.intscheme == 0 || coulForce == 1)  //Nearest Cell Interpolation
         {
            //#pragma omp parallel private(i,j,k) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            for(k=0; k<s_part[0].nsp;k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;


               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].pt[i];
                  //#pragma omp atomic // OMP
                  s_cont[k].N[s_index] += 1;
                  for(j=0;j<s_msh.vecdims;j++)
                  {
                     //#pragma omp atomic // OMP
                     s_cont[k].U[s_index*s_msh.vecdims+j] += s_part[k].vel[i*s_msh.vecdims+j];
                     s_cont[k].Tvec[s_index*s_msh.vecdims+j] += s_part[k].vel[i*s_msh.vecdims+j]*s_part[k].vel[i*s_msh.vecdims+j];
                  }
                  //s_cont[k].Qvec[s_index*s_msh.vecdims] += s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]; //Hard Code for Q1D, parallel
                  //s_cont[k].Qvec[s_index*s_msh.vecdims+1] += s_part[k].vel[i*s_msh.vecdims]*(s_part[k].vel[i*s_msh.vecdims+1]*s_part[k].vel[i*s_msh.vecdims+1]+s_part[k].vel[i*s_msh.vecdims+2]*s_part[k].vel[i*s_msh.vecdims+2]); //Hard Code for Q1D, perp
                  //s_cont[k].Qvec[s_index*s_msh.vecdims+2] += 0; //Hard Code for Q1D
                  //#pragma omp atomic // OMP
                  s_cont[k].En[s_index] += s_part[k].en[i];
               }
               //#pragma omp master //OMP
               //{
               if(s_msh.perflag==0) // non-periodic boundary
               {
                  s_cont[k].N[0] += s_cont[k].N[0];
                  s_cont[k].N[s_msh.pmesh.size()-1] += s_cont[k].N[s_msh.pmesh.size()-1];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[0+j] += s_cont[k].U[0+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+j] += s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[0+j] += s_cont[k].Tvec[0+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[0+j] += s_cont[k].Qvec[0+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j] += s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j] += s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  s_cont[k].En[0] += s_cont[k].En[0];
                  s_cont[k].En[s_msh.pmesh.size()-1] += s_cont[k].En[s_msh.pmesh.size()-1];
               }
               if(s_msh.perflag==1) // period boundary
               {
                  s_cont[k].N[0] += s_cont[k].N[s_msh.pmesh.size()-1];
                  s_cont[k].N[s_msh.pmesh.size()-1] = s_cont[k].N[0];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[0+j] += s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+j] = s_cont[k].U[0+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[0+j] += s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[0+j] += s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j] = s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j] = s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+j];
                  s_cont[k].En[0] += s_cont[k].En[s_msh.pmesh.size()-1];
                  s_cont[k].En[s_msh.pmesh.size()-1] = s_cont[k].En[0];
               }
               //}  //OMP
               //#pragma omp for schedule(static)   //OMP

               mok = s_cont[k].mass/kb; 

               for(i=0; i<(s_msh.pmesh.size()); i++)
               {
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[i*s_msh.vecdims+j] = s_cont[k].U[i*s_msh.vecdims+j]/s_cont[k].N[i];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[i*s_msh.vecdims+j] = mok*s_cont[k].Tvec[i*s_msh.vecdims+j]/s_cont[k].N[i];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[i*s_msh.vecdims+j] = s_cont[k].Qvec[i*s_msh.vecdims+j]/s_cont[k].N[i];
                  s_cont[k].En[i] = s_cont[k].En[i]/s_cont[k].N[i];  //DC
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.pointvolume(i);

                  if(s_cont[k].N[i]==0) 
                  {
                     s_cont[k].En[i] =0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].U[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Tvec[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Qvec[s_msh.vecdims*i+j] = 0.0;
                  }
               }
            }
            //}  //OMP
         }
         else if(s_msh.intscheme == 1)  //Linear Interpolation
         {
            //#pragma omp parallel private(i,j,k,m) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            std::vector<int> s_neighbors;
            //std::vector<double> s_arearatio;
            std::vector<double> s_pos;
            int size=2;
            double s_arearatio[2];

            for(j=0;j<2; j++) s_neighbors.push_back(0);
            //for(j=0;j<2; j++) s_arearatio.push_back(0.0);
            for(j=0;j<2; j++) s_pos.push_back(0.0);

            for(k=0; k<s_part[0].nsp; k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;

               //#pragma omp for schedule (static) //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].cell[i];
                  for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];
                  //s_cellvolume = s_msh.pointvolume(s_index);
                  //s_arearatio = s_msh.pt_lineweight(s_part[k].pos[i*s_msh.meshdims],s_index);
                  //s_msh.pt_lw(s_arearatio,s_part[k].pos[i*s_msh.meshdims],s_index);
                  s_msh.cell_lw(s_arearatio,s_part[k].pos[i*s_msh.meshdims],s_index);

                  for(j=0;j<size; j++) s_arearatio[j] = s_arearatio[j]/s_msh.deltax;//s_cellvolume;

                  for(j=0;j<s_neighbors.size();j++) 
                  {
                     //#pragma omp atomic // OMP
                     s_cont[k].N[s_neighbors[j]] += s_arearatio[j];

                     for(m=0;m<s_msh.vecdims;m++) 
                     {
                        //#pragma omp atomic // OMP
                        s_cont[k].U[s_neighbors[j]*s_msh.vecdims+m] += s_part[k].vel[i*s_msh.vecdims+m]*s_arearatio[j];
                        s_cont[k].Tvec[s_neighbors[j]*s_msh.vecdims+m] += s_part[k].vel[i*s_msh.vecdims+m]*s_part[k].vel[i*s_msh.vecdims+m]*s_arearatio[j];
                     }
                     //s_cont[k].Qvec[s_index*s_msh.vecdims] += s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]*s_arearatio[j]; //Hard Code for Q1D, parallel
                     //s_cont[k].Qvec[s_index*s_msh.vecdims+1] += s_part[k].vel[i*s_msh.vecdims]*(s_part[k].vel[i*s_msh.vecdims+1]*s_part[k].vel[i*s_msh.vecdims+1]+s_part[k].vel[i*s_msh.vecdims+2]*s_part[k].vel[i*s_msh.vecdims+2])*s_arearatio[j]; //Hard Code for Q1D, perp
                     //s_cont[k].Qvec[s_index*s_msh.vecdims+2] += 0; //Hard Code for Q1D

                     //#pragma omp atomic
                     s_cont[k].En[s_neighbors[j]] += s_arearatio[j]*s_part[k].en[i];
                  }
               }

               if(s_msh.perflag==0)
               {
                  s_cont[k].En[0] += s_cont[k].En[0] ;  // DC:  Want total energy?  s_cont[k].N[0];
                  s_cont[k].N[0] += s_cont[k].N[0];

                  s_cont[k].En[s_msh.pmesh.size()-1] += s_cont[k].En[s_msh.pmesh.size()-1];  // DC:  Want total energy?  s_cont[k].N[s_msh.cmesh.size()-1];
                  s_cont[k].N[s_msh.pmesh.size()-1] += s_cont[k].N[s_msh.pmesh.size()-1];

                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[m] += s_cont[k].U[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+m] += s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+m];

                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[m] += s_cont[k].Tvec[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[m] += s_cont[k].Qvec[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m] += s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m] += s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m];
               }
               else if(s_msh.perflag==1)
               {
                  s_cont[k].N[0] += s_cont[k].N[s_msh.pmesh.size()-1];
                  s_cont[k].N[s_msh.pmesh.size()-1] = s_cont[k].N[0];
                  s_cont[k].En[0] += s_cont[k].En[s_msh.pmesh.size()-1] ;  // DC:  Want total energy?  s_cont[k].N[0];
                  s_cont[k].En[s_msh.pmesh.size()-1] = s_cont[k].En[0];  // DC:  Want total energy?  s_cont[k].N[s_msh.cmesh.size()-1];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[m] += s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[s_msh.vecdims*(s_msh.pmesh.size()-1)+m] = s_cont[k].U[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[m] += s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[m] += s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m] = s_cont[k].Tvec[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.pmesh.size()-1)+m] = s_cont[k].Qvec[m];
               }

               mok = s_cont[k].mass/kb; 

               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<(s_msh.pmesh.size()); i++)
               {
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[i*s_msh.vecdims+j] = s_cont[k].U[i*s_msh.vecdims+j]/(s_cont[k].N[i]);
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[i*s_msh.vecdims+j] = mok*s_cont[k].Tvec[i*s_msh.vecdims+j]/(s_cont[k].N[i]);
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[i*s_msh.vecdims+j] = s_cont[k].Qvec[i*s_msh.vecdims+j]/(s_cont[k].N[i]);
                  s_cont[k].En[i] = s_cont[k].En[i]/s_cont[k].N[i];  // DC
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.pointvolume(i);

                  if(s_cont[k].N[i]==0.0) 
                  {
                     s_cont[k].En[i] =0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].U[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Tvec[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Qvec[s_msh.vecdims*i+j] = 0.0;
                  }
               }
               //#pragma omp master //OMP
               //{

               //if(s_cont[k].N[0] != 0.0)  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[0+j] = s_cont[k].U[0+j]/(s_cont[k].N[0]);
               //if(s_cont[k].N[s_msh.cmesh.size()-1] != 0.0) for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[(s_msh.cmesh.size()-1)*s_msh.vecdims+j] = s_cont[k].U[(s_msh.cmesh.size()-1)*s_msh.vecdims+j]/(s_cont[k].N[s_msh.cmesh.size()-1]);

               //s_cont[k].N[0] = s_part[k].pwght*s_cont[k].N[0]/s_msh.cellvolume(1);
               //s_cont[k].N[s_msh.cmesh.size()-1] = s_part[k].pwght*s_cont[k].N[s_msh.cmesh.size()-1]/s_msh.cellvolume(s_msh.cmesh.size()-2);
               //}  //OMP
            }
            //}  //OMP
         }
      }
      else // at edge
      {
         if(s_msh.intscheme == 0 || coulForce == 1)  //Nearest Cell Interpolation
         {
            //#pragma omp parallel private(i,j,k) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            for(k=0; k<s_part[0].nsp;k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;


               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].cell[i];
                  //#pragma omp atomic // OMP
                  s_cont[k].N[s_index] += 1;
                  for(j=0;j<s_msh.vecdims;j++)
                  {
                     //#pragma omp atomic // OMP
                     s_cont[k].U[s_index*s_msh.vecdims+j] += s_part[k].vel[i*s_msh.vecdims+j];
                     s_cont[k].Tvec[s_index*s_msh.vecdims+j] += s_part[k].vel[i*s_msh.vecdims+j]*s_part[k].vel[i*s_msh.vecdims+j];
                  }
                  //s_cont[k].Qvec[s_index*s_msh.vecdims] += s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]; //Hard Code for Q1D, parallel
                  //s_cont[k].Qvec[s_index*s_msh.vecdims+1] += s_part[k].vel[i*s_msh.vecdims]*(s_part[k].vel[i*s_msh.vecdims+1]*s_part[k].vel[i*s_msh.vecdims+1]+s_part[k].vel[i*s_msh.vecdims+2]*s_part[k].vel[i*s_msh.vecdims+2]); //Hard Code for Q1D, perp
                  //s_cont[k].Qvec[s_index*s_msh.vecdims+2] += 0; //Hard Code for Q1D
                  //#pragma omp atomic // OMP
                  s_cont[k].En[s_index] += s_part[k].en[i];
               }
               //#pragma omp master //OMP
               //{
               if(s_msh.perflag==1)
               {
                  s_cont[k].N[0] = 0.0;
                  s_cont[k].N[s_msh.cmesh.size()-1] = 0.0;
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[0+j] =0.0;
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[s_msh.vecdims*(s_msh.cmesh.size()-1)+j] = 0.0;
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[0+j] =0.0;
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[0+j] =0.0;
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.cmesh.size()-1)+j] = 0.0;
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.cmesh.size()-1)+j] = 0.0;
                  s_cont[k].En[0] = 0.0;
                  s_cont[k].En[s_msh.cmesh.size()-1] = 0.0;
               }
               //}  //OMP
               //#pragma omp for schedule(static)   //OMP

               mok = s_cont[k].mass/kb; 

               for(i=s_msh.numghost; i<(s_msh.cmesh.size()-s_msh.numghost); i++)
               {
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[i*s_msh.vecdims+j] = s_cont[k].U[i*s_msh.vecdims+j]/s_cont[k].N[i];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[i*s_msh.vecdims+j] = mok*s_cont[k].Tvec[i*s_msh.vecdims+j]/s_cont[k].N[i];
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[i*s_msh.vecdims+j] = mok*s_cont[k].Qvec[i*s_msh.vecdims+j]/s_cont[k].N[i];
                  s_cont[k].En[i] = s_cont[k].En[i]/s_cont[k].N[i];  //DC
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);

                  if(s_cont[k].N[i]==0) 
                  {
                     s_cont[k].En[i] =0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].U[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Tvec[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Qvec[s_msh.vecdims*i+j] = 0.0;
                  }
               }
            }
            //}  //OMP
         }
         else if(s_msh.intscheme == 1)  //Linear Interpolation
         {
            //#pragma omp parallel private(i,j,k,m) //OMP
            //{
            long int npg,np;
            long int s_index;
            double s_cellvolume;

            std::vector<int> s_neighbors;
            //std::vector<double> s_arearatio;
            std::vector<double> s_pos;
            int size=2;
            double s_arearatio[2];

            for(j=0;j<2; j++) s_neighbors.push_back(0);
            //for(j=0;j<2; j++) s_arearatio.push_back(0.0);
            for(j=0;j<2; j++) s_pos.push_back(0.0);

            for(k=0; k<s_part[0].nsp; k++)
            {
               np = s_part[k].pos.size()/s_msh.meshdims; 
               npg = np - s_part[k].gp;

               //#pragma omp for schedule (static) //OMP
               for(i=0; i<npg; i++)
               {
                  s_index= s_part[k].pt[i];
                  for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];
                  //s_cellvolume = s_msh.pointvolume(s_index);
                  //s_arearatio = s_msh.pt_lineweight(s_part[k].pos[i*s_msh.meshdims],s_index);
                  s_msh.pt_lw(s_arearatio,s_part[k].pos[i*s_msh.meshdims],s_index);

                  for(j=0;j<size; j++) s_arearatio[j] = s_arearatio[j]/s_msh.deltax;//s_cellvolume;

                  for(j=0;j<s_neighbors.size();j++) 
                  {
                     //#pragma omp atomic // OMP
                     s_cont[k].N[s_neighbors[j]] += s_arearatio[j];

                     for(m=0;m<s_msh.vecdims;m++) 
                     {
                        //#pragma omp atomic // OMP
                        s_cont[k].U[s_neighbors[j]*s_msh.vecdims+m] += s_part[k].vel[i*s_msh.vecdims+m]*s_arearatio[j];
                        s_cont[k].Tvec[s_neighbors[j]*s_msh.vecdims+m] += s_part[k].vel[i*s_msh.vecdims+m]*s_part[k].vel[i*s_msh.vecdims+m]*s_arearatio[j];
                     }
                     //s_cont[k].Qvec[s_index*s_msh.vecdims] += s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]*s_part[k].vel[i*s_msh.vecdims]*s_arearatio[j]; //Hard Code for Q1D, parallel
                     //s_cont[k].Qvec[s_index*s_msh.vecdims+1] += s_part[k].vel[i*s_msh.vecdims]*(s_part[k].vel[i*s_msh.vecdims+1]*s_part[k].vel[i*s_msh.vecdims+1]+s_part[k].vel[i*s_msh.vecdims+2]*s_part[k].vel[i*s_msh.vecdims+2])*s_arearatio[j]; //Hard Code for Q1D, perp
                     //s_cont[k].Qvec[s_index*s_msh.vecdims+2] += 0; //Hard Code for Q1D

                     //#pragma omp atomic
                     s_cont[k].En[s_neighbors[j]] += s_arearatio[j]*s_part[k].en[i];
                  }
               }

               if(s_msh.perflag==0)
               {
                  s_cont[k].En[1] += s_cont[k].En[0] ;  // DC:  Want total energy?  s_cont[k].N[0];
                  s_cont[k].N[1] += s_cont[k].N[0];

                  s_cont[k].En[s_msh.cmesh.size()-2] += s_cont[k].En[s_msh.cmesh.size()-1];  // DC:  Want total energy?  s_cont[k].N[s_msh.cmesh.size()-1];
                  s_cont[k].N[s_msh.cmesh.size()-2] += s_cont[k].N[s_msh.cmesh.size()-1];

                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[s_msh.vecdims+m] += s_cont[k].U[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[s_msh.vecdims*(s_msh.cmesh.size()-2)+m] += s_cont[k].U[s_msh.vecdims*(s_msh.cmesh.size()-1)+m];

                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[s_msh.vecdims+m] += s_cont[k].Tvec[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[s_msh.vecdims+m] += s_cont[k].Qvec[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.cmesh.size()-2)+m] += s_cont[k].Tvec[s_msh.vecdims*(s_msh.cmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.cmesh.size()-2)+m] += s_cont[k].Qvec[s_msh.vecdims*(s_msh.cmesh.size()-1)+m];
               }
               else if(s_msh.perflag==1)
               {
                  s_cont[k].N[1] += s_cont[k].N[s_msh.cmesh.size()-1];
                  s_cont[k].N[s_msh.cmesh.size()-2] += s_cont[k].N[0];
                  s_cont[k].En[1] += s_cont[k].En[s_msh.cmesh.size()-1] ;  // DC:  Want total energy?  s_cont[k].N[0];
                  s_cont[k].En[s_msh.cmesh.size()-2] += s_cont[k].En[0];  // DC:  Want total energy?  s_cont[k].N[s_msh.cmesh.size()-1];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[s_msh.vecdims+m] += s_cont[k].U[s_msh.vecdims*(s_msh.cmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].U[s_msh.vecdims*(s_msh.cmesh.size()-2)+m] += s_cont[k].U[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[s_msh.vecdims+m] += s_cont[k].Tvec[s_msh.vecdims*(s_msh.cmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[s_msh.vecdims+m] += s_cont[k].Qvec[s_msh.vecdims*(s_msh.cmesh.size()-1)+m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Tvec[s_msh.vecdims*(s_msh.cmesh.size()-2)+m] += s_cont[k].Tvec[m];
                  for(m=0;m<s_msh.vecdims;m++) s_cont[k].Qvec[s_msh.vecdims*(s_msh.cmesh.size()-2)+m] += s_cont[k].Qvec[m];
               }

               mok = s_cont[k].mass/kb; 

               //#pragma omp for schedule(static)  //OMP
               for(i=0; i<(s_msh.cmesh.size()); i++)
               {
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[i*s_msh.vecdims+j] = s_cont[k].U[i*s_msh.vecdims+j]/(s_cont[k].N[i]);
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Tvec[i*s_msh.vecdims+j] = mok*s_cont[k].Tvec[i*s_msh.vecdims+j]/(s_cont[k].N[i]);
                  for(j=0;j<s_msh.vecdims;j++) s_cont[k].Qvec[i*s_msh.vecdims+j] = s_cont[k].Qvec[i*s_msh.vecdims+j]/(s_cont[k].N[i]);
                  s_cont[k].En[i] = s_cont[k].En[i]/s_cont[k].N[i];  // DC
                  s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);

                  if(s_cont[k].N[i]==0.0) 
                  {
                     s_cont[k].En[i] =0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].U[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Tvec[s_msh.vecdims*i+j] = 0.0;
                     for(j=0;j<s_msh.vecdims; j++) s_cont[k].Qvec[s_msh.vecdims*i+j] = 0.0;
                  }
               }
               //#pragma omp master //OMP
               //{

               //if(s_cont[k].N[0] != 0.0)  for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[0+j] = s_cont[k].U[0+j]/(s_cont[k].N[0]);
               //if(s_cont[k].N[s_msh.cmesh.size()-1] != 0.0) for(j=0;j<s_msh.vecdims;j++) s_cont[k].U[(s_msh.cmesh.size()-1)*s_msh.vecdims+j] = s_cont[k].U[(s_msh.cmesh.size()-1)*s_msh.vecdims+j]/(s_cont[k].N[s_msh.cmesh.size()-1]);

               //s_cont[k].N[0] = s_part[k].pwght*s_cont[k].N[0]/s_msh.cellvolume(1);
               //s_cont[k].N[s_msh.cmesh.size()-1] = s_part[k].pwght*s_cont[k].N[s_msh.cmesh.size()-1]/s_msh.cellvolume(s_msh.cmesh.size()-2);
               //}  //OMP
            }
            //}  //OMP
         }
      }
      /*
   }
   else if(s_msh.meshdims==2)  //FIX,GS 2D
   {
      if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
      {
         long int npg,np;
         long int s_index;
         for(k=0; k<s_part[0].nsp;k++)
         {
            np = s_part[k].pos.size()/s_msh.meshdims; 
            npg = np - s_part[k].gp;

            for(i=0; i<npg; i++)
            {
               s_index= s_part[k].cell[i];
               s_cont[k].N[s_index] = s_cont[k].N[s_index] + 1;
            }

            for(i=0; i<(s_msh.cmesh.size()/s_msh.meshdims); i++)
            {
               s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);
            }
         }
      }
      else if(s_msh.intscheme == 1)  //Linear Interpolation
      {
         long int s_index;
         double s_cellvolume;
         std::vector<int> s_neighbors;
         std::vector<double> s_arearatio;

         for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
         for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);

         for(k=0;k<s_part[0].nsp;k++)
         {
            for(i=0; i<(s_part[k].pos.size()/s_msh.meshdims); i++)
            {
               s_index= s_part[k].pt[i];
               for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];
               s_cellvolume = s_msh.pointvolume(s_index);
               //s_arearatio = s_msh.pt_areaweight(s_part[k].pos[i*s_msh.meshdims],s_part[k].pos[i*s_msh.meshdims+1],s_index);

               for(j=0;j<s_arearatio.size(); j++) s_arearatio[j] = s_arearatio[j]/s_cellvolume;

               for(j=0;j<s_neighbors.size();j++) 
               {
                  s_cont[k].N[s_neighbors[j]] = s_cont[k].N[s_neighbors[j]] + s_arearatio[j];
               }
            }

            for(i=0; i<(s_msh.cmesh.size()/s_msh.meshdims); i++)
            {
               s_cont[k].N[i] = s_part[k].pwght*s_cont[k].N[i]/s_msh.cellvolume(i);
            }
         }
      }
   }
   */


   int numprocs,procid; //MPI

   /*std::vector<double> Ntotal;
     std::vector<double> Tvectotal;
     std::vector<double> Utotal;
     std::vector<double> Entotal;
     */

   /*std::vector<double> N_all;
     std::vector<double> Tvec_all;
     std::vector<double> U_all;
     std::vector<double> En_all;
     */

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   if(numprocs>1)
   {
      /*for(i=0;i<s_cont[0].N.size()*s_part[0].nsp;i++) Ntotal.push_back(0.0);
        for(i=0;i<s_cont[0].Tvec.size()*s_part[0].nsp;i++) Tvectotal.push_back(0.0);
        for(i=0;i<s_cont[0].U.size()*s_part[0].nsp;i++) Utotal.push_back(0.0);
        for(i=0;i<s_cont[0].En.size()*s_part[0].nsp;i++) Entotal.push_back(0.0);
        */
      for(i=0;i<s_cont[0].N.size()*s_part[0].nsp;i++) s_mpiv.Ntotal[i] = 0.0;
      for(i=0;i<s_cont[0].Tvec.size()*s_part[0].nsp;i++) s_mpiv.Tvectotal[i] = (0.0);
      for(i=0;i<s_cont[0].Qvec.size()*s_part[0].nsp;i++) s_mpiv.Qvectotal[i] = (0.0);
      for(i=0;i<s_cont[0].U.size()*s_part[0].nsp;i++) s_mpiv.Utotal[i] = (0.0);
      for(i=0;i<s_cont[0].En.size()*s_part[0].nsp;i++) s_mpiv.Entotal[i] = (0.0);

      for(i=0;i<s_part[0].nsp;i++)  //..MPI
      {
         //for(j=0;j<s_cont[0].N.size();j++) N_all.push_back(s_cont[i].N[j]);
         //for(j=0;j<s_cont[0].En.size();j++) En_all.push_back(s_cont[i].N[j]*s_cont[i].En[j]);
         for(j=0;j<s_cont[0].N.size();j++) s_mpiv.N_all[i*s_cont[0].N.size()+j] = s_cont[i].N[j];
         for(j=0;j<s_cont[0].En.size();j++) s_mpiv.En_all[i*s_cont[0].En.size()+j] = s_cont[i].En[j];

         for(j=0;j<s_cont[0].N.size();j++) 
         {
            //for(k=0;k<s_msh.vecdims;k++)  Tvec_all.push_back(s_cont[i].N[j]*s_cont[i].Tvec[j*s_msh.vecdims+k]);
            for(k=0;k<s_msh.vecdims;k++)  s_mpiv.Tvec_all[s_cont[0].Tvec.size()*i+(s_msh.vecdims*j+k)] = s_cont[i].N[j]*s_cont[i].Tvec[j*s_msh.vecdims+k];
            for(k=0;k<s_msh.vecdims;k++)  s_mpiv.Qvec_all[s_cont[0].Qvec.size()*i+(s_msh.vecdims*j+k)] = s_cont[i].N[j]*s_cont[i].Qvec[j*s_msh.vecdims+k];
         }
         for(j=0;j<s_cont[0].N.size();j++) 
         {
            //for(k=0;k<s_msh.vecdims;k++)  U_all.push_back(s_cont[i].N[j]*s_cont[i].U[j*s_msh.vecdims+k]);
            for(k=0;k<s_msh.vecdims;k++)  s_mpiv.U_all[s_cont[0].U.size()*i+(s_msh.vecdims*j+k)] = s_cont[i].N[j]*s_cont[i].U[j*s_msh.vecdims+k];
         } 
      }

      MPI_Allreduce(&s_mpiv.N_all.front(),&s_mpiv.Ntotal.front(),s_cont[0].N.size()*s_part[0].nsp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&s_mpiv.Tvec_all.front(),&s_mpiv.Tvectotal.front(),s_cont[0].Tvec.size()*s_part[0].nsp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      //MPI_Allreduce(&s_mpiv.Qvec_all.front(),&s_mpiv.Qvectotal.front(),s_cont[0].Qvec.size()*s_part[0].nsp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&s_mpiv.U_all.front(),&s_mpiv.Utotal.front(),s_cont[0].U.size()*s_part[0].nsp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
      MPI_Allreduce(&s_mpiv.En_all.front(),&s_mpiv.Entotal.front(),s_cont[0].En.size()*s_part[0].nsp, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

      for(i=0;i<s_part[0].nsp;i++) //MPI
      { 
         for(j=0;j<s_cont[0].N.size();j++) s_cont[i].N[j] = s_mpiv.Ntotal[i*s_cont[0].N.size()+j];
         for(j=0;j<s_cont[0].En.size();j++) s_cont[i].En[j] = s_mpiv.Entotal[i*s_cont[0].En.size()+j];
         for(j=0;j<s_cont[0].Tvec.size();j++) s_cont[i].Tvec[j] = s_mpiv.Tvectotal[i*s_cont[0].Tvec.size()+j];
         for(j=0;j<s_cont[0].Qvec.size();j++) s_cont[i].Qvec[j] = s_mpiv.Qvectotal[i*s_cont[0].Qvec.size()+j];
         for(j=0;j<s_cont[0].U.size();j++) s_cont[i].U[j] = s_mpiv.Utotal[i*s_cont[0].U.size()+j];
      }
      for(i=0;i<s_part[0].nsp;i++) //MPI
      {
         for(j=0;j<s_cont[0].N.size();j++)
         {   
            if(s_cont[i].N[j]!=0)  
            {
               s_cont[i].En[j] = s_cont[i].En[j]/s_cont[i].N[j];       
               for(k=0;k<s_msh.vecdims;k++)  s_cont[i].Tvec[j*s_msh.vecdims+k] = s_cont[i].Tvec[j*s_msh.vecdims+k]/s_cont[i].N[j];
               for(k=0;k<s_msh.vecdims;k++)  s_cont[i].Qvec[j*s_msh.vecdims+k] = s_cont[i].Qvec[j*s_msh.vecdims+k]/s_cont[i].N[j];
               for(k=0;k<s_msh.vecdims;k++)  s_cont[i].U[j*s_msh.vecdims+k] = s_cont[i].U[j*s_msh.vecdims+k]/s_cont[i].N[j];
            }
            else 
            {  
               s_cont[i].En[j] = 0.0;
               for(k=0;k<s_msh.vecdims;k++)  s_cont[i].Tvec[j*s_msh.vecdims+k] = 0.0;
               for(k=0;k<s_msh.vecdims;k++)  s_cont[i].Qvec[j*s_msh.vecdims+k] = 0.0;
               for(k=0;k<s_msh.vecdims;k++)  s_cont[i].U[j*s_msh.vecdims+k] = 0.0;
            }
         }       
      } 
   }

   for(k=0;k<s_part[0].nsp;k++) //..Temperature
   {
      for(i=0;i<s_cont[0].N.size();i++)
      {   
         //s_cont[k].T[i] = s_cont[k].En[i];
         s_cont[k].T[i] = 0.0;
         for(j=0;j<s_msh.vecdims;j++) s_cont[k].T[i] += 0.5*kb*s_cont[k].Tvec[i*s_msh.vecdims+j] -  0.5*s_cont[k].mass*(s_cont[k].U[s_msh.vecdims*i+j]*s_cont[k].U[s_msh.vecdims*i+j]);
         s_cont[k].T[i] = s_cont[k].T[i]*2.0/(3.0*kb); 
         if(s_cont[k].T[i]<1e-10) s_cont[k].T[i] = 1e-10; 
      }
   }

   //std::cout << "x";
}

//..Interpolate electric field to particles..//

void solver::weighEfield(std::vector<double> &int_E, const particles &s_part, const fields &s_EM,const mesh &s_msh,int s_pindex)
{
   int i,j,s_index;
   double arearatio,pos;
   int size = s_msh.meshdims*2;
   /*std::vector<int> s_neighbors;
     std::vector<double> s_arearatio;

     for(i=0;i<(2*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
     for(i=0;i<(2*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);*/

   int s_neighbors[2];
   double s_arearatio[2];

   /* only 1D mesh
   if(s_msh.meshdims==1)
   {
   */
      if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
      {
         if(s_msh.Eloc==0)
         {
            s_index= s_part.cell[s_pindex];
         }
         else if(s_msh.Eloc==1)
         {
            s_index= s_part.pt[s_pindex];
         }

         for(i=0; i<s_msh.vecdims; i++)  int_E[i] = s_EM.E[s_msh.vecdims*s_index+i];
      }
      if(s_msh.intscheme == 1)  //Linear Interpolation
      {
         for(i=0;i<s_msh.vecdims;i++) int_E[i] = 0.0;

         if(s_msh.Eloc==0)
         {
            s_index= s_part.pt[s_pindex];
            for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];

            //s_arearatio = s_msh.cell2part_lineweight(s_part.pos[s_pindex],s_neighbors);
            pos = s_part.pos[s_pindex];
            s_msh.c2p_lineweight(pos,s_arearatio,s_neighbors);

            for(i=0;i<size;i++) 
            {
               if(s_neighbors[i]>=0) for(j=0;j<s_msh.vecdims;j++) int_E[j] = int_E[j] + s_EM.E[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
            }
         }
         else if(s_msh.Eloc==1)
         {
            s_index= s_part.cell[s_pindex];
            for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];

            //s_arearatio = s_msh.pt2part_lineweight(s_part.pos[s_pindex],s_neighbors);
            pos = s_part.pos[s_pindex];
            s_msh.p2p_lineweight(pos,s_arearatio,s_neighbors);

            for(i=0;i<size;i++) 
            {
               if(s_neighbors[i] >=0) for(j=0;j<s_msh.vecdims;j++) int_E[j] = int_E[j] + s_EM.E[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
            }
         }
      }
/*
   }
   else if(s_msh.meshdims==2)
   {
      if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
      {
         if(s_msh.Eloc==0) s_index = s_msh.nearc(s_part, s_pindex);
         else if(s_msh.Eloc==1) s_index = s_msh.nearp(s_part, s_pindex);
         for(i=0; i<s_msh.vecdims; i++)  int_E[i] = s_EM.E[s_msh.vecdims*s_index+i];
      }
      if(s_msh.intscheme == 1)  //Linear Interpolation
      {
         for(i=0;i<s_msh.vecdims;i++) int_E[i] = 0.0;

         if(s_msh.Eloc==0)
         { 
            double s_cellvolume;

            s_index= s_msh.nearp(s_part,s_pindex);
            for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];
            s_cellvolume = s_msh.pointvolume(s_index);
            //s_arearatio = s_msh.pt_areaweight(s_part.pos[s_pindex*s_msh.meshdims],s_part.pos[s_pindex*s_msh.meshdims+1],s_index); //2D

            for(i=0;i<size;i++) s_arearatio[i] = s_arearatio[i]/s_cellvolume;

            for(i=0;i<size;i++) 
            {
               for(j=0;j<s_msh.vecdims;j++) int_E[j] = int_E[j] + s_EM.E[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
            }
         }
         else if(s_msh.Eloc==1)
         {
            double s_cellvolume;
            std::vector<int> s_neighbors;
            std::vector<double> s_arearatio;

            for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
            for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);

            s_index= s_part.cell[s_pindex];
            for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];
            s_cellvolume = s_msh.cellvolume(s_index);
            //s_arearatio = s_msh.cell_areaweight(s_part.pos[s_pindex*s_msh.meshdims],s_part.pos[s_pindex*s_msh.meshdims+1],s_index); //2D

            for(i=0;i<size;i++) s_arearatio[i] = s_arearatio[i]/s_cellvolume;

            for(i=0;i<size;i++) 
            {
               //if(s_neighbors[i]>=0)
               //{
               for(j=0;j<s_msh.vecdims;j++) int_E[j] = int_E[j] + s_EM.E[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
               //}
            }
         }
      }
   }
   */
}

//..Interpolate magnetic field to particles..//

void solver::weighBfield(std::vector<double> &int_B, particles &s_part, const fields &s_EM,const mesh &s_msh,int s_pindex)
{
   int i,j,s_index;
   double moq = fabs(s_part.mass/s_part.charge);  //For rL
   double s_gradBcoef = 0.0;
   //std::vector<int> s_neighbors;
   //std::vector<double> s_arearatio;

   //for(i=0;i<(2*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
   //for(i=0;i<(2*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);
   int size = s_msh.meshdims*2;
   int s_neighbors[2];
   double s_arearatio[2];
   double pos;

   if(s_msh.meshdims==1)
   {

      if(s_msh.q1dflag == 0)  // Normal Interpolation
      { 
         if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
         {
            if(s_msh.Bloc==0)
            {
               s_index= s_part.cell[s_pindex];
            }
            else if(s_msh.Bloc==1)
            {
               s_index= s_part.pt[s_pindex];
            }
            for(i=0; i<s_msh.vecdims; i++)  int_B[i] = s_EM.B[s_msh.vecdims*s_index+i];
         }
         if(s_msh.intscheme == 1)  //Linear Interpolation
         {
            for(i=0;i<s_msh.vecdims;i++) int_B[i] = 0.0;

            if(s_msh.Bloc==0)
            {
               s_index= s_part.pt[s_pindex];
               //s_neighbors = s_msh.cneighofp(s_index);
               for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];

               //s_arearatio = s_msh.cell2part_lineweight(s_part.pos[s_pindex],s_neighbors);
               pos = s_part.pos[s_pindex];
               s_msh.c2p_lineweight(pos,s_arearatio,s_neighbors);

               for(i=0;i<size;i++) 
               {
                  if(s_neighbors[i]>=0) for(j=0;j<s_msh.vecdims;j++) int_B[j] = int_B[j] + s_EM.B[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
                  //else if(i==0) for(j=0;j<s_msh.vecdims;j++) int_B[j] = s_EM.B[s_msh.vecdims*s_neighbors[1]+j];
                  //else if(i==1) for(j=0;j<s_msh.vecdims;j++) int_B[j] = s_EM.B[s_msh.vecdims*s_neighbors[0]+j];
               }
            }
            else if(s_msh.Bloc==1)
            {
               s_index= s_part.cell[s_pindex];
               for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];

               //s_arearatio = s_msh.pt2part_lineweight(s_part.pos[s_pindex],s_neighbors);
               pos = s_part.pos[s_pindex];
               s_msh.p2p_lineweight(pos,s_arearatio,s_neighbors);

               for(i=0;i<size;i++) 
               {
                  if(s_neighbors[i]>=0) for(j=0;j<s_msh.vecdims;j++) int_B[j] = int_B[j] + s_EM.B[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
                  //else if(i==0) for(j=0;j<s_msh.vecdims;j++) int_B[j] = s_EM.B[s_msh.vecdims*s_neighbors[1]+j];
                  //else if(i==1) for(j=0;j<s_msh.vecdims;j++) int_B[j] = s_EM.B[s_msh.vecdims*s_neighbors[0]+j];
               }
            }
         }
      }
      else if(s_msh.q1dflag == 1)  //Q1D Interpolation
      { 
         if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
         {
            if(s_msh.Bloc==0)
            {
               s_index= s_part.cell[s_pindex];
            }
            else if(s_msh.Bloc==1)
            {
               s_index= s_part.pt[s_pindex];
            }
            int_B[0] = s_EM.B[s_msh.vecdims*s_index];
            int_B[1] = s_EM.gradBcoef[s_index];
            int_B[2] = 0.0;
         }
         if(s_msh.intscheme == 1)  //Linear Interpolation
         {
            for(i=0;i<s_msh.vecdims;i++) int_B[i] = 0.0;

            if(s_msh.Bloc==0)
            {
               s_index= s_part.pt[s_pindex];
               for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];

               //s_arearatio = s_msh.cell2part_lineweight(s_part.pos[s_pindex],s_neighbors);
               pos = s_part.pos[s_pindex];
               s_msh.c2p_lineweight(pos,s_arearatio,s_neighbors);

               for(i=0;i<size;i++) 
               {
                  if(s_neighbors[i]>=0) 
                  {
                     int_B[0] = int_B[0] + s_EM.B[s_msh.vecdims*s_neighbors[i]]*s_arearatio[i];
                     s_gradBcoef = s_gradBcoef + s_EM.gradBcoef[s_neighbors[i]]*s_arearatio[i];
                  }
               }
               int_B[1] = s_gradBcoef;
               //std::cout << std::endl << "\tBz: " << int_B[0] << "\tBr: " << int_B[1] << std::endl;
               int_B[2] = 0.0;
            }
            else if(s_msh.Bloc==1)
            {
               s_index= s_part.cell[s_pindex];
               for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];

               //s_arearatio = s_msh.pt2part_lineweight(s_part.pos[s_pindex],s_neighbors);
               pos = s_part.pos[s_pindex];
               s_msh.p2p_lineweight(pos,s_arearatio,s_neighbors);

               for(i=0;i<size;i++) 
               {
                  if(s_neighbors[i]>=0)
                  {
                     int_B[0] = int_B[0] + s_EM.B[s_msh.vecdims*s_neighbors[i]]*s_arearatio[i];
                     s_gradBcoef = s_gradBcoef + s_EM.gradBcoef[s_neighbors[i]]*s_arearatio[i];
                  }
               }
               int_B[1] = s_gradBcoef;
               //std::cout << std::endl << "\tBz: " << int_B[0] << "\tBr: " << int_B[1] << std::endl;
               int_B[2] = 0.0;
            }
         }
      }
   } 
   else if(s_msh.meshdims==2)
   {
      if(s_msh.intscheme == 0)  //Nearest Cell Interpolation
      {
         if(s_msh.Bloc==0) s_index = s_msh.nearc(s_part, s_pindex);
         if(s_msh.Bloc==1) s_index = s_msh.nearp(s_part, s_pindex);
         for(i=0; i<s_msh.vecdims; i++) int_B[i] = s_EM.B[s_msh.vecdims*s_index+i];
      }
      if(s_msh.intscheme == 1)  //Linear Interpolation
      {
         if(s_msh.Bloc==0)
         {
            double s_cellvolume;
            std::vector<int> s_neighbors;
            std::vector<double> s_arearatio;

            for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
            for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);

            s_index= s_msh.nearp(s_part,s_pindex);
            for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.cofpneigh[2*s_index+j];
            s_cellvolume = s_msh.pointvolume(s_index);
            //s_arearatio = s_msh.pt_areaweight(s_part.pos[s_pindex*s_msh.meshdims],s_part.pos[s_pindex*s_msh.meshdims+1],s_index);

            for(i=0;i<size;i++) s_arearatio[i] = s_arearatio[i]/s_cellvolume;

            for(i=0;i<size;i++) 
            {
               //if(s_neighbors[i]>=0)
               //{
               for(j=0;j<s_msh.vecdims;j++) int_B[j] = int_B[j] + s_EM.B[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
               //}
            }
         } 
         else if(s_msh.Bloc==1)
         {
            double s_cellvolume;
            std::vector<int> s_neighbors;
            std::vector<double> s_arearatio;

            for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_neighbors.push_back(0);
            for(i=0;i<((s_msh.meshdims)*(s_msh.meshdims)); i++) s_arearatio.push_back(0.0);

            s_index= s_part.cell[s_pindex];
            for(j=0;j<2*s_msh.meshdims;j++) s_neighbors[j] = s_msh.pofcneigh[2*s_index+j];
            s_cellvolume = s_msh.cellvolume(s_index);
            s_arearatio = s_msh.cell_areaweight(s_part.pos[s_pindex*s_msh.meshdims],s_part.pos[s_pindex*s_msh.meshdims+1],s_index);

            for(i=0;i<size;i++) s_arearatio[i] = s_arearatio[i]/s_cellvolume;

            for(i=0;i<size;i++) 
            {
               //if(s_neighbors[i]>=0)
               //{
               for(j=0;j<s_msh.vecdims;j++) int_B[j] = int_B[j] + s_EM.B[s_msh.vecdims*s_neighbors[i]+j]*s_arearatio[i];
               //}
            }
         }
      }
   }
}

//..Clears all particles..//

void solver::clearallParticles(particles &s_part)
{
   //std::cout << "\n\tClearing all particles...";
   s_part.pos.clear();
   s_part.vel.clear();
   s_part.en.clear();
   s_part.pt.clear();
   s_part.cell.clear();
   //std::cout << "x";
}

//..Solve 1D Poisson's equation..//

void solver::poisson1D(mesh m_msh, std::vector<double> &m_phi, const std::vector<contnm> &m_cont, std::vector<boundvars> &s_bdv, const std::vector<particles> &s_part, double s_dt, double s_time)
{
   int i,j,k;
   int neumind;
   const double ieps0 = 1.0/eps0;
   const double tol = 1e-10;
   int numprocs,procid; //MPI

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI


#if COMPILER==0
   int n,nrhs,ldb,info,nm,nw; 
#elif COMPILER==1
   MKL_INT n,nrhs,ldb,info,nm,nw;
#endif

   std::vector<double> m_rho; 
   //std::vector<double> m_sigma; 
   std::vector<double> m_E; 
   std::vector<double> meshpt;
   std::vector<double> area;
   bool zeroflag;
   int numneum=0;
   int s_nsp = m_cont[0].nsp; 

   for(i=0;i<m_msh.meshdims;i++) meshpt.push_back(0.0);

   for(i=0;i<2;i++)
   {
      if(s_bdv[i].boundtype[3*s_nsp] == "NEUMANN")  numneum += 1;
   }

   if(m_msh.meshdims==1) //..Interpolating areas for surface charge located at cell interface 
   {
      area.push_back(m_msh.parea[s_bdv[0].pboundcells[0]]);
      area.push_back(m_msh.parea[s_bdv[1].pboundcells[0]]);
   }

   //std::cout << std::endl << area[0] << "\t" << area[1] << std::endl;


   if(m_msh.philoc==1)   //..For non-staggerred grids, i.e. phi on cell corners  //CHECK ALL THESE BCS
   {
      n = m_msh.pintcells.size();
      nrhs = 1;
      ldb = n;

      //std::cout << "\n\tCalculating Potential...";

      if(m_msh.perflag==0)     // Solve Poisson's Equation for Dirichlet Bounday Conditions
      {
         if(numneum==0)  // Standard Dirichlet Boundary Conditions
         {
            double dl[n-1], du[n-1],d[n],b[n];
            double J_cond,J_cond_tot,J_tot;

            for(i=0;i<(n-1);i++)  //Set upper and lower diagonal
            {
               dl[i] = 1.0;
               du[i] = 1.0;
            } 

            for(i=0;i<n;i++) d[i] = -2.0;  //Set diagonal

            for(i=0;i<m_msh.pmesh.size();i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);

            for(i=1;i<m_cont[0].nsp;i++)
            {
               for(j=0;j<m_msh.pmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];  //Collect charge density
            } 

            if(s_bdv[0].wallcap==1) m_phi[0] = m_phi[n+1]-s_bdv[0].sigma*area[0]/s_bdv[0].bdphi_val;    //Calculate the potentials based on capacitances and collected charge
            else if(s_bdv[1].wallcap==1) m_phi[n+1] = m_phi[0]+s_bdv[1].sigma*area[1]/s_bdv[1].bdphi_val;


            if(s_bdv[0].wallvolt==1) m_phi[0] = s_bdv[0].bdphi_val*sin(2*pi*s_bdv[0].Vfreq*s_time);
            else if(s_bdv[1].wallvolt==1) m_phi[n+1] = s_bdv[0].bdphi_val*sin(2*pi*s_bdv[0].Vfreq*s_time);

            //std::cout << std::endl << m_phi[0] << "\t" << s_bdv[0].sigma << "\t" << s_bdv[0].bdphi_val;


            b[0] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
            for(i=1;i<(n-1);i++) b[i] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[i+1]));  //Set RHS
            b[n-1] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS

            // matrix solver
#if COMPILER==0
            dgtsv_(&n,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
            info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
            //pdptsv(n,nrhs,d,dl,0,b,ldb);  //Solve matrix equation
#endif

            for(i=0;i<n;i++) m_phi[i+1] = b[i];

            //if(s_bdv[0].wallcap==1) s_bdv[1].sigma = eps0*(m_phi[n+1]-m_phi[n])/(m_msh.deltax)-0.5*m_rho[n+1]*m_msh.deltax;  //Calculate charge densities for other wall 
            //else if(s_bdv[1].wallcap==1) s_bdv[0].sigma = eps0*(m_phi[0]-m_phi[1])/(m_msh.deltax)-0.5*m_rho[0]*m_msh.deltax;
            s_bdv[0].Ebd = (m_phi[0]-m_phi[1])/(m_msh.deltax)-0.5*ieps0*m_rho[0]*m_msh.deltax;   //Calculate boundary electric fields
            s_bdv[1].Ebd = -((m_phi[n+1]-m_phi[n])/(m_msh.deltax)-0.5*ieps0*m_rho[n+1]*m_msh.deltax);  //
            //std::cout << "\t" << s_bdv[0].Ebd << "\t" << s_bdv[1].Ebd << "\t" ;

         }
         else if(numneum==2)    // Double Neumann Boundary Conditions....Double check all of this
         {
            nw = n;
            ldb = nw; //DC

            double dl[nw-1], du[nw-1],d[nw],b[nw];

            for(i=0;i<(nw-1);i++)  //Set upper and lower diagonals
            {
               dl[i] = 1.0;
               du[i] = 1.0;
            } 

            for(i=0;i<nw-1;i++) d[i] = -2.0;   // Set diagonal
            d[nw-1] = -1.0;  // Negative of final equation for gradient

            for(i=0;i<m_msh.pmesh.size();i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);
            //for(i=0;i<s_bdv[0].nbound;i++) m_sigma.push_back(m_cont[0].wcharge*s_bdv[i].partcount[0]/area[i]); //Calculating surface charge density

            for(i=1;i<m_cont[0].nsp;i++)
            {
               for(j=0;j<m_msh.pmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];  //Collecting charge density
               //for(j=0;j<s_bdv[0].nbound;j++) m_sigma[j] += (m_cont[i].wcharge*s_bdv[j].partcount[i]/area[j]);  //Collecting surface charge density
            }

            //std::cout << " sigma: " << m_sigma[0] << "\t" << m_sigma[1] << ".." ;

            for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[i];   //Setting left mesh point

            m_phi[0] = 0.0;   //Set wall cell as reference potential
            if(s_bdv[0].wallfloat==1)  m_phi[1] = m_phi[0]-(m_msh.deltax*ieps0)*(s_bdv[0].sigma+0.5*m_rho[0]*m_msh.deltax); //Calculate potential at i=1 based on gradient (with charge collection)
            else  m_phi[1] = m_phi[0]-(m_msh.deltax*s_bdv[0].bdphi_val+0.5*m_rho[0]*m_msh.deltax*m_msh.deltax*ieps0); //Calculate potential at i=1 domain based on gradient

            b[0] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[2]))-m_phi[1];  // Modified RHS for first point, note Poisson's equation is not solved for i=1;
            for(i=1;i<(nw-1);i++) b[i] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[i+2]));  //Set RHS

            for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[m_msh.meshdims*(m_msh.pmesh.size()-1)+i]; //Setting right mesh point

            if(s_bdv[1].wallfloat==1)  b[nw-1] = -(m_msh.deltax*ieps0)*(s_bdv[1].sigma+0.5*m_rho[nw-1]*m_msh.deltax); // Neumann BC equation and surface charge density (invert sign for convenience)
            else  b[nw-1] = (m_msh.deltax*s_bdv[1].bdphi_val+0.5*m_rho[nw-1]*m_msh.deltax*m_msh.deltax*ieps0); //Neumann BC only

#if COMPILER==0
            dgtsv_(&nw,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
            info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,nw,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
#endif

            for(i=0;i<nw;i++) m_phi[i+2] = b[i];

            //std::cout << std::endl << m_phi[0] << std::endl;

         }
         else // Single NEUMANN BC
         {
            nw = n+1;
            ldb = nw;
            double J_cond,J_cond_tot,J_tot;
            J_cond = 0.0;
            J_cond_tot = 0.0; 
            J_tot = 0.0;

            if(s_bdv[0].boundtype[3*s_nsp] == "NEUMANN") neumind = 0;  //Set left as NEUMANN BC
            else if(s_bdv[1].boundtype[3*s_nsp] == "NEUMANN") neumind = 1;  //Set right as NEUMANN BC

            double dl[nw-1], du[nw-1],d[nw],b[nw];

            for(i=0;i<(nw-1);i++)  //Set Upper and lower diagonals
            {
               dl[i] = 1.0;
               du[i] = 1.0;
            } 

            for(i=0;i<nw;i++) d[i] = -2.0;  //Set diagonal

            if(neumind==0) d[0] = -1.0; 
            else if(neumind==1) d[nw-1] = -1.0; 

            for(i=0;i<m_msh.pmesh.size();i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);
            //for(i=0;i<s_bdv[0].nbound;i++) m_sigma.push_back(m_cont[0].wcharge*s_bdv[i].partcount[0]/area[i]); //Calculating surface charge density

            for(i=1;i<m_cont[0].nsp;i++)
            {
               for(j=0;j<m_msh.pmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];  //Collecting charge density
               //for(j=0;j<s_bdv[0].nbound;j++) m_sigma[j] += (m_cont[i].wcharge*s_bdv[j].partcount[i]/area[j]);  //Collecting surface charge density
            }

            //std::cout << "Neumann: "  << neumind << " sigma: " << m_sigma[0] << "\t" << m_sigma[1] << ".." ;
            //std::cout << "Neumann2: "  << neumind << " sigma: " << s_bdv[0].sigma << "\t" << s_bdv[1].sigma << ".." ;

            if(s_bdv[0].wallcur==1) J_tot = s_bdv[0].bdphi_val*sin(2*pi*s_bdv[0].Jfreq*s_time);
            else if(s_bdv[1].wallcur==1) J_tot = s_bdv[1].bdphi_val*sin(2*pi*s_bdv[1].Jfreq*s_time);

            //if(s_bdv[0].wallcur==1) J_tot = s_bdv[0].bdphi_val/(2*pi*s_bdv[0].Jfreq)*(cos(2*pi*s_bdv[0].Jfreq*(s_time))-cos(2*pi*s_bdv[0].Jfreq*(s_time+s_dt)));
            //else if(s_bdv[1].wallcur==1) J_tot = s_bdv[1].bdphi_val/(2*pi*s_bdv[1].Jfreq)*(cos(2*pi*s_bdv[1].Jfreq*(s_time))-cos(2*pi*s_bdv[1].Jfreq*(s_time+s_dt)));

            //std::cout << "\tJ_tot: " << J_tot << "\t";

            if(neumind==0)
            {
               for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[i];  //Set left mesh point

               if(s_bdv[0].wallcur==1)  s_bdv[0].sigma += J_tot*s_dt;  // For conducting wall with driven current

               //if(s_bdv[0].wallfloat==1)  b[0] = -(m_msh.cmesh[1]-m_msh.cmesh[0])*(ieps0*m_sigma[0]+s_bdv[0].bdphi_val);  // For conducting wall include surface charge and input
               if(s_bdv[0].wallfloat==1)  b[0] = -(m_msh.deltax*ieps0)*(s_bdv[0].sigma+0.5*m_rho[0]*m_msh.deltax);  // For conducting wall include surface charge
               else if(s_bdv[0].wallcur==1)  b[0] = -(m_msh.deltax*ieps0)*(s_bdv[0].sigma+0.5*m_rho[0]*m_msh.deltax);  // For conducting wall include surface charge and currents
               else  b[0] = -(m_msh.deltax)*(s_bdv[0].bdphi_val);  // Only input
               //b[0] = -(m_msh.pmesh[1]-m_msh.pmesh[0])*ieps0*m_sigma[0];

               for(i=0;i<(nw-2);i++) b[i+1] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[i+1]));  //Set RHS
               b[nw-1] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[nw-1])-m_phi[nw]);  // Set last term RHS
            }
            else if(neumind==1)
            {
               if(s_bdv[1].wallcur==1)  s_bdv[1].sigma += J_tot*s_dt;  // For conducting wall include surface charge and input

               for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[(m_msh.pmesh.size()-1)*m_msh.meshdims+i];  //Set right mesh point.

               b[0] = (-(m_msh.deltax*m_msh.deltax)*ieps0*(m_rho[1])-m_phi[0]);   //Set first term RHS
               for(i=1;i<(nw-1);i++) b[i] = (-m_msh.deltax*m_msh.deltax*ieps0*(m_rho[i+1]));  //Set RHS
               //b[nw-1] = -(m_msh.cmesh[nw]-m_msh.cmesh[nw-1])*ieps0*m_sigma[0];
               //if(s_bdv[1].wallfloat==1)  b[nw-1] = (m_msh.cmesh[nw]-m_msh.cmesh[nw-1])*(s_bdv[1].bdphi_val-ieps0*m_sigma[1]);  //For conducting wall include surface charge and input
               if(s_bdv[1].wallfloat==1)  b[nw-1] = -(m_msh.deltax*ieps0)*(s_bdv[1].sigma+0.5*m_rho[nw]*m_msh.deltax);  //For conducting wall include surface charge 
               else if(s_bdv[1].wallcur==1)  b[nw-1] = -(m_msh.deltax*ieps0)*(s_bdv[1].sigma + 0.5*m_rho[nw]*m_msh.deltax);  //For conducting wall include surface charge and currents
               else  b[nw-1] = -0.5*m_rho[nw]*m_msh.deltax*m_msh.deltax*ieps0;  // Only input
               //else  b[nw-1] = -(m_msh.deltax)*(s_bdv[1].bdphi_val);  // Only input
            }

            //std::cout << std::setprecision(5);
            /*std::cout << std::endl << "\nrho:     " << std::endl; 
              for(j=0;j<(nw+1);j++) std::cout << m_rho[j]  << "\t";
              std::cout << std::endl << "\nbd:     " << std::endl; 
              for(j=0;j<(nw);j++) std::cout << b[j]  << "\t" << d[j] << "\t";
              std::cout << std::endl << "\nd's:     " << std::endl; 
              for(j=0;j<(nw-1);j++) std::cout << dl[j]  << "\t" << du[j] << "\t";*/

#if COMPILER==0
            dgtsv_(&nw,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
            info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,nw,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
#endif

            //std::cout << std::endl << "ba:    "  << std::endl;
            //for(j=0;j<(nw);j++) std::cout << b[j]  << "\t";

            if(neumind==0) for(i=0;i<nw;i++) m_phi[i] = b[i];
            else if(neumind==1) for(i=0;i<nw;i++) m_phi[i+1] = b[i];

            //if(neumind==0) s_bdv[1].sigma = eps0*(m_phi[nw]-m_phi[nw-1])/(m_msh.deltax)-0.5*m_rho[nw]*m_msh.deltax;  //Calculate charge densities for other wall 
            //else if(neumind==1) s_bdv[0].sigma = eps0*(m_phi[0]-m_phi[1])/(m_msh.deltax)-0.5*m_rho[0]*m_msh.deltax;

            s_bdv[0].Ebd = (m_phi[0]-m_phi[1])/(m_msh.deltax)-0.5*m_rho[0]*m_msh.deltax*ieps0;   //Calculate boundary electric fields
            s_bdv[1].Ebd = -((m_phi[nw]-m_phi[nw-1])/(m_msh.deltax)-0.5*m_rho[nw]*m_msh.deltax*ieps0);  //

            //std::cout << "\t" << s_bdv[0].sigma << "\t" << s_bdv[1].sigma << "\t" ;
            //std::cout << "\t" << s_bdv[0].Ebd << "\t" << s_bdv[1].Ebd << "\t" ;
         }

      }
      else if(m_msh.perflag==1)   // Solve Poisson's Equation for Periodic Boundary
      {
         nm = n;
         ldb = nm;
         double dl[nm-1], du[nm-1],d[nm],b[nm];
         //std::vector<double> b;
         double bdb[nm];

         //std::cout << std::endl << "n:  "<< n << std::endl;

         for(i=0;i<(nm-1);i++)
         {
            dl[i] = 1.0;
            du[i] = 1.0;
         } 

         for(i=0;i<(nm);i++) d[i] = -2.0;

         for(i=0;i<(m_msh.pmesh.size());i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);

         for(i=1;i<m_cont[0].nsp;i++)
         {
            for(j=0;j<m_msh.pmesh.size();j++) m_rho[j] += m_cont[i].charge*m_cont[i].N[j];
         }

         //std::cout << std::endl << "RHO:     " << std::endl; 
         //for(j=0;j<m_msh.pmesh.size();j++) std::cout << m_rho[j] << "\t";
         //std::cout << std::endl; 

         b[0] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[1])-m_phi[0]);
         for(i=1;i<(nm-1);i++) b[i] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[i+1]));
         b[nm-1] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[nm])-m_phi[nm+1]);

         //for(j=0;j<nm;j++) std::cout << b[j] << "\t";

#if COMPILER==0
         dgtsv_(&nm,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
         info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,nm,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
         /*auto tdg1 = std::chrono::high_resolution_clock::now();
           for(k=0;k<1000;k++) 
           {
           b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
           for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
           b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
           info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
           }
           auto tdg2 = std::chrono::high_resolution_clock::now();

           auto tdt1 = std::chrono::high_resolution_clock::now();
           for(k=0;k<1000;k++)
           {
           b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
           for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
           b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
           info = LAPACKE_dptsv(LAPACK_COL_MAJOR,n,nrhs,d,dl,b,ldb);  //Solve matrix equation
           }
           auto tdt2 = std::chrono::high_resolution_clock::now();

           std::cout << "\n\nAnnnnnnd dgtsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdg2-tdg1).count() << " ms\n\n";
           std::cout << "\n\nAnnnnnnd dptsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdt2-tdt1).count() << " ms\n\n";
           */
#endif

         //std::cout << std::setprecision(20);
         //std::cout << std::endl << "bdb:     " << std::endl; 
         //for(j=0;j<(nm);j++) std::cout << bdb[j]  << "\t";
         //std::cout << std::endl << "b:    "  << std::endl;
         //for(j=0;j<(nm);j++) std::cout << b[j]  << "\t";

         for(i=0;i<(nm);i++)
         {  
            m_phi[i+1] = b[i];
            if(std::isnan(b[i])==true) 
            {
               std::cout << "\n.....FOUND NAN......\n";
               exit(EXIT_FAILURE);
            }
         }
         //m_phi[0] = 0.0;
         //m_phi[nm] = 0.0;

         //for(i=0;i<(n+2);i++)  m_phi[i] = m_phi[i] - m_phi[n];

         //std::cout << std::endl << "phi:     " << std::endl; 
         //for(i=0;i<(n+2);i++) std::cout << m_phi[i] << "\t";
      }
      else
      {
         std::cout << "\n.....Boundary not supported in Poisson......\n";
         exit(EXIT_FAILURE);
      }
   }
   else      //Solve for staggered mesh (i.e. phi on cell centers)  //CHECK ALL THESE BCS
   {
      n = m_msh.cintcells.size();
      nrhs = 1;
      ldb = n;

      //std::cout << "\n\tCalculating Potential...";

      if(m_msh.perflag==0)     // Solve Poisson's Equation for Dirichlet Bounday Conditions
      {
         if(numneum==0)  // Standard Dirichlet Boundary Conditions
         {
            double dl[n-1], du[n-1],d[n],b[n];
            double J_cond,J_cond_tot,J_tot;

            for(i=0;i<(n-1);i++)  //Set upper and lower diagonal
            {
               dl[i] = 1.0;
               du[i] = 1.0;
            } 

            for(i=0;i<n;i++) d[i] = -2.0;  //Set diagonal
            d[0] = -3.0;  //Set first cell 
            d[n-1] = -3.0;  //Set last cell

            for(i=0;i<m_msh.cmesh.size();i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);

            for(i=1;i<m_cont[0].nsp;i++)
            {
               for(j=0;j<m_msh.cmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];  //Collect charge density
            } 


            /*if(s_bdv[1].wallcur==1 || s_bdv[0].wallcur==1)
              {
              for(j=0;j<s_part[0].nsp;j++)  //CHG?   Does this below in Dirichlet??
              {
              for(i=0;i<s_part[j].pos.size();i++) //Calculate Conducting Current 1D
              {
              J_cond += s_part[j].wcharge*s_part[j].vel[m_msh.vecdims*i]/m_msh.area[s_part[j].cell[i]];
              }
              }

              MPI_Allreduce(&J_cond,&J_cond_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  //CHG?
              J_cond = J_cond_tot;
              }*/


            /*if(s_bdv[1].wallcur==1)  J_tot = J_cond + s_bdv[1].bdphi_val*cos(2*pi*s_bdv[1].Jfreq*s_time); //Calculate Current
              else if(s_bdv[0].wallcur==1)  J_tot = J_cond + s_bdv[0].bdphi_val*cos(2*pi*s_bdv[0].Jfreq*s_time); 


              if(s_bdv[1].wallcur==1)  m_phi[n+1] = m_phi[0] + (J_tot)*s_dt/eps0; //Calculate potential at i=n+1 based on current
              else if(s_bdv[0].wallcur==1)  m_phi[0] = m_phi[n+1] - (J_tot)*s_dt/eps0; //Calculate potential at i=0 based on current
              */
            //std::cout << "J_cond: " << J_cond/m_msh.meshlength[0] << ".." ;

            if(s_bdv[0].wallcap==1) m_phi[0] = m_phi[n+1]+s_bdv[0].sigma*area[0]/s_bdv[0].bdphi_val;
            else if(s_bdv[1].wallcap==1) m_phi[n+1] = m_phi[0]-s_bdv[1].sigma*area[1]/s_bdv[1].bdphi_val;


            b[0] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[m_msh.numghost])-2*s_bdv[0].bdphi_val);  //Set first RHS
            for(i=1;i<(n-1);i++) b[i] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[i+1]));  //Set RHS
            b[n-1] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[n])-2*s_bdv[1].bdphi_val);  //Set last RHS

#if COMPILER==0
            dgtsv_(&n,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
            info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
            //pdptsv(n,nrhs,d,dl,0,b,ldb);  //Solve matrix equation
            /*auto tdg1 = std::chrono::high_resolution_clock::now();
              for(k=0;k<1000;k++) 
              {
              b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
              for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
              b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
              info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
              }
              auto tdg2 = std::chrono::high_resolution_clock::now();

              auto tdt1 = std::chrono::high_resolution_clock::now();
              for(k=0;k<1000;k++)
              {
              b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
              for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
              b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
              info = LAPACKE_dptsv(LAPACK_COL_MAJOR,n,nrhs,d,dl,b,ldb);  //Solve matrix equation
              }
              auto tdt2 = std::chrono::high_resolution_clock::now();

              std::cout << "\n\nAnnnnnnd dgtsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdg2-tdg1).count() << " ms\n\n";
              std::cout << "\n\nAnnnnnnd dptsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdt2-tdt1).count() << " ms\n\n";
              */
#endif

            for(i=0;i<n;i++) m_phi[i+1] = b[i];

         }
         else if(numneum==2)    // Double Neumann Boundary Conditions....Double check all of this
         {
            nw = n;
            ldb = nw; //DC

            double dl[nw-1], du[nw-1],d[nw],b[nw];

            for(i=0;i<(nw-1);i++)  //Set upper and lower diagonals
            {
               dl[i] = 1.0;
               du[i] = 1.0;
            } 

            for(i=0;i<nw-1;i++) d[i] = -2.0;   // Set diagonal
            d[nw-1] = -1.0;  // Negative of final equation for gradient

            for(i=0;i<m_msh.cmesh.size();i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);
            //for(i=0;i<s_bdv[0].nbound;i++) m_sigma.push_back(m_cont[0].wcharge*s_bdv[i].partcount[0]/area[i]); //Calculating surface charge density

            for(i=1;i<m_cont[0].nsp;i++)
            {
               for(j=0;j<m_msh.cmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];  //Collecting charge density
               //for(j=0;j<s_bdv[0].nbound;j++) m_sigma[j] += (m_cont[i].wcharge*s_bdv[j].partcount[i]/area[j]);  //Collecting surface charge density
            }

            //std::cout << " sigma: " << m_sigma[0] << "\t" << m_sigma[1] << ".." ;

            for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[i];   //Setting left mesh point

            m_phi[0] = 0.0;   //Set ghost cell as reference potential
            if(s_bdv[0].wallfloat==1)  m_phi[1] = m_phi[0]-(m_msh.cmesh[1]-m_msh.cmesh[0])*(ieps0*s_bdv[0].sigma+s_bdv[0].bdphi_val); //Calculate potential at i=1 based on gradient (with charge collection)
            else  m_phi[1] = m_phi[0]-(m_msh.cmesh[1]-m_msh.cmesh[0])*s_bdv[0].bdphi_val; //Calculate potential at i=1 domain based on gradient

            b[0] = (-(0.25)*(m_msh.cmesh[3]-m_msh.cmesh[1])*(m_msh.cmesh[3]-m_msh.cmesh[1])*ieps0*(m_rho[2]))-m_phi[1];  // Modified RHS for first point, not Poisson's equation is not solved for i=1;
            for(i=1;i<(nw-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+2]));  //Set RHS

            for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[m_msh.pmesh.size()]; //Setting right mesh point

            if(s_bdv[1].wallfloat==1)  b[nw-1] = (m_msh.cmesh[nw]-m_msh.cmesh[nw-1])*(s_bdv[1].bdphi_val-ieps0*s_bdv[1].sigma); // Neumann BC equation and surface charge density (invert sign for convenience)
            else  b[nw-1] = (m_msh.cmesh[nw]-m_msh.cmesh[nw-1])*s_bdv[1].bdphi_val; //Neumann BC only

#if COMPILER==0
            dgtsv_(&nw,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
            info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,nw,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
            /*auto tdg1 = std::chrono::high_resolution_clock::now();
              for(k=0;k<1000;k++) 
              {
              b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
              for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
              b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
              info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
              }
              auto tdg2 = std::chrono::high_resolution_clock::now();

              auto tdt1 = std::chrono::high_resolution_clock::now();
              for(k=0;k<1000;k++)
              {
              b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
              for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
              b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
              info = LAPACKE_dptsv(LAPACK_COL_MAJOR,n,nrhs,d,dl,b,ldb);  //Solve matrix equation
              }
              auto tdt2 = std::chrono::high_resolution_clock::now();

              std::cout << "\n\nAnnnnnnd dgtsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdg2-tdg1).count() << " ms\n\n";
              std::cout << "\n\nAnnnnnnd dptsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdt2-tdt1).count() << " ms\n\n";
              */
#endif

            for(i=0;i<nw;i++) m_phi[i+2] = b[i];

            //std::cout << std::endl << m_phi[0] << std::endl;

         }
         else // Single NEUMANN BC
         {
            nw = n+1;
            ldb = nw;
            double J_cond,J_cond_tot,J_tot;
            J_cond = 0.0;
            J_cond_tot = 0.0; 
            J_tot = 0.0;

            if(s_bdv[0].boundtype[3*s_nsp] == "NEUMANN") neumind = 0;  //Set left as NEUMANN BC
            else if(s_bdv[1].boundtype[3*s_nsp] == "NEUMANN") neumind = 1;  //Set right as NEUMANN BC

            double dl[nw-1], du[nw-1],d[nw],b[nw];

            for(i=0;i<(nw-1);i++)  //Set Upper and lower diagonals
            {
               dl[i] = 1.0;
               du[i] = 1.0;
            } 

            for(i=0;i<nw;i++) d[i] = -2.0;  //Set diagonal

            if(neumind==0) 
            { 
               d[0] = -1.0; 
               d[nw-1] = -3.0;
            } 
            else if(neumind==1) 
            {
               d[0] = -3.0; 
               d[nw-1] = -1.0; 
            }

            for(i=0;i<m_msh.cmesh.size();i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);
            //for(i=0;i<s_bdv[0].nbound;i++) m_sigma.push_back(m_cont[0].wcharge*s_bdv[i].partcount[0]/area[i]); //Calculating surface charge density

            for(i=1;i<m_cont[0].nsp;i++)
            {
               for(j=0;j<m_msh.cmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];  //Collecting charge density
               //for(j=0;j<s_bdv[0].nbound;j++) m_sigma[j] += (m_cont[i].wcharge*s_bdv[j].partcount[i]/area[j]);  //Collecting surface charge density
            }

            //std::cout << "Neumann: "  << neumind << " sigma: " << m_sigma[0] << "\t" << m_sigma[1] << ".." ;
            //std::cout << "Neumann2: "  << neumind << " sigma: " << s_bdv[0].sigma << "\t" << s_bdv[1].sigma << ".." ;

            /*if(s_bdv[1].wallcur==1 || s_bdv[0].wallcur==1)
              {
              for(j=0;j<s_part[0].nsp;j++)  
              {
              for(i=0;i<s_part[j].pos.size();i++) //Calculate Conducting Current 1D
              {
              J_cond += s_part[j].wcharge*s_part[j].vel[m_msh.vecdims*i]/m_msh.area[s_part[j].cell[i]];
              }
              }

              MPI_Allreduce(&J_cond,&J_cond_tot,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);  //CHG?
              J_cond = J_cond_tot;
              }

              J_cond = J_cond/m_msh.meshlength[0];*/
            //std::cout << "J_cond: " << J_cond << "\n" ;

            if(s_bdv[0].wallcur==1) J_tot = s_bdv[0].bdphi_val*sin(2*pi*s_bdv[0].Jfreq*s_time);
            else if(s_bdv[1].wallcur==1) J_tot = s_bdv[1].bdphi_val*sin(2*pi*s_bdv[1].Jfreq*s_time);

            //std::cout << "J_tot: " << J_tot << "\n";

            if(neumind==0)
            {
               for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[0];  //Set left mesh point

               if(s_bdv[0].wallcur==1)  s_bdv[0].sigma += J_tot*s_dt;  // For conducting wall with driven current

               //if(s_bdv[0].wallfloat==1)  b[0] = -(m_msh.cmesh[1]-m_msh.cmesh[0])*(ieps0*m_sigma[0]+s_bdv[0].bdphi_val);  // For conducting wall include surface charge and input
               if(s_bdv[0].wallfloat==1)  b[0] = -(m_msh.deltax)*(ieps0*s_bdv[0].sigma);  // For conducting wall include surface charge
               else if(s_bdv[0].wallcur==1)  b[0] = -(m_msh.deltax)*(ieps0*s_bdv[0].sigma);  // For conducting wall include surface charge and currents
               else  b[0] = -(m_msh.deltax)*(s_bdv[0].bdphi_val);  // Only input
               //b[0] = -(m_msh.cmesh[1]-m_msh.cmesh[0])*ieps0*m_sigma[0];

               for(i=0;i<(nw-1);i++) b[i+1] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
               b[nw-1] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[nw-1])-2*s_bdv[1].bdphi_val);  // Set last term RHS
            }
            else if(neumind==1)
            {
               if(s_bdv[1].wallcur==1)  s_bdv[1].sigma += J_tot*s_dt;  // For conducting wall include surface charge and input

               for(i=0;i<m_msh.meshdims;i++) meshpt[i] = m_msh.pmesh[m_msh.pmesh.size()];  //Set right mesh point.
               b[0] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[m_msh.numghost])-2*s_bdv[0].bdphi_val);   //Set first term RHS
               for(i=1;i<(nw-1);i++) b[i] = (-(m_msh.deltax)*(m_msh.deltax)*ieps0*(m_rho[i+1]));  //Set RHS
               //b[nw-1] = -(m_msh.cmesh[nw]-m_msh.cmesh[nw-1])*ieps0*m_sigma[0];
               //if(s_bdv[1].wallfloat==1)  b[nw-1] = (m_msh.cmesh[nw]-m_msh.cmesh[nw-1])*(s_bdv[1].bdphi_val-ieps0*m_sigma[1]);  //For conducting wall include surface charge and input
               if(s_bdv[1].wallfloat==1)  b[nw-1] = (m_msh.deltax)*(-ieps0*s_bdv[1].sigma);  //For conducting wall include surface charge 
               else if(s_bdv[1].wallcur==1)  b[nw-1] = (m_msh.deltax)*(-ieps0*s_bdv[1].sigma);  //For conducting wall include surface charge and currents
               else  b[nw-1] = -(m_msh.deltax)*(s_bdv[1].bdphi_val);  // Only input
            }

            //std::cout << std::setprecision(5);
            //std::cout << std::endl << "bb:     " << std::endl; 
            //for(j=0;j<(nw);j++) std::cout << b[j]  << "\t";

#if COMPILER==0
            dgtsv_(&nw,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
            info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,nw,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
            /*auto tdg1 = std::chrono::high_resolution_clock::now();
              for(k=0;k<1000;k++) 
              {
              b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
              for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
              b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
              info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
              }
              auto tdg2 = std::chrono::high_resolution_clock::now();

              auto tdt1 = std::chrono::high_resolution_clock::now();
              for(k=0;k<1000;k++)
              {
              b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
              for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
              b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
              info = LAPACKE_dptsv(LAPACK_COL_MAJOR,n,nrhs,d,dl,b,ldb);  //Solve matrix equation
              }
              auto tdt2 = std::chrono::high_resolution_clock::now();

              std::cout << "\n\nAnnnnnnd dgtsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdg2-tdg1).count() << " ms\n\n";
              std::cout << "\n\nAnnnnnnd dptsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdt2-tdt1).count() << " ms\n\n";
              */
#endif

            //std::cout << std::endl << "ba:    "  << std::endl;
            //for(j=0;j<(nw);j++) std::cout << b[j]  << "\t";

            if(neumind==0) for(i=0;i<nw;i++) m_phi[i] = b[i];
            else if(neumind==1) for(i=0;i<nw;i++) m_phi[i+1] = b[i];

         }

      }
      else if(m_msh.perflag==1)   // Solve Poisson's Equation for Periodic Boundary
      {
         nm = n-1;
         ldb = nm;
         double dl[nm-1], du[nm-1],d[nm],b[nm];
         //std::vector<double> b;
         double bdb[nm];

         //std::cout << std::endl << "n:  "<< n << std::endl;

         for(i=0;i<(nm-1);i++)
         {
            dl[i] = 1.0;
            du[i] = 1.0;
         } 

         for(i=0;i<(nm);i++) d[i] = -2.0;

         for(i=0;i<(m_msh.cmesh.size());i++) m_rho.push_back(m_cont[0].rho_back[i] + m_cont[0].charge*m_cont[0].N[i]);

         for(i=1;i<m_cont[0].nsp;i++)
         {
            for(j=0;j<m_msh.cmesh.size();j++) m_rho[j] = m_rho[j] + m_cont[i].charge*m_cont[i].N[j];
         }

         //std::cout << std::endl << "RHO:     " << std::endl; 
         //for(j=0;j<m_msh.cmesh.size();j++) std::cout << m_rho[j] << "\t";
         //std::cout << std::endl; 

         b[0] = (-(0.25)*(m_msh.cmesh[3]-m_msh.cmesh[1])*(m_msh.cmesh[3]-m_msh.cmesh[1])*ieps0*(m_rho[m_msh.numghost+1])-m_phi[m_msh.numghost]);
         for(i=1;i<(nm-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[3+i]-m_msh.cmesh[1+i])*(m_msh.cmesh[3+i]-m_msh.cmesh[1+i])*ieps0*(m_rho[i+1+m_msh.numghost]));
         b[nm-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);


#if COMPILER==0
         dgtsv_(&nm,&nrhs,dl,d,du,b,&ldb,&info);
#elif COMPILER==1
         info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,nm,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
         /*auto tdg1 = std::chrono::high_resolution_clock::now();
           for(k=0;k<1000;k++) 
           {
           b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
           for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
           b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
           info = LAPACKE_dgtsv(LAPACK_COL_MAJOR,n,nrhs,dl,d,du,b,ldb);  //Solve matrix equation
           }
           auto tdg2 = std::chrono::high_resolution_clock::now();

           auto tdt1 = std::chrono::high_resolution_clock::now();
           for(k=0;k<1000;k++)
           {
           b[0] = (-(0.25)*(m_msh.cmesh[2]-m_msh.cmesh[0])*(m_msh.cmesh[2]-m_msh.cmesh[0])*ieps0*(m_rho[m_msh.numghost])-m_phi[0]);  //Set first RHS
           for(i=1;i<(n-1);i++) b[i] = (-(0.25)*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*(m_msh.cmesh[2+i]-m_msh.cmesh[i])*ieps0*(m_rho[i+1]));  //Set RHS
           b[n-1] = (-(0.25)*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*(m_msh.cmesh[n+1]-m_msh.cmesh[n-1])*ieps0*(m_rho[n])-m_phi[n+1]);  //Set last RHS
           info = LAPACKE_dptsv(LAPACK_COL_MAJOR,n,nrhs,d,dl,b,ldb);  //Solve matrix equation
           }
           auto tdt2 = std::chrono::high_resolution_clock::now();

           std::cout << "\n\nAnnnnnnd dgtsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdg2-tdg1).count() << " ms\n\n";
           std::cout << "\n\nAnnnnnnd dptsv: " << std::chrono::duration_cast<std::chrono::milliseconds>(tdt2-tdt1).count() << " ms\n\n";
           */
#endif

         //std::cout << std::setprecision(20);
         //std::cout << std::endl << "bdb:     " << std::endl; 
         //for(j=0;j<(nm);j++) std::cout << bdb[j]  << "\t";
         //std::cout << std::endl << "b:    "  << std::endl;
         //for(j=0;j<(nm);j++) std::cout << b[j]  << "\t";

         for(i=0;i<(nm);i++)
         {  
            m_phi[i+2] = b[i];
            if(std::isnan(b[i])==true) 
            {
               std::cout << "\n.....FOUND NAN......\n";
               exit(EXIT_FAILURE);
            }
         }
         m_phi[0] = m_phi[n];

         //for(i=0;i<(n+2);i++)  m_phi[i] = m_phi[i] - m_phi[n];

         //std::cout << std::endl << "phi:     " << std::endl; 
         //for(i=0;i<(n+2);i++) std::cout << m_phi[i] << "\t";
      }
      else
      {
         std::cout << "\n.....Boundary not supported in Poisson......\n";
         exit(EXIT_FAILURE);
      }
   }

   //std::cout << "x";
}


//..Solve 1D Poisson's equation with Thomas algorithm..//

void solver::thomas1Dpoisson(mesh m_msh, std::vector<double> &m_phi, const std::vector<contnm> &m_cont)
{
   int i,j,k;
   double m;
   std::vector<double> m_c,m_d,m_dp,m_cp,m_x;
   std::vector<double> m_x1, m_x2;
   const double ieps0 = 1.0/eps0;
   std::vector<double> m_rho; 

   for(i=0;i<m_msh.pmesh.size();i++) m_rho.push_back(m_cont[0].charge*m_cont[0].N[i]);

   for(i=1;i<m_cont[0].nsp;i++)
   {
      for(j=0;j<m_msh.pmesh.size();j++) m_rho[j] = m_cont[i].charge*m_cont[i].N[j];
   } 


   if(m_msh.perflag==0) //Dirichlet Thomas Algorithm
   {

      //std::cout << "\nDirichlet Thomas:  " << std::endl;

      m_d.push_back(-(0.25)*(m_msh.pmesh[2]-m_msh.pmesh[0])*(m_msh.pmesh[2]-m_msh.pmesh[0])*ieps0*(m_rho[1])-m_phi[0]);
      for(i=1;i<m_msh.cintcells.size()-1;i++) m_d.push_back(-(0.25)*(m_msh.pmesh[2+i]-m_msh.pmesh[i])*(m_msh.pmesh[2+i]-m_msh.pmesh[i])*ieps0*(m_rho[i+1]));
      m_d.push_back(-(0.25)*(m_msh.pmesh[m_msh.cintcells.size()+1]-m_msh.pmesh[m_msh.cintcells.size()-1])*(m_msh.pmesh[m_msh.cintcells.size()+1]-m_msh.pmesh[m_msh.cintcells.size()-1])*ieps0*(m_rho[m_msh.cintcells.size()])-m_phi[m_msh.cintcells.size()+1]);

      m_cp.push_back(-0.5);
      for(i=1;i<m_msh.cintcells.size()-1;i++) m_cp.push_back(1.0/(-2.0-m_cp[i-1]));

      m_dp.push_back(-0.5*m_d[0]);
      for(i=1;i<m_msh.cintcells.size();i++) m_dp.push_back((m_d[i]-m_dp[i-1])/(-2.0-m_cp[i-1]));   //OPTIMIZE HERE

      m_phi[m_msh.cintcells.size()]  = m_dp[m_msh.cintcells.size()-1] ; // CAUTION ::: Only works for one Ghost Cell
      for(i=(m_msh.cintcells.size()-2);i>=0;i--) m_phi[i+1]  = m_dp[i]-m_cp[i]*m_phi[i+2];

   }
   else if(m_msh.perflag==1) //Periodic Modified Thomas Algorithm:   http://www.cfm.brown.edu/people/gk/chap6/node14.html, See Paper, Think about just using above with same BC/Bias as stated in Birdsall.
   {
      //std::cout << "\nPeriodic Thomas:  " << std::endl;

      for(i=0;i<m_msh.cintcells.size();i++)
      {
         m_x1.push_back(0.0);
         m_x2.push_back(0.0);
      }

      //First Matrix Solution

      m_d.push_back(-(0.25)*(m_msh.pmesh[2]-m_msh.pmesh[0])*(m_msh.pmesh[2]-m_msh.pmesh[0])*ieps0*(m_rho[1]));
      for(i=1;i<m_msh.cintcells.size()-1;i++) m_d.push_back(-(0.25)*(m_msh.pmesh[2+i]-m_msh.pmesh[i])*(m_msh.pmesh[2+i]-m_msh.pmesh[i])*ieps0*(m_rho[i+1]));

      m_cp.push_back(-0.5);
      for(i=1;i<m_msh.cintcells.size()-2;i++) m_cp.push_back(1.0/(-2.0-m_cp[i-1]));

      m_dp.push_back(-0.5*m_d[0]);
      for(i=1;i<m_msh.cintcells.size()-1;i++) m_dp.push_back((m_d[i]-m_dp[i-1])/(-2.0-m_cp[i-1]));   //OPTIMIZE HERE

      m_x1[m_msh.cintcells.size()-1]  = m_dp[m_msh.cintcells.size()-2] ;
      for(i=(m_msh.cintcells.size()-3);i>=0;i--) m_x1[i+1]  = m_dp[i]-m_cp[i]*m_x1[i+2];

      //Second Matrix Solution

      m_d[0] = -1.0;
      for(i=1;i<m_msh.cintcells.size()-2;i++) m_d[i] = 0.0;
      m_d[m_msh.cintcells.size()-2] = -1.0;

      m_cp[0] = -0.5;
      for(i=1;i<m_msh.cintcells.size()-2;i++) m_cp[i] = 1.0/(-2.0-m_cp[i-1]);

      m_dp[0] = -0.5*m_d[0];
      for(i=1;i<m_msh.cintcells.size()-1;i++) m_dp[i] = (m_d[i]-m_dp[i-1])/(-2.0-m_cp[i-1]);   //OPTIMIZE HERE

      m_x2[m_msh.cintcells.size()-1]  = m_dp[m_msh.cintcells.size()-2] ;
      for(i=(m_msh.cintcells.size()-3);i>=0;i--) m_x2[i+1]  = m_dp[i]-m_cp[i]*m_x2[i+2];

      //Back Substitution into Original System

      m_phi[m_msh.cintcells.size()-1+m_msh.numghost] =  ((-(0.25)*(m_msh.pmesh[m_msh.cintcells.size()+1]-m_msh.pmesh[m_msh.cintcells.size()-1])*(m_msh.pmesh[m_msh.cintcells.size()+1]-m_msh.pmesh[m_msh.cintcells.size()-1])*ieps0*(m_rho[m_msh.cintcells.size()])) - m_x1[0]-m_x1[m_msh.cintcells.size()-1])/(-2.0+m_x2[0]+m_x2[m_msh.cintcells.size()-1]);

      for(i=0;i<m_msh.cintcells.size();i++) m_phi[i] = m_x1[i]+m_x2[i]*m_phi[m_msh.cintcells.size()-1+m_msh.numghost];

   }


}


//..Update applied electric field..//

void solver::updateAppliedEfield(std::vector<std::string> s_initE, std::vector<double> &s_EM_E, const mesh &s_msh)
{
   int i,j,k;
   int tot_cells;
   double temp1,temp2;
   std::vector<double> s_pos;

   for(i=0;i<s_msh.meshdims;i++) s_pos.push_back(0.0);


   if(s_msh.Eloc==0)
   {
      tot_cells = s_msh.cmesh.size()/s_msh.meshdims;
      for(i=0; i<tot_cells; i++)
      {
         for(j=0;j<s_msh.meshdims;j++) s_pos[j] = s_msh.cmesh[s_msh.meshdims*i+j];
         for(j=0;j<s_msh.vecdims; j++)   s_EM_E[i*s_msh.vecdims+j] = eqnparser(s_initE[j],s_pos);
      }
   }
   else if(s_msh.Eloc==1)
   {
      tot_cells = s_msh.pmesh.size()/s_msh.meshdims;
      for(i=0; i<tot_cells; i++)
      {
         for(j=0;j<s_msh.meshdims;j++) s_pos[j] = s_msh.pmesh[s_msh.meshdims*i+j];
         for(j=0;j<s_msh.vecdims; j++)   s_EM_E[i*s_msh.vecdims+j] = eqnparser(s_initE[j],s_pos);
      }
   }
}

//..Update applied magnetic field..//

void solver::updateAppliedBfield(std::vector<std::string> s_initB, std::vector<double> &s_EM_B, const mesh &s_msh)
{
   int i,j,k;
   int tot_cells;
   double temp1,temp2;
   std::vector<double> s_pos;

   for(i=0;i<s_msh.meshdims;i++) s_pos.push_back(0.0);


   if(s_msh.Bloc==0)
   {
      tot_cells = s_msh.cmesh.size()/s_msh.meshdims;
      for(i=0; i<tot_cells; i++)
      {
         for(j=0;j<s_msh.meshdims;j++) s_pos[j] = s_msh.cmesh[s_msh.meshdims*i+j];
         for(j=0;j<s_msh.vecdims; j++)   s_EM_B[i*s_msh.vecdims+j] = eqnparser(s_initB[j],s_pos);
      }
   }
   else if(s_msh.Bloc==1)
   {
      tot_cells = s_msh.pmesh.size()/s_msh.meshdims;
      for(i=0; i<tot_cells; i++)
      {
         for(j=0;j<s_msh.meshdims;j++) s_pos[j] = s_msh.pmesh[s_msh.meshdims*i+j];
         for(j=0;j<s_msh.vecdims; j++)   s_EM_B[i*s_msh.vecdims+j]= eqnparser(s_initB[j],s_pos);
      }
   }
}

//..Update applied potential..//

void solver::updateAppliedPotential(std::string s_initphi, std::vector<double> &s_EM_phi, const mesh &s_msh)
{
   int i,j,k;
   int tot_cells;
   double temp1,temp2;
   std::vector<double> s_pos;

   for(i=0;i<s_msh.meshdims;i++) s_pos.push_back(0.0);


   if(s_msh.philoc==0)
   {
      tot_cells = s_msh.cmesh.size()/s_msh.meshdims;
      for(i=0; i<tot_cells; i++)
      {
         for(j=0;j<s_msh.meshdims;j++) s_pos[j] = s_msh.cmesh[s_msh.meshdims*i+j];
         s_EM_phi[i] = eqnparser(s_initphi,s_pos);
      }
   }
   else if(s_msh.philoc==1)
   {
      tot_cells = s_msh.pmesh.size()/s_msh.meshdims;
      for(i=0; i<tot_cells; i++)
      {
         for(j=0;j<s_msh.meshdims;j++) s_pos[j] = s_msh.pmesh[s_msh.meshdims*i+j];
         s_EM_phi[i] = eqnparser(s_initphi,s_pos);
      }
   }
}

//..Calculate electric field from potential..//

void solver::phitoE(std::vector<double> &s_EM_phi, std::vector<double> &s_EM_E, const mesh &s_msh)
{
   int i,j,k;
   int num_cells;
   double temp1,temp2;

   mathFunctions mth;

   /* Use 1 dimension always
   if(s_msh.meshdims==1)
   {
   */
      //std::cout << "\n\tCalculating E-Field...";

      if(s_msh.philoc==0)
      {

         num_cells = s_msh.pmesh.size(); // I guess this automatically assumes that the E field is on the edges

         for(i=0;i<num_cells;i++)
         {
            //s_EM_E[s_msh.vecdims*i] = -mth.ddx1Dpwc(s_msh,s_EM_phi,i);
            s_EM_E[s_msh.vecdims*i] = -mth.ddx1Dpwc(s_msh,s_msh.deltax,s_EM_phi,i); //EFF
            //s_EM_E[s_msh.vecdims*i+1] = 0.0;
            //s_EM_E[s_msh.vecdims*i+2] = 0.0;
         }
      }
      else if(s_msh.philoc==1)
      {
         num_cells = s_msh.pmesh.size();

         for(i=1;i<(num_cells-1);i++)
         {
            s_EM_E[s_msh.vecdims*i] = -mth.ddx1Dpwp(s_msh,s_EM_phi,i);
            //s_EM_E[s_msh.vecdims*i+1] = 0.0;
            //s_EM_E[s_msh.vecdims*i+2] = 0.0;
         }
      /* Use 1 dimension always
      }
      */
   }

   //for(i=0;i<s_EM_E.size();i++) std::cout << s_EM_E[i];

   //std::cout << "x";

}

//..Equation parser..//

double solver::eqnparser(std::string expression_string, std::vector<double> s_point)
{
   int s_size = s_point.size();
   double s_x, s_y, s_z;
   exprtk::symbol_table<double> symbol_table;

   switch (s_size)
   {
      case 1:
         s_x = s_point[0];
         symbol_table.add_variable("x",s_x);
         break;
      case 2:
         s_x = s_point[0];
         s_y = s_point[1];
         symbol_table.add_variable("x",s_x);
         symbol_table.add_variable("y",s_y);
         symbol_table.add_variable("z",s_z);
         break;
      case 3:
         s_x = s_point[0];
         s_y = s_point[1];
         s_z = s_point[2];
         symbol_table.add_variable("x",s_x);
         symbol_table.add_variable("y",s_y);
         symbol_table.add_variable("z",s_z);
         break;
   }

   symbol_table.add_variable("t",totalTime);

   exprtk::expression<double> expression;
   expression.register_symbol_table(symbol_table);

   exprtk::parser<double> parser;
   parser.compile(expression_string,expression);

   return expression.value();
}


//..Flux source of particles..//

void solver::particlefluxsource1D(particles &s_part, const mesh &s_msh, std::string s_ddist, std::string s_thdist, std::vector<double> s_vel, double s_dens, double s_temp,double s_dt, int s_cell, int fdir,int extrapart)
{ 
   //.....Maxwellian set based on "Loading and Injection of Maxwellian Distributions in Particle Simulations" by Cartwright, Verboncoeur, and Birdsall....//

   int i,j,k,l;
   int npctot,npc,npo,npcrem,s_cindex,s_pindex;
   long double temp1,temp2,temp3,tempen;
   long double vupper,vlower,vtherm,i_dens,area,dnpc,vlowerf,vlfsq,vupperf,vufsq;
   double tol=0.0e-10;

   std::vector<int> neigh;
   std::vector<double> s_pos;
   std::vector<double> dx;
   std::vector<std::string> s_initvel;
   std::vector<long double> s_R,s_Rtot;

   int numprocs,procid; //MPI

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   npo = s_part.pos.size()/s_msh.meshdims;

   for(i=0;i<s_msh.meshdims;i++)  s_pos.push_back(s_msh.cmesh[s_msh.meshdims*s_cell+i]);
   for(i=0;i<s_msh.vecdims;i++)  s_initvel.push_back("0.0");
   for(i=0;i<s_msh.meshdims;i++)  dx.push_back(0.0);
   for(i=0;i<2*s_msh.meshdims;i++)  neigh.push_back(0);

   area = s_msh.carea[s_cell];
   for(j=0;j<2*s_msh.meshdims;j++) neigh[j] = s_msh.pofcneigh[2*s_cell+j];
   s_cindex = s_cell+fdir; 
   if(neigh[0]<0) s_pindex = neigh[1];
   else if(neigh[1]<0) s_pindex = neigh[0];

   mathFunctions s_mth;

   for(j=0;j<s_msh.meshdims;j++) dx[j] =  (fabs(s_msh.cmesh[s_msh.meshdims*s_cindex+j]-s_msh.cmesh[s_msh.meshdims*s_cell+j]));

   vtherm = sqrt(s_temp*kb/(s_part.mass));
   npc =  s_dt*(0.25)*sqrt(8.0/pi)*vtherm*s_dens*area/s_part.pwght;
   dnpc =  s_dt*(0.25)*sqrt(8.0/pi)*vtherm*s_dens*area/s_part.pwght;
   //npc =  s_dt*s_dens*s_vel[0]*area/s_part.pwght; //TMP: Bulk flux
   //dnpc =  s_dt*s_dens*s_vel[0]*area/s_part.pwght; //TMP: Bulk Flux

   //std::cout << std::endl << npc << std::endl;

   s_part.cseedcount[s_cell] += (dnpc-npc);

   if(s_part.cseedcount[s_cell]>1.0) 
   {
      npc += 1;
      s_part.cseedcount[s_cell] += -1.0;
   }

   s_part.gp = 0; 

   //npc += extrapart; 

   npctot = npc;
   npc = npctot/numprocs;
   npcrem = npctot - npc*numprocs;

   if(procid==0) npc = npc + npcrem;

   for(j=0;j<npc;j++) s_R.push_back((j+0.5)/npc);
   for(j=0;j<(s_msh.vecdims-1)*npc;j++) s_Rtot.push_back(0.0);

   auto seed = time(NULL)+procid;

   for(j=0;j<(s_msh.vecdims-1);j++)
   {
      seed = (j+1)*(time(NULL)+procid);
      std::shuffle(s_R.begin(),s_R.end(),std::default_random_engine(seed));
      for(k=0;k<npc;k++) s_Rtot[(s_msh.vecdims-1)*k+j] = s_R[k];
   }
   std::shuffle(s_R.begin(),s_R.end(),std::default_random_engine(seed));

   vupper = 5.0;
   vlower = -5.0;

   vupperf = 5.0;
   vufsq = vupperf*vupperf;
   vlowerf = 0.01;
   vlfsq = vlowerf*vlowerf;

   for(j=0;j<npc;j++)
   {
      tempen=0.0;

      if(s_thdist == "MAXWELLIAN")
      {
         //Flux direction goes here//

         if(s_vel[0]==0)
         {
            temp1 = rand();
            temp1 = temp1/RAND_MAX;
            temp1 = fdir*sqrt(vlfsq+vufsq - log(temp1*exp(vlfsq)+(1-temp1)*exp(vufsq)));
            s_part.vel.push_back(sqrt(2)*vtherm*temp1);
            tempen= tempen + temp1*temp1;
         }
         else
         {
            s_part.vel.push_back(s_vel[0]); //TMP: Use only bulk flux
         }

         for(k=0;k<(s_msh.vecdims-1);k++)
         {
            temp1=rand();
            temp1=temp1/RAND_MAX;
            temp1=s_mth.erfinv(temp1*erf(vupper)+(1-temp1)*erf(vlower));
            s_part.vel.push_back(sqrt(2.0)*vtherm*temp1+s_vel[k+1]);
            tempen= tempen + (sqrt(2.0)*vtherm*temp1+s_vel[k+1])*(sqrt(2.0)*vtherm*temp1+s_vel[k+1]);
         }
         s_part.en.push_back((0.5)*tempen*s_part.mass);   //Energy of particle in group, not of collection
      }
      else if(s_thdist == "QUIET")
      {
         //Flux direction goes here//

         if(s_vel[0]==0)
         {
            temp1 = s_R[j];
            temp1 = sqrt(vlfsq+vufsq - log(temp1*exp(vlfsq)+(1-temp1)*exp(vufsq)));
            s_part.vel.push_back(fdir*sqrt(2)*vtherm*temp1);
            tempen= tempen + temp1*temp1;
         }

         for(k=0;k<(s_msh.vecdims-1);k++)
         {
            temp1=s_Rtot[j*(s_msh.vecdims-1)+k];
            temp1=s_mth.erfinv(temp1*erf(vupper)+(1-temp1)*erf(vlower));
            s_part.vel.push_back(sqrt(2)*vtherm*temp1+s_vel[k+1]);
            tempen= tempen + (sqrt(2)*vtherm*temp1+s_vel[k+1])*(sqrt(2)*vtherm*temp1+s_vel[k+1]);
         }
         s_part.en.push_back((0.5)*tempen*s_part.mass); //Energy of particle in group, not collection
      }

      if(s_ddist == "RANDOM")
      {
         temp1 =  rand();
         temp1 =  s_dt*s_part.vel[s_msh.vecdims*(npo+j)]*temp1/RAND_MAX;
         temp1 =  fdir*temp1 + s_msh.pmesh[s_msh.meshdims*s_pindex];
         s_part.pos.push_back(temp1);
      }
      else if(s_ddist == "UNIFORM")
      {
         temp1 =  fdir*(j+0.5)*(s_dt*s_part.vel[s_msh.vecdims*(npo+j)]/npc) + s_msh.pmesh[s_msh.meshdims*s_pindex];
         s_part.pos.push_back(temp1);
      }

      s_part.cell.push_back(s_cindex);
      s_part.pt.push_back(s_msh.nearp1Dnew(s_part,npo+j));

   }

}

//..Particle source at cell..//

void solver::particlecellsource(particles &s_part, const mesh &s_msh, std::string s_ddist, std::string s_thdist, std::vector<double> s_vel, double s_dens, double s_temp, int s_cell)
{ 

   //.....Maxwellian set based on "Loading and Injection of Maxwellian Distributions in Particle Simulations" by Cartwright, Verboncoeur, and Birdsall....//

   int i,j,k,l;
   int npc,npctot,npcrem,npo,nearestp;
   long double temp1,temp2,temp3,tempen;
   long double vupper,vlower,vtherm,i_dens,area,dnpc;
   double tol=0.0e-10;
   int numpart,s_seed_proc;
   int numprocs,procid;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   std::vector<int> neigh;
   std::vector<double> s_pos;
   std::vector<double> dx;
   std::vector<std::string> s_initvel;
   std::vector<long double> s_R,s_Rtot;

   for(i=0;i<s_msh.meshdims;i++)  s_pos.push_back(s_msh.cmesh[s_msh.meshdims*s_cell+i]);
   for(i=0;i<s_msh.vecdims;i++)  s_initvel.push_back("0.0");
   for(i=0;i<s_msh.meshdims;i++)  dx.push_back(0.0);
   for(i=0;i<s_msh.meshdims;i++)  neigh.push_back(0);

   npo = s_part.pos.size()/s_msh.meshdims;

   area = s_msh.carea[s_cell];

   mathFunctions s_mth;

   for(j=0;j<2*s_msh.meshdims;j++) neigh[j] = s_msh.pofcneigh[2*s_cell+j];
   for(j=0;j<s_msh.meshdims;j++) dx[j] =  2*(fabs(s_msh.pmesh[s_msh.meshdims*neigh[2*j]+j]-s_msh.cmesh[s_msh.meshdims*s_cell+j]));
   for(j=0;j<s_msh.meshdims;j++) area =  area*dx[j];

   vtherm = sqrt(s_temp*kb/(s_part.mass));

#if SIMTYPE!=1
   npc = s_dens*area/s_part.pwght;
   dnpc = s_dens*area/s_part.pwght;
#else
   npc = s_dens/s_part.pwght;
   dnpc = s_dens/s_part.pwght;
#endif

   //std::cout << npc << "\t" << dnpc <<  "\t" << vtherm << std::endl;

   s_part.cseedcount[s_cell] += (dnpc-npc);

   if(s_part.cseedcount[s_cell]>1.0) 
   {
      npc += 1;
      s_part.cseedcount[s_cell] += -1.0;
   }

   npctot = npc;
   npc = npctot/numprocs;
   npcrem = npctot - npc*numprocs;

   if(procid==0) npc = npc + npcrem;

   for(j=0;j<npc;j++) s_R.push_back((j+0.5)/npc);

   //s_part.gp += npc; 
   s_part.gp = 0;  // LTR: Above for BC version 

   for(j=0;j<npc;j++) s_R.push_back((j+0.5)/npc);
   for(j=0;j<s_msh.vecdims*npc;j++) s_Rtot.push_back(0.0);

   auto seed = time(NULL)+procid;

   for(j=0;j<s_msh.vecdims;j++)
   {
      seed = (j+1)*(time(NULL)+procid);
      //std::random_shuffle(i_R.begin(),i_R.end());
      std::shuffle(s_R.begin(),s_R.end(),std::default_random_engine(seed));
      for(k=0;k<npc;k++) s_Rtot[s_msh.vecdims*k+j] = s_R[k];
   }

   vupper = 5.0;
   vlower = -5.0;

   temp1 = erf(1.5);
   temp1 = s_mth.erfinv(temp1);

   for(j=0;j<npc;j++)
   {
      if(s_ddist == "RANDOM")
      {
         for(k=0;k<s_msh.meshdims;k++)
         {
            temp1 =  rand();
            if(neigh[0]<0)  temp1 =  -dx[k]*temp1/(RAND_MAX)+s_msh.pmesh[s_msh.meshdims*neigh[1]+k];
            else  temp1 =  dx[k]*temp1/(RAND_MAX)+s_msh.pmesh[s_msh.meshdims*neigh[0]+k];
#if SIMTYPE!=1
            s_part.pos.push_back(temp1);
#else
            s_part.pos.push_back(0.0);
#endif

            temp1 = (temp1+s_msh.hdeltax)/s_msh.deltax;    //Forced rounding using truncation
            nearestp = temp1;
            s_part.pt.push_back(nearestp);
         }
      }
      else if(s_ddist == "UNIFORM")
      {
         for(k=0;k<s_msh.meshdims;k++)
         {
            if(neigh[0]<0)  temp1 =  -(j+0.5)*(dx[k]/npc)+s_msh.pmesh[s_msh.meshdims*neigh[1]+k];
            else   temp1 =  (j+0.5)*(dx[k]/npc)+s_msh.pmesh[s_msh.meshdims*neigh[0]+k];
#if SIMTYPE!=1
            s_part.pos.push_back(temp1);
#else
            s_part.pos.push_back(0.0);
#endif
            temp1 = (temp1+s_msh.hdeltax)/s_msh.deltax;    //Forced rounding using truncation
            nearestp = temp1;
            s_part.pt.push_back(nearestp);
         }
      }

      tempen=0.0;

      if(s_thdist == "MAXWELLIAN")
      {
         for(k=0;k<s_msh.vecdims;k++)
         {
            temp1=rand();
            temp1=temp1/RAND_MAX;
            temp1=s_mth.erfinv(temp1*erf(vupper)+(1-temp1)*erf(vlower));
            s_part.vel.push_back(sqrt(2.0)*vtherm*temp1+s_vel[k]);
            tempen= tempen + (sqrt(2.0)*vtherm*temp1+s_vel[k])*(sqrt(2.0)*vtherm*temp1+s_vel[k]);
         }
         s_part.en.push_back((0.5)*tempen*s_part.mass);   //Energy of particle in group, not collection
      }
      else if(s_thdist == "QUIET")
      {
         for(k=0;k<s_msh.vecdims;k++)
         {
            temp1=s_Rtot[j*s_msh.vecdims+k];
            temp1=s_mth.erfinv(temp1*erf(vupper)+(1-temp1)*erf(vlower));
            s_part.vel.push_back(sqrt(2)*vtherm*temp1+s_vel[k]);
            tempen= tempen + (sqrt(2)*vtherm*temp1+s_vel[k])*(sqrt(2)*vtherm*temp1+s_vel[k]);
         }
         s_part.en.push_back((0.5)*tempen*s_part.mass); //Energy of particle in group, not collection
      }

      s_part.cell.push_back(s_cell);
      s_part.pt.push_back(s_msh.nearp1Dnew(s_part,npo+j));
   }

}


//..Check particles and field for NAN's..//

void  solver::checkNAN(const particles s_part, const fields s_EM,const mesh s_msh)
{
   int i,j,k;
   bool nancheck;

   nancheck = false;

   for(i=0;i<(s_part.pos.size()/s_msh.meshdims);i++)
   {
      for(j=0;j<s_msh.vecdims;j++)  
      {
         if(isnan(s_part.vel[i*s_msh.vecdims+j])==true)
         {
            std::cout << "ERROR::NAN found in velocity in " << s_part.name << ":  "<< i << "\t" << j << std::endl;
            nancheck = true;
         }
      }
   }


   if(s_msh.Eloc==0)
   {
      for(i=0;i<(s_msh.cmesh.size()/s_msh.meshdims);i++)
      {
         for(j=0;j<s_msh.vecdims;j++)  
         {
            if(isnan(s_EM.E[i*s_msh.vecdims+j])==true)
            {
               std::cout << "ERROR::NAN found in E-Field" << i << "\t" << j << std::endl;
               nancheck = true;
            }
         }
      }
   }
   else if(s_msh.Eloc==1)
   {
      for(i=0;i<(s_msh.pmesh.size()/s_msh.meshdims);i++)
      {
         for(j=0;j<s_msh.vecdims;j++)  
         {
            if(isnan(s_EM.E[i*s_msh.vecdims+j])==true)
            {
               std::cout << "ERROR::NAN found in E-Field" << i << "\t" << j << std::endl;
               nancheck = true;
            }
         }
      }
   }

   if(s_msh.philoc==0)
   {
      for(i=0;i<(s_msh.cmesh.size()/s_msh.meshdims);i++)
      {
         if(isnan(s_EM.phi[i])==true)
         {
            std::cout << "ERROR::NAN found in Potential" << i << std::endl;
            nancheck = true;
         }
      }
   }
   else if(s_msh.philoc==1)
   {
      for(i=0;i<(s_msh.pmesh.size()/s_msh.meshdims);i++)
      {
         if(isnan(s_EM.phi[i])==true)
         {
            std::cout << "ERROR::NAN found in Potential" << i << std::endl;
            nancheck = true;
         }
      }
   }

   if(nancheck == true) exit(EXIT_FAILURE);

}


//...............Collide Particles..................//

//..According to Birdsall, 1991, "Particle-in-cell Charged-Particle Simulations, Plus Monte Carlo Collisions with Neutral Atoms, PIC-MCC"..//
//..Null collisions according to Vahedi, 1995, "A Monte Carlo collision model for the particle-in-cell method: applications to argon and oxygen discharges..//
//..Also a lot from Sydorenko, 2006, "PIC Simulations of Electron Dynamics in Low Pressure Discharges with Magnetic Fields"..//

void solver::collideParticles(std::vector<particles> &s_prt,const mesh &s_msh, contnm &s_neut, double s_time, int id)
{
   int j,k,iii,kk;
   int numpart,numcol;
   int s_ind,s_ind3,s_cell,s_pt,s_loc,s_tabind;
   int tempi;
   double B,weight;
   long double nu,numax;
   double maxen,maxvmag,maxneutN,maxneutT;
   long double P1,Pmax;
   long double R1,R2,R3,R4;
   double vmag_en,vmag;
   double en_inc,en_ex,en_ej,en_ion,en_eV,delta_en;
   long double tempd;
   long double chi,phi_ang,alpha;
   long double vupper,vlower,vth;
   double dt = deltaT;//s_prt[id].cycIT*deltaT;

   vupper = 5.0;
   vlower = -5.0;

   mathFunctions s_mth;

   std::vector<int> colarray;
   std::vector<double> vel;
   std::vector<double> velneut;
   std::vector<double> nu_array;
   std::vector<double> pos;
   std::vector<double> veldrift;

   for(j=0;j<s_msh.vecdims;j++) vel.push_back(0.0);
   for(j=0;j<s_msh.vecdims;j++) velneut.push_back(0.0);
   for(j=0;j<s_msh.vecdims;j++) veldrift.push_back(0.0);
   for(j=0;j<s_msh.meshdims;j++) pos.push_back(0.0);
   for(j=0;j<3;j++) nu_array.push_back(0.0);


   for(iii=0;iii<s_prt[id].pnct;iii++)
   {
      maxen = *std::max_element(s_prt[id].en.begin(),s_prt[id].en.end());
      maxvmag = sqrt(maxen*2.0/(s_prt[id].mass));
      maxneutN = *std::max_element(s_neut.N.begin(),s_neut.N.end());
      maxneutT= *std::max_element(s_neut.T.begin(),s_neut.T.end());
      numpart = s_prt[id].en.size();

      if(s_prt[id].cllsntype[iii]=="NEUTRAL" && s_prt[id].name=="electron")   //..Handles electron-neutral collisions
      {
         //if(s_prt[id].crsstype[iii]=="ARGON") B = 10*qe;
         //else if(s_prt[id].crsstype[iii]=="XENON") B = 8.7*qe;
         B = s_prt[id].cllsnB;

         if(s_prt[id].crsstype[iii]=="CONSTANT")
         {
            numax = 2.0*maxvmag*maxneutN*s_prt[id].crsssctn[iii];
            Pmax = 1.0 - exp(-numax*dt);

            s_prt[id].ncolcount += Pmax*numpart;
            numcol = s_prt[id].ncolcount;
            s_prt[id].ncolcount -= numcol;

#if SIMTYPE!=2
            randomselect(colarray,numpart,numcol);
#endif        

#if SIMTYPE==2
            for(j=0;j<numpart;j++)
#else
               for(j=0;j<numcol;j++)
#endif
               {
#if SIMTYPE==2
                  s_ind = j;
#else
                  s_ind = colarray[j];
#endif
                  s_ind3 = 3*s_ind;

                  s_cell = s_prt[id].cell[s_ind];
                  s_pt = s_prt[id].pt[s_ind];
                  if(s_msh.philoc==1) s_loc = s_pt;
                  else if(s_msh.philoc==0) s_loc = s_cell;

                  en_inc = s_prt[id].en[s_ind]; 
                  for(k=0;k<s_msh.vecdims;k++) vel[k] = s_prt[id].vel[s_ind3+k];

                  vmag_en = sqrt(2.0*en_inc/s_prt[id].mass);

                  tempd = 0.0;
                  for(k=0;k<s_msh.vecdims;k++)  tempd += vel[k]*vel[k];
                  vmag = sqrt(tempd);


                  nu = s_neut.N[s_loc]*vmag*s_prt[id].crsssctn[iii];
                  P1 = nu/numax;
                  R1 = (long double) rand()/RAND_MAX;

                  //std::cout << nu << "\t" << P1 << "\t" << R1 << std::endl;

#if SIMTYPE!=2
                  if(R1<P1)
                  {
#endif
                     en_inc = vmag*vmag*0.5*s_prt[id].mass;
                     en_eV = en_inc/qe;

                     /*....Elastic Type Collisions.....*/

                     R2 = (long double) rand()/RAND_MAX;
                     chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                     R3 = (long double) rand()/RAND_MAX;
                     phi_ang = R3*2.0*pi;

                     //std::cout  << "Ang:" << chi << "\t" << phi_ang << std::endl;

                     scatterParticle(vel,vmag,chi,phi_ang);

                     //delta_en = 2.0*s_prt[id].mass/s_prt[id]cllsnmass*(1.0-cos(chi))*en_inc;
                     delta_en = 0.0;
                     alpha = sqrt(1.0-delta_en/en_inc);

                     for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k];
                     s_prt[id].en[s_ind] = en_inc - delta_en;

                     //if(s_prt[id].en[s_ind] >1.6e-13) std::cout << "\n En1:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";
                     /*....Inelastic Collisions.....*/

                     /*delta_en = 3.3e-20;

                       if(en_inc > delta_en)
                       {
                       en_ex = en_inc-delta_en;
                       en_eV = en_ex/qe;

                       R2 = (long double) rand()/RAND_MAX;
                       chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                       R3 = (long double) rand()/RAND_MAX;
                       phi_ang = R3*2.0*pi;

                       scatterParticle(vel,vmag,chi,phi_ang);

                       alpha = sqrt(1-delta_en/en_inc);

                       for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k];
                       s_prt[id].en[s_ind] = en_ex;
                       }
                       */
                     /*....Ionization Collisions.....*/

                     /*en_ion = 3.3e-20;

                       if(en_inc > en_ion)
                       {
                       R2 = (long double) rand()/RAND_MAX;
                       en_ej = B*tan(R2*atan((en_inc-en_ion)/(2.0*B)));
                       en_ex = en_inc - en_ej - en_ion;

                       en_eV = en_ej;
                       R2 = (long double) rand()/RAND_MAX;
                       chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                       R3 = (long double) rand()/RAND_MAX;
                       phi_ang = R3*2.0*pi;

                       scatterParticle(vel,vmag,chi,phi_ang);

                       alpha = sqrt(en_ej/en_inc);

                       s_prt[id].pos.push_back(s_prt[id].pos[s_ind]);
                       for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel.push_back(alpha*vel[k]);
                       s_prt[id].en.push_back(en_ej);
                       s_prt[id].cell.push_back(s_prt[id].cell[s_ind]);
                       s_prt[id].pt.push_back(s_prt[id].pt[s_ind]);


                       for(k=0;k<s_msh.vecdims;k++) vel[k] = s_prt[id].vel[s_ind3+k];

                       en_eV = en_ex;
                       R2 = (long double) rand()/RAND_MAX;
                       chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                       R3 = (long double) rand()/RAND_MAX;
                       phi_ang = R3*2.0*pi;

                       scatterParticle(vel,vmag,chi,phi_ang);

                       alpha = sqrt(en_ex/en_inc);

                       for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k];
                       s_prt[id].en[s_ind] = en_ex ;
                       }*/
#if SIMTYPE!=2
                  }
#endif
               }
         }
         else
         {
            numax = 2.0*maxvmag*maxneutN*s_prt[id].crsssctn_max;
            Pmax = 1.0 - exp(-numax*dt);

            s_prt[id].ncolcount += Pmax*numpart;
            numcol = s_prt[id].ncolcount;
            s_prt[id].ncolcount -= numcol;

            randomselect(colarray,numpart,numcol);

            //std::cout << "\nCol:  " << numax << "\t" << Pmax << "\t" << numpart << "\t" << numcol;

            for(j=0;j<numcol;j++)
            {
               s_ind = colarray[j];
               s_ind3 = 3*s_ind;
               s_cell = s_prt[id].cell[s_ind];
               s_pt = s_prt[id].pt[s_ind];
               if(s_msh.philoc==1) s_loc = s_pt;
               else if(s_msh.philoc==0) s_loc = s_cell;

               en_inc = s_prt[id].en[s_ind]; 
               for(k=0;k<s_msh.vecdims;k++) vel[k] = s_prt[id].vel[s_ind3+k];
               vmag_en = sqrt(2.0*en_inc/s_prt[id].mass);

               tempd = 0.0;
               for(k=0;k<s_msh.vecdims;k++)  tempd += vel[k]*vel[k];
               vmag = sqrt(tempd);

               en_inc = vmag*vmag*0.5*s_prt[id].mass;
               en_eV = en_inc/qe;

               s_tabind = sqrt(log(en_eV+1.0)*1000.0)-1;

               if(s_tabind>81) 
               { 
                  s_tabind = 81;
                  weight = 1.0;
               }
               else if(s_tabind<0)
               {
                  s_tabind=0;
                  weight=1.0;
               }
               else 
               {
                  weight = (en_inc-s_prt[id].en_crsstable[s_tabind])/(s_prt[id].en_crsstable[s_tabind+1]-s_prt[id].en_crsstable[s_tabind]);
               }

               nu_array[0] = s_neut.N[s_loc]*vmag*(s_prt[id].el_crsstable[s_tabind]*(1.0-weight) + s_prt[id].el_crsstable[s_tabind+1]*(weight)); 
               nu_array[1] = nu_array[0]+s_neut.N[s_loc]*vmag*(s_prt[id].inel_crsstable[s_tabind]*(1.0-weight) + s_prt[id].inel_crsstable[s_tabind+1]*(weight)); 
               nu_array[2] = nu_array[1]+s_neut.N[s_loc]*vmag*(s_prt[id].ion_crsstable[s_tabind]*(1.0-weight) + s_prt[id].ion_crsstable[s_tabind+1]*(weight)); 

               R1 = (long double) rand()/RAND_MAX;

               //std::cout << "\nNu:  " << nu_array[0] << "\t" << nu_array[1] << "\t" << nu_array[2] << "\t" << numcol;

               if(R1<nu_array[0]/numax)
               {
                  s_prt[id].elcount += 1;

                  R2 = (long double) rand()/RAND_MAX;
                  chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                  R3 = (long double) rand()/RAND_MAX;
                  phi_ang = R3*2.0*pi;

                  scatterParticle(vel,vmag,chi,phi_ang);

                  delta_en = 2.0*s_prt[id].mass/s_prt[id].cllsnmass*(1.0-cos(chi))*en_inc;
                  //delta_en = 0.0;
                  alpha = sqrt(1.0-delta_en/en_inc);

                  for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k];
                  s_prt[id].en[s_ind] = en_inc - delta_en; 

                  for(k=0;k<s_msh.vecdims;k++)  if(std::isnan(s_prt[id].vel[s_ind3+k])) std::cout << "ElV:  " << s_ind << "\t" << k; 
                  if(std::isnan(s_prt[id].en[s_ind])) std::cout << "ElEn:  " << s_ind << "\t" << k;
                  //if(s_prt[id].en[s_ind] >1.6e-13) std::cout << "\n En2:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";
               }
               else if(R1<nu_array[1]/numax && en_inc > s_prt[id].cllsnenergy[1])
               {
                  s_prt[id].exccount += 1;

                  delta_en = s_prt[id].cllsnenergy[1];
                  en_ex = en_inc-delta_en;
                  en_eV = en_ex/qe;

                  R2 = (long double) rand()/RAND_MAX;
                  chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                  R3 = (long double) rand()/RAND_MAX;
                  phi_ang = R3*2.0*pi;

                  scatterParticle(vel,vmag,chi,phi_ang);

                  alpha = sqrt(1.0-delta_en/en_inc);

                  for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k];
                  s_prt[id].en[s_ind] = en_ex;

                  for(k=0;k<s_msh.vecdims;k++) if(std::isnan(s_prt[id].vel[s_ind3+k])) std::cout << "\nExcV:  " << s_ind << "\t" << k << "\t" << chi << "\t" << phi_ang << "\t" << alpha ;
                  if(std::isnan(s_prt[id].en[s_ind])) std::cout << "\nExcEn:  " << s_ind << "\t" << k << "\t" << chi << "\t" << phi_ang << "\t" << alpha << "\t" << en_ex;
                  //if(s_prt[id].en[s_ind] >1.6e-13) std::cout << "\n En3:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";
               }
               else if(R1<nu_array[2]/numax && en_inc > s_prt[id].cllsnenergy[2])
               {
                  s_prt[id].ioncount += 1;
                  s_neut.ioncount[s_loc] += s_prt[id].pwght;  //NEW

                  for(k=0;k<s_msh.meshdims;k++) pos[k] = s_prt[id].pos[s_msh.meshdims*s_ind+k];

                  en_ion = s_prt[id].cllsnenergy[2];
                  R2 = (long double) rand()/RAND_MAX;
                  en_ej = B*tan(R2*atan((en_inc-en_ion)/(2.0*B)));
                  en_ex = en_inc - en_ej - en_ion;

                  //std::cout << en_ex;

                  en_eV = en_ej/qe;
                  R2 = (long double) rand()/RAND_MAX;
                  chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                  R3 = (long double) rand()/RAND_MAX;
                  phi_ang = R3*2.0*pi;

                  scatterParticle(vel,vmag,chi,phi_ang);

                  alpha = sqrt(en_ej/en_inc);

                  s_prt[id].pos.push_back(s_prt[id].pos[s_ind]);
                  for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel.push_back(alpha*vel[k]);
                  s_prt[id].en.push_back(en_ej);
                  s_prt[id].cell.push_back(s_prt[id].cell[s_ind]);
                  s_prt[id].pt.push_back(s_prt[id].pt[s_ind]);

                  for(k=0;k<s_msh.vecdims;k++) if(std::isnan(s_prt[id].vel[s_ind3+3+k])) std::cout << "IonVEj:  " << s_ind << "\t" << k << "\t" << chi << "\t" << phi_ang << "\t" << alpha ;
                  if(std::isnan(s_prt[id].en[s_ind+1])) std::cout << "IonEnEj:  " << s_ind << "\t" << k << "\t" << chi << "\t" << phi_ang << "\t" << alpha << "\t" << en_ex;
                  //if(en_ej >1.6e-13) std::cout << "\n En4:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";


                  for(k=0;k<s_msh.vecdims;k++) vel[k] = s_prt[id].vel[s_ind3+k];

                  en_eV = en_ex/qe;
                  R2 = (long double) rand()/RAND_MAX;
                  chi = acos((2.0+en_eV-2.0*pow(1.0+en_eV,R2))/en_eV);

                  R3 = (long double) rand()/RAND_MAX;
                  phi_ang = R3*2.0*pi;

                  scatterParticle(vel,vmag,chi,phi_ang);

                  alpha = sqrt(en_ex/en_inc);

                  for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k];
                  s_prt[id].en[s_ind] = en_ex ;

                  for(k=0;k<s_msh.vecdims;k++) if(std::isnan(s_prt[id].vel[s_ind3+k])) std::cout << "IonVEx:  " << s_ind << "\t" << k << "\t" << chi << "\t" << phi_ang << "\t" << alpha ;
                  if(std::isnan(s_prt[id].en[s_ind])) std::cout << "IonEnEx:  " << s_ind << "\t" << chi << "\t" << phi_ang << "\t" << alpha << "\t" << en_ex;
                  //if(en_ej >1.6e-13) std::cout << "\n En5:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";

                  for(k=0;k<s_prt[id].nsp;k++)      //Seed Ion
                  {
                     if(s_prt[id].crsstype[iii]==s_prt[k].name)  seedSingleParticle(s_prt[k],s_msh,veldrift,pos,s_neut.T[s_loc],s_cell,s_pt); //CHK  s_cell
                  }
               }
            }
         }
      }
      else if(s_prt[id].cllsntype[iii]=="NEUTRAL" && s_prt[id].name!="electron")   //..Handles ion-neutral collisions
      {
         if(s_prt[id].crsstype[iii]=="CONSTANT")
         {
            numax = 2.0*maxvmag*maxneutN*s_prt[id].crsssctn[iii];
            Pmax = 1.0 - exp(-numax*dt);

            s_prt[id].ncolcount += Pmax*numpart;
            numcol = s_prt[id].ncolcount;
            s_prt[id].ncolcount -= numcol;

#if SIMTYPE!=2
            randomselect(colarray,numpart,numcol);
#endif        

#if SIMTYPE==2
            for(j=0;j<numpart;j++)
#else
               for(j=0;j<numcol;j++)
#endif
               {
#if SIMTYPE==2
                  s_ind = j;
#else
                  s_ind = colarray[j];
#endif
                  s_ind3 = s_ind*3;

                  //....Particle Info...//

                  s_cell = s_prt[id].cell[s_ind];
                  s_pt = s_prt[id].pt[s_ind];
                  if(s_msh.philoc==1) s_loc = s_pt;
                  else if(s_msh.philoc==0) s_loc = s_cell;

                  en_inc = s_prt[id].en[s_ind]; 
                  vmag_en = sqrt(2.0*en_inc/s_prt[id].mass);
                  for(k=0;k<s_msh.vecdims;k++) vel[k] = s_prt[id].vel[s_ind3+k];

                  //..Shift frames...//

                  vth = sqrt(kb*s_neut.T[s_loc]/s_prt[id].mass);
                  for(k=0;k<s_msh.vecdims;k++)
                  {
                     tempd= (long double) rand()/RAND_MAX;
                     tempd=s_mth.erfinv(tempd*erf(vupper)+(1.0-tempd)*erf(vlower));
                     tempd = sqrt(2.0)*vth*tempd;
                     velneut[k] = tempd;
                     vel[k] -= velneut[k];
                  }

                  tempd=0.0;
                  for(k=0;k<s_msh.vecdims;k++)  tempd += vel[k]*vel[k];
                  vmag = sqrt(tempd);

                  //..Calculate collision freq...//

                  nu = s_neut.N[s_loc]*vmag*s_prt[id].crsssctn[iii];

                  P1 = nu/numax;
                  R1 = (long double) rand()/RAND_MAX;

#if SIMTYPE!=2
                  if(R1<P1)
                  {
#endif
                     en_inc = vmag*vmag*0.5*s_prt[id].mass;
                     en_eV = en_inc/qe;

                     /*.........Ion Elastic.........*/

                     R2 = (long double) rand()/RAND_MAX;
                     chi = acos(sqrt(1.0-R2));

                     R3 = (long double) rand()/RAND_MAX;
                     phi_ang = 2.0*pi*R3;

                     en_ex = en_inc*cos(chi)*cos(chi);
                     delta_en = en_inc-en_ex;

                     alpha = sqrt(en_ex/en_inc);
                     scatterParticle(vel,vmag,chi,phi_ang);

                     for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k] + velneut[k];

                     tempd = 0.0;
                     for(k=0;k<s_msh.vecdims;k++) tempd += s_prt[id].vel[s_ind3+k]*s_prt[id].vel[s_ind3+k];

                     s_prt[id].en[s_ind] = tempd*s_prt[id].mass*0.5;

                     /*.............CEX............*/

                     /*for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = velneut[k];

                       tempd=0.0;
                       for(k=0;k<s_msh.vecdims;k++) tempd += velneut[k]*velneut[k];

                       s_prt[id].en[s_ind] = tempd*0.5*s_prt[id].mass;
                       */
#if SIMTYPE!=2
                  }
#endif 
               }

         }
         else
         {
            vth = sqrt(kb*maxneutT/s_prt[id].mass);

            numax = 2.0*(maxvmag+vth)*maxneutN*s_prt[id].crsssctn_max;
            Pmax = 1.0 - exp(-numax*dt);

            s_prt[id].ncolcount += Pmax*numpart;
            numcol = s_prt[id].ncolcount;
            s_prt[id].ncolcount -= numcol;

            randomselect(colarray,numpart,numcol);

            //std::cout << "\nCol:  " << numax << "\t" << Pmax << "\t" << numpart << "\t" << numcol;

            for(j=0;j<numcol;j++)
            {
               //....Particle Info...//

               s_ind = colarray[j];
               s_ind3 = 3*s_ind;
               s_cell = s_prt[id].cell[s_ind];
               s_pt = s_prt[id].pt[s_ind];
               if(s_msh.philoc==1) s_loc = s_pt;
               else if(s_msh.philoc==0) s_loc = s_cell;

               en_inc = s_prt[id].en[s_ind]; 
               vmag_en = sqrt(2.0*en_inc/s_prt[id].mass);
               for(k=0;k<s_msh.vecdims;k++) vel[k] = s_prt[id].vel[s_ind3+k];

               //..Shift frames...//

               vth = sqrt(kb*s_neut.T[s_loc]/s_prt[id].mass);
               for(k=0;k<s_msh.vecdims;k++)
               {
                  tempd= (long double) rand()/RAND_MAX;
                  if(tempd==0) tempd=1.0/RAND_MAX;
                  if(tempd==1) tempd=(RAND_MAX-1)/RAND_MAX;
                  tempd=s_mth.erfinv(tempd*erf(vupper)+(1.0-tempd)*erf(vlower));
                  tempd = sqrt(2.0)*vth*tempd;
                  velneut[k] = tempd;
                  vel[k] -= velneut[k];
               }

               tempd=0.0;
               for(k=0;k<s_msh.vecdims;k++)  tempd += vel[k]*vel[k];
               vmag = sqrt(tempd);

               //..Calculate collision freq...//

               en_inc = vmag*vmag*0.5*s_prt[id].mass;
               en_eV = en_inc/qe;

               s_tabind = sqrt(log(en_eV+1.0)*1000)-1;

               if(s_tabind>81) 
               { 
                  s_tabind = 81;
                  weight = 1.0;
               }
               else if(s_tabind<0)
               {
                  s_tabind=0;
                  weight=1.0;
               }
               else 
               {
                  weight = (en_inc-s_prt[id].en_crsstable[s_tabind])/(s_prt[id].en_crsstable[s_tabind+1]-s_prt[id].en_crsstable[s_tabind]);
               }

               nu_array[0] = s_neut.N[s_loc]*vmag*(s_prt[id].el_crsstable[s_tabind]*(1.0-weight) + s_prt[id].el_crsstable[s_tabind+1]*(weight)); 
               nu_array[1] = nu_array[0]+s_neut.N[s_loc]*vmag*(s_prt[id].inel_crsstable[s_tabind]*(1.0-weight) + s_prt[id].inel_crsstable[s_tabind+1]*(weight)); 
               nu_array[2] = nu_array[1]+s_neut.N[s_loc]*vmag*(s_prt[id].ion_crsstable[s_tabind]*(1.0-weight) + s_prt[id].ion_crsstable[s_tabind+1]*(weight)); 

               R1 = (long double) rand()/RAND_MAX;

               //std::cout << "\nNu:  " << nu_array[0] << "\t" << nu_array[1] << "\t" << nu_array[2] << "\t" << numcol;

               if(R1<nu_array[0]/numax)
               {
                  /*.........Ion Elastic.........*/
                  s_prt[id].elcount += 1;

                  R2 = (long double) rand()/RAND_MAX;
                  chi = acos(sqrt(1.0-R2));

                  R3 = (long double) rand()/RAND_MAX;
                  phi_ang = 2.0*pi*R3;

                  en_ex = en_inc*cos(chi)*cos(chi);
                  delta_en = en_inc-en_ex;
                  alpha = sqrt(en_ex/en_inc);
                  scatterParticle(vel,vmag,chi,phi_ang);

                  //std::cout << "\nElastic:  " << s_ind << "\t" << en_inc << "\t" << en_ex << "\t" << alpha;

                  for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = alpha*vel[k] + velneut[k];

                  tempd = 0.0;
                  for(k=0;k<s_msh.vecdims;k++) tempd += s_prt[id].vel[s_ind3+k]*s_prt[id].vel[s_ind3+k];

                  s_prt[id].en[s_ind] = tempd*s_prt[id].mass*0.5;
                  //std::cout << "\nVel:  " << vel[0] << "\t" << vel[1] << "\t" << vel[2] << "\t" << tempd;
                  //if(s_prt[id].en[s_ind] >1.6e-13) std::cout << "\n En6:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";
               }
               else if(R1<nu_array[1]/numax)
               {
                  s_prt[id].exccount += 1;
                  /*.............CEX............*/

                  for(k=0;k<s_msh.vecdims;k++) s_prt[id].vel[s_ind3+k] = velneut[k];

                  tempd=0.0;
                  for(k=0;k<s_msh.vecdims;k++) tempd += velneut[k]*velneut[k];
                  s_prt[id].en[s_ind] = tempd*0.5*s_prt[id].mass;
                  //if(s_prt[id].en[s_ind] >1.6e-13) std::cout << "\n En7:  " << s_prt[id].en[s_ind] << "\t" << s_prt[id].vel[s_ind3] << "\t" << s_prt[id].vel[s_ind3+1] << "\t" << s_prt[id].vel[s_ind3+2] << "\t";
               }
            }
         }
      }
      //std::cout << "\nParticle: "   << id << "   Elastic:  " << s_prt[id].elcount*s_prt[0].pwght/s_time << "\t" << "Inelastic: " << s_prt[id].exccount*s_prt[0].pwght/s_time << "\tIonization: "<< s_prt[id].ioncount*s_prt[0].pwght/s_time  << std::endl; //CV
   }

   //std::cout << "\nION Elastic: "   << s_prt[1].elcount*s_prt[1].pwght/s_time << "\t" << "Inelastic: "  << s_prt[1].exccount*s_prt[1].pwght/s_time << std::endl; //CV

}


/*.................OLD COLLISION ALGORITHM..................*/


/*void solver::collideParticles(std::vector<particles> &s_part,const mesh &s_msh, contnm &s_neutrals, double s_time, int id)
  {
  int i,j,k;
  int s_cell, s_pt, s_loc;
  int s_pind,s_pindT,s_pind3;
  int s_tabind;
  int rand_proc;
  int numpart,numcol,numpart_total,numrem,numavg;
  double numcold;
  double maxen,maxvmag,maxNneut;
  double Pmax,nu_max;
  double coef,inv_hmass,inv_sqrt,vmag,vmagsq;
  long double P,nu;
  long double P1,P2,P3,P4;
  long double R1,R2,R3,R4,R5;
  double phi_ang,chi_ang;
  double delta_en,alpha,weight;
  double en_inc,en_ej,en_eV,en_new,B;
  double temp,vtherm,vupper,vlower;
  double tol = 1.0e-10;
  double invsqrt2 = 1.0/sqrt(2.0);
  double colpar;
  double dt = deltaT;//s_prt[id].cycIT*deltaT;

  std::vector<double> vel;
  std::vector<double> velneut;
  std::vector<int> col_array;
  std::vector<double> nu_array;
  std::vector<double> s_pos;
  std::vector<double> s_drift;

  int numprocs,procid; //MPI
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

  for(i=0;i<3;i++) nu_array.push_back(0.0);
  for(i=0;i<s_msh.meshdims;i++) s_pos.push_back(0.0);
  for(i=0;i<s_msh.vecdims;i++) s_drift.push_back(0.0);

  if(s_msh.vecdims != 3)
  {
  std::cout << "\n.....Must have 3 velocity dimensions for collisions......\n";
  exit(EXIT_FAILURE);
  }

  for(i=0;i<3;i++) vel.push_back(0.0); //..Must have 3 dimensions!
  for(i=0;i<3;i++) velneut.push_back(0.0); //..Must have 3 dimensions!

  mathFunctions s_mth;  //..Declare math function class

  numpart_total=0;
  numpart = s_part[id].pos.size()/s_msh.meshdims;
  MPI_Allreduce(&numpart,&numpart_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

  maxen = *std::max_element(s_part[id].en.begin(), s_part[id].en.end());
  maxvmag = sqrt(maxen*2.0/s_part[id].mass);

//std::cout << "\nmaxnen:  " << maxen << "\t" << maxvmag << std::endl;

//std::cout << "\n\tColliding Particles...";

for(i=0;i<s_part[id].pnct;i++)
{
//..............Electron Neutral Collisions.....................//
if(s_part[id].cllsntype[i]=="NEUTRAL" && s_part[id].name=="electron")   //..Electron Neutral Collisions
{
maxNneut = *std::max_element(s_neutrals.N.begin(), s_neutrals.N.end());

//if(s_part[id].crsstype[i]=="ARGON") B = 10.0*qe;      //..Set integration constants
//else if(s_part[id].crsstype[i]=="XENON") B = 8.7*qe;  
B = s_part[id].cllsnB;  

//................CONSTANT Cross Section.......................///

if(s_part[id].crsstype[i]=="CONSTANT")  //..Constant Cross section
{
   coef = s_part[id].crsssctn[i]*dt;
   inv_hmass = 2.0/s_part[id].mass;
   //nu_max = 5.0*maxvmag*maxNneut*s_part[id].crsssctn[i];   //..Calculate Maximum Collision Frequency
   nu_max = 2.0*maxvmag*maxNneut*s_part[id].crsssctn[i];   //..Calculate Maximum Collision Frequency

   Pmax = 1.0 - exp(-nu_max*dt);

   s_part[id].ncolcount += Pmax*numpart;
   numcol = s_part[id].ncolcount;
   s_part[id].ncolcount -= numcol;

   /*numavg = numcol/numprocs;
     numrem = numcol - numavg*numprocs;
     temp = rand();
     rand_proc = temp/RAND_MAX*numprocs;
     if(procid==rand_proc) numcol = numavg+numrem;
     else numcol = numavg;*/

   /*         #if SIMTYPE!=2
              randomselect(col_array,numpart,numcol); //..CV
#endif


#if SIMTYPE==2
for(j=0;j<numpart;j++)//..CV
#else
for(j=0;j<col_array.size();j++)
#endif
{
#if SIMTYPE==2
s_pind = j;  //..CV
#else
s_pind = col_array[j]; //..CV
#endif

s_cell = s_part[id].cell[s_pind];
s_pt = s_part[id].pt[s_pind];
if(s_msh.philoc==1) s_loc = s_pt;
else if(s_msh.philoc==0) s_loc = s_cell;

en_inc = s_part[id].en[s_pind];

vmag =   sqrt(inv_hmass*en_inc);
nu = s_neutrals.N[s_loc]*vmag*s_part[id].crsssctn[i];

R1 = rand();
R1 = (long double) R1/RAND_MAX;

P1 = nu/nu_max;

#if SIMTYPE!=2
if(R1<P1) //CCV
{         //CCV
#endif        
   //std::cout << "\n........COLLISION.........." << std::endl;

   s_pind3 = 3*s_pind;
   en_eV = en_inc/qe;
   alpha =1.0;
   delta_en=0.0;
   //inv_sqrt = 1.0/sqrt(s_part.vel[s_pind3+1]*s_part.vel[s_pind3+1] + s_part.vel[s_pind3+2]*s_part.vel[s_pind3+2]);
   //inv_sqrt = 1.0/sqrt(s_part[id].vel[s_pind3]*s_part[id].vel[s_pind3] + s_part[id].vel[s_pind3+1]*s_part[id].vel[s_pind3+1]);

   R3 = (long double) rand(); 
   phi_ang = 2.0*pi*R3/RAND_MAX;
   R4 = rand(); 
   R4 = (long double) R4/RAND_MAX; 
   chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

   for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

   scatterParticle(vel,vmag,chi_ang,phi_ang);

   //vel[0] = s_part.vel[s_pind3]*cos(chi_ang)    +  sin(chi_ang)*cos(phi_ang)/inv_sqrt; 
   //vel[1] = s_part.vel[s_pind3+1]*cos(chi_ang)  +  sin(chi_ang)*sin(phi_ang)*vmag*s_part.vel[s_pind3+2]*inv_sqrt - sin(chi_ang)*cos(phi_ang)*s_part.vel[s_pind3]*s_part.vel[s_pind3+1]*inv_sqrt; 
   //vel[2] = s_part.vel[s_pind3+2]*cos(chi_ang)  -  sin(chi_ang)*sin(phi_ang)*vmag*s_part.vel[s_pind3+1]*inv_sqrt - sin(chi_ang)*cos(phi_ang)*s_part.vel[s_pind3]*s_part.vel[s_pind3+2]*inv_sqrt;

   //vel[0] = s_part[id].vel[s_pind3]*cos(chi_ang)    +  s_part[id].vel[s_pind3+1]*vmag*sin(chi_ang)*sin(phi_ang)*inv_sqrt   +   s_part[id].vel[s_pind3]*s_part[id].vel[s_pind3+2]*sin(chi_ang)*cos(phi_ang)*inv_sqrt; 
   //vel[1] = s_part[id].vel[s_pind3+1]*cos(chi_ang)  -  s_part[id].vel[s_pind3]*vmag*sin(chi_ang)*sin(phi_ang)*inv_sqrt     +   s_part[id].vel[s_pind3+1]*s_part[id].vel[s_pind3+2]*sin(chi_ang)*cos(phi_ang)*inv_sqrt; 
   //vel[2] = s_part[id].vel[s_pind3+2]*cos(chi_ang)  - sin(chi_ang)*cos(phi_ang)/inv_sqrt; 

   s_part[id].vel[s_pind3] = vel[0];
   s_part[id].vel[s_pind3+1] = vel[1];
   s_part[id].vel[s_pind3+2] = vel[2];


   /*------------Check INELASTIC-----------------*//*
                                                      s_pind3 = 3*s_pind;
                                                      en_eV = en_inc/qe;

                                                      delta_en = 3.3e-20;
                                                      alpha = sqrt(1-delta_en/en_inc);
                                                      en_eV = en_eV-delta_en/qe;

                                                      if(en_eV>=0)   //..Check energy, may be negative due to interpolation
                                                      {
                                                      R3 = rand(); 
                                                      phi_ang = 2.0*pi*R3/RAND_MAX;
                                                      R4 = rand(); 
                                                      R4 = R4/RAND_MAX; 
                                                      chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

                                                      for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

                                                      scatterParticle(vel,vmag,chi_ang,phi_ang);

                                                      s_part[id].vel[s_pind3] = alpha*vel[0];
                                                      s_part[id].vel[s_pind3+1] = alpha*vel[1];
                                                      s_part[id].vel[s_pind3+2] = alpha*vel[2];
                                                      s_part[id].en[s_pind] = en_inc-delta_en;
                                                      }

*//*------------Check INELASTIC-----------------*/


   /*------------Check IONIZATION-----------------*//*

                                                       s_pind3 = 3*s_pind;
                                                       en_eV = en_inc/qe;

                                                       R2 = rand();
                                                       R2 = R2/RAND_MAX;

                                                       B = 0.01*3.197e-20;      //..Set integration constants
   //en_ej = B*tan(R2*atan((en_inc-3.3e-20)/(2.0*B)));
   en_ej = 3.197e-20;
   delta_en = 3.197e-20+en_ej;

   //..Calculation for incident electron

   alpha = sqrt(1.0-delta_en/en_inc);
   en_eV = en_eV-delta_en/qe;

   if(en_eV>=0)   //..Check energy, may be negative due to interpolation
   {
   R3 = rand(); 
   phi_ang = 2.0*pi*R3/RAND_MAX;
   R4 = rand(); 
   R4 = R4/RAND_MAX; 
   chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

   for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];


   scatterParticle(vel,vmag,chi_ang,phi_ang);

   s_part[id].vel[s_pind3]   = alpha*vel[0];
   s_part[id].vel[s_pind3+1] = alpha*vel[1];
   s_part[id].vel[s_pind3+2] = alpha*vel[2];
   s_part[id].en[s_pind] = en_inc-delta_en;

   //std::cout << "\nen1: " << en_inc-delta_en << std::endl;

   //..Calculation for ejected electron

   alpha = sqrt(1-(en_inc-en_ej)/en_inc);
   en_eV = en_ej/qe;

   R3 = rand(); 
   phi_ang = 2.0*pi*R3/RAND_MAX;
   R4 = rand(); 
   R4 = R4/RAND_MAX; 
   chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

   //for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

   scatterParticle(vel,vmag,chi_ang,phi_ang);

   //..Make new electron

   s_part[id].pos.push_back(s_part[id].pos[s_pind]); // 1D
   s_part[id].vel.push_back(alpha*vel[0]);
   s_part[id].vel.push_back(alpha*vel[1]);
   s_part[id].vel.push_back(alpha*vel[2]);
   s_part[id].en.push_back(en_ej);
   s_part[id].cell.push_back(s_cell);
   s_part[id].pt.push_back(s_pt);

   //std::cout << "\ne2:  " << alpha << "\t" << en_eV << "\t" << en_ej << "\t" << s_loc<< std::endl ; 


   for(k=0;k<s_msh.meshdims;k++) s_pos[k] = s_part[id].pos[s_pind+k];

   for(k=0;k<s_part[id].nsp;k++)   //.....Check other species and seed ion if matching species exists
   {
   if("electron"==s_part[k].name) seedSingleParticle(s_part[k],s_msh,s_drift,s_pos,s_neutrals.T[s_loc],s_cell,s_pt); 
} 

}
*//*------------Check IONIZATION-----------------*/


/*           #if SIMTYPE!=2
             } //..CCV
#endif
}
}
else
{
//................Tabular Cross Section.......................///

inv_hmass = 2.0/s_part[id].mass;
nu_max = 2.0*maxvmag*s_part[id].crsssctn_max*maxNneut; //..Calculate maximum collision frequency
//nu_max = 5.0*s_part[id].invmfp_max*maxNneut; //..Calculate maximum collision frequency

Pmax = 1.0 - exp(-nu_max*dt);

s_part[id].ncolcount += Pmax*numpart;
numcol = s_part[id].ncolcount;
s_part[id].ncolcount -= numcol;

/*numcold = Pmax*numpart_total;
numcold = numcold-numcol;
temp = rand();
if(temp/RAND_MAX>numcold) numcol += 1;
*/

/*numavg = numcol/numprocs;
  numrem = numcol - numavg*numprocs;
  temp = rand();
  rand_proc = temp/RAND_MAX*numprocs;
  if(procid==rand_proc) numcol = numavg+numrem;
  else numcol = numavg;*/

//std::cout <<"\nelec:\t" << procid << "\t" << s_part[id].ncolcount << "\t" <<  numcol << "\n";

/*         randomselect(col_array,numpart,numcol); 

           for(j=0;j<col_array.size();j++)
           {
           s_pind = col_array[j]; //..CV
           s_cell = s_part[id].cell[s_pind];
           s_pt = s_part[id].pt[s_pind];
           if(s_msh.philoc==1) s_loc = s_pt;
           else if(s_msh.philoc==0) s_loc = s_cell;

           en_inc = s_part[id].en[s_pind];
           vmag =   sqrt(inv_hmass*en_inc);
           en_eV = en_inc/qe;

           s_tabind = sqrt(log(en_eV+1.0)*1000)-1;

           if(s_tabind>81) 
           { 
//std::cout << "\ns_tabind: \t" <<  s_tabind << std::endl;
//std::cout << s_pind << "\t" << s_part[id].en[s_pind] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+1] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+2] << s_part[id].pos[s_pind] <<std::endl; 
s_tabind = 81;
weight = 1.0;
}
else if(s_tabind<0)
{
//std::cout << "\ns_tabind: \t" <<  s_tabind << std::endl; 
//std::cout << s_pind << "\t" << s_part[id].en[s_pind] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+1] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+2] << s_part[id].pos[s_pind] <<std::endl; 
s_tabind=0;
weight=1;
}
else 
{
weight = (en_inc-s_part[id].en_crsstable[s_tabind])/(s_part[id].en_crsstable[s_tabind+1]-s_part[id].en_crsstable[s_tabind]);
}

//std::cout << "en_inc,tabind,weight:  " << en_inc << "\t" << s_tabind << "\t" << weight;

nu_array[0] = s_neutrals.N[s_loc]*vmag*(s_part[id].el_crsstable[s_tabind]*(1.0-weight) + s_part[id].el_crsstable[s_tabind+1]*(weight)); 
nu_array[1] = nu_array[0]+s_neutrals.N[s_loc]*vmag*(s_part[id].inel_crsstable[s_tabind]*(1.0-weight) + s_part[id].inel_crsstable[s_tabind+1]*(weight)); 
nu_array[2] = nu_array[1]+s_neutrals.N[s_loc]*vmag*(s_part[id].ion_crsstable[s_tabind]*(1.0-weight) + s_part[id].ion_crsstable[s_tabind+1]*(weight)); 

R1 = rand();
R1 = R1/RAND_MAX;

P1 = nu_array[2]/nu_max;


//std::cout << "\nPind,en_inc,vmag,cell:   " << s_pind << "\t" << en_inc << "\t" << vmag << "\t" <<  s_loc << std::endl;
//std::cout << "\nNu:   " << nu_max << "\t" << nu_array[0] << "\t" << nu_array[1] << "\t" <<  nu_array[2] << std::endl;
//std::cout << "\nR1,P1,en_eV,ind: " << R1 << "\t" << P1 << "\t" << en_eV << "\t" << s_tabind << std::endl;

if(R1<P1)   //...Undergo collision
{                 
#if SIMTYPE==2
std::cout << "\n........COLLISION.........." << std::endl;
#endif

nu_array[0] = nu_array[0]/nu_max; 
nu_array[1] = nu_array[1]/nu_max; 

if(R1<nu_array[0]) //...Elastic Collision
{
s_part[id].elcount += 1; //CV
#if SIMTYPE==2
std::cout << "\n........ELASTIC.........." << std::endl;
std::cout << "\n Number of Elastic Collisions:"   << s_part[i].elcount << "\t" << procid << std::endl; //CV
std::cout << "\n Collision Frequency:"   << nu_array[0]*nu_max << std::endl; //CV
#endif

s_pind3 = 3*s_pind;

R3 = rand(); 
phi_ang = 2.0*pi*R3/RAND_MAX;
R4 = rand(); 
R4 = R4/RAND_MAX; 
chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);
delta_en = 2.0*s_part[id].mass/s_part[id].cllsnmass*(1.0-cos(chi_ang))*en_inc;
alpha = sqrt(1.0-delta_en/en_inc);

for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

scatterParticle(vel,vmag,chi_ang,phi_ang);

s_part[id].vel[s_pind3]   = alpha*vel[0];
s_part[id].vel[s_pind3+1] = alpha*vel[1];
s_part[id].vel[s_pind3+2] = alpha*vel[2];

s_part[id].en[s_pind] = en_inc-delta_en;

//std::cout << "\nEl:  " << en_inc << "\t" << s_part[id].en[s_pind] << "\t" << 0.5*s_part[id].mass*(s_part[id].vel[s_pind3]*s_part[id].vel[s_pind3]+s_part[id].vel[s_pind3+1]*s_part[id].vel[s_pind3+1]+s_part[id].vel[s_pind3+2]*s_part[id].vel[s_pind3+2]); 
}
//else if(R1<nu_array[1]) //...Inelastic Collision
else if(R1<nu_array[1] && en_inc>s_part[id].cllsnenergy[1])      //Check if energy enough for ionizing
{
   s_part[id].exccount += 1;  //CV
#if SIMTYPE==2
   std::cout << "\n........IN-ELASTIC.........." << std::endl;
   std::cout << "\n Number of In-elastic Collisions:"   << s_part[i].exccount << "\t" << procid << std::endl;  //CV
#endif

   s_pind3 = 3*s_pind;

   delta_en = s_part[id].cllsnenergy[1];
   alpha = sqrt(1-delta_en/en_inc);
   //std::cout << "\nen_eV1:  " << en_eV << std::endl;
   en_eV = en_eV-delta_en/qe;
   //std::cout << "\nen_eV2:  " << en_eV << std::endl;



   //if(en_eV>=0)   //..Check energy, may be negative due to interpolation
   //{
   R3 = rand(); 
   phi_ang = 2.0*pi*R3/RAND_MAX;
   R4 = rand(); 
   R4 = R4/RAND_MAX; 
   chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

   for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

   scatterParticle(vel,vmag,chi_ang,phi_ang);

   s_part[id].vel[s_pind3] = alpha*vel[0];
   s_part[id].vel[s_pind3+1] = alpha*vel[1];
   s_part[id].vel[s_pind3+2] = alpha*vel[2];
   s_part[id].en[s_pind] = en_inc-delta_en;
   //std::cout << "\nExc:  " << en_inc << "\t" << s_part[id].en[s_pind] << "\t" << 0.5*s_part[id].mass*(s_part[id].vel[s_pind3]*s_part[id].vel[s_pind3]+s_part[id].vel[s_pind3+1]*s_part[id].vel[s_pind3+1]+s_part[id].vel[s_pind3+2]*s_part[id].vel[s_pind3+2]); 
   //}
}
//else //...Ionization Collision
else if(en_inc>s_part[id].cllsnenergy[2])      //Check if energy enough for ionizing
{
   s_part[id].ioncount += 1;  //CV
#if SIMTYPE==2
   std::cout << "\n........IONIZATION........." << std::endl;
   std::cout << "\n Number of Ionization collisions:"   << s_part[id].ioncount << "\t" << procid << std::endl;  //CV
#endif
   //std::cout << "\n........IONIZATION........." << std::endl;
   s_neutrals.ioncount[s_loc] += s_part[id].pwght;  //NEW

   s_pind3 = 3*s_pind;

   R2 = rand();
   R2 = R2/RAND_MAX;

   //if(en_inc>s_part[id].cllsnenergy[2])      //Check if energy enough for ionizing
   //{
   en_ej = B*tan(R2*atan((en_inc-s_part[id].cllsnenergy[2])/(2.0*B)));
   delta_en = s_part[id].cllsnenergy[2]+en_ej;

   //..Calculation for incident electron

   alpha = sqrt(1.0-delta_en/en_inc);
   //std::cout << "\nIen_eV1:  " << en_eV << std::endl;
   en_eV = en_eV-delta_en/qe;
   //std::cout << "\nIen_eV2:  " << en_eV << std::endl;

   //if(en_eV>=0)   //..Check energy, may be negative due to interpolation
   //{
   R3 = rand(); 
   phi_ang = 2.0*pi*R3/RAND_MAX;
   R4 = rand(); 
   R4 = R4/RAND_MAX; 
   chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

   for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

   scatterParticle(vel,vmag,chi_ang,phi_ang);

   s_part[id].vel[s_pind3]   = alpha*vel[0];
   s_part[id].vel[s_pind3+1] = alpha*vel[1];
   s_part[id].vel[s_pind3+2] = alpha*vel[2];
   s_part[id].en[s_pind] = en_inc-delta_en;

   //std::cout << "\nIon1:  " << s_part[id].en[s_pind] << "\t" << 0.5*s_part[id].mass*(s_part[id].vel[s_pind3]*s_part[id].vel[s_pind3]+s_part[id].vel[s_pind3+1]*s_part[id].vel[s_pind3+1]+s_part[id].vel[s_pind3+2]*s_part[id].vel[s_pind3+2]); 

   //std::cout << "\nen1: " << en_inc-delta_en << std::endl;

   //..Calculation for ejected electron

   alpha = sqrt(1-(en_inc-en_ej)/en_inc);
   en_eV = en_ej/qe;

   R3 = rand(); 
   phi_ang = 2.0*pi*R3/RAND_MAX;
   R4 = rand(); 
   R4 = R4/RAND_MAX; 
   chi_ang = acos((2+en_eV-2*pow((1.0+en_eV),R4))/en_eV);

   //for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];   //FIX

   scatterParticle(vel,vmag,chi_ang,phi_ang);

   //..Make new electron

   s_part[id].pos.push_back(s_part[id].pos[s_pind]); // 1D
   s_part[id].vel.push_back(alpha*vel[0]);
   s_part[id].vel.push_back(alpha*vel[1]);
   s_part[id].vel.push_back(alpha*vel[2]);
   s_part[id].en.push_back(en_ej);
   s_part[id].cell.push_back(s_cell);
   s_part[id].pt.push_back(s_pt);

   s_pindT = s_part[id].en.size()-1;
   //std::cout << "\ne2:  " << alpha << "\t" << en_eV << "\t" << en_ej << "\t" << s_loc<< std::endl ; 
   //std::cout << "\nIon2:  " << en_inc<< "\t" <<  en_ej << "\t" << 0.5*s_part[id].mass*alpha*alpha*(vel[0]*vel[0]+vel[1]*vel[1]+vel[2]*vel[2]); 
   //std::cout << "\nIon2*:  " << s_part[id].en[s_pindT] << "\t" << 0.5*s_part[id].mass*(s_part[id].vel[s_pindT*3]*s_part[id].vel[s_pindT*3]+s_part[id].vel[s_pindT*3+1]*s_part[id].vel[s_pindT*3+1]+s_part[id].vel[s_pindT*3+2]*s_part[id].vel[s_pindT*3+2]); 


   for(k=0;k<s_msh.meshdims;k++) s_pos[k] = s_part[id].pos[s_pind+k];

   for(k=0;k<s_part[id].nsp;k++)   //.....Check other species and seed ion if matching species exists
   {
      if(s_part[id].crsstype[i]==s_part[k].name) seedSingleParticle(s_part[k],s_msh,s_drift,s_pos,s_neutrals.T[s_loc],s_cell,s_pt); //CHK s_cell 
   } 
   //}
   //}
}
} //..CV
}
}
}
else if(s_part[id].cllsntype[i]=="NEUTRAL")   //..Ion Neutral Collisions
{
   maxNneut = *std::max_element(s_neutrals.N.begin(), s_neutrals.N.end());
   vupper = 6.0;
   vlower = -6.0;

   if(s_part[id].crsstype[i]=="CONSTANT")  //..Constant Cross section
   {
      coef = s_part[id].crsssctn[i]*dt;
      inv_hmass = 2.0/s_part[id].mass;
      nu_max = 2.0*maxvmag*maxNneut*s_part[id].crsssctn[i];   //..Calculate Maximum Collision Frequency

      Pmax = 1.0 - exp(-nu_max*dt);
      //numcol = Pmax*numpart;

      s_part[id].ncolcount += Pmax*numpart;
      numcol = s_part[id].ncolcount;
      s_part[id].ncolcount -= numcol;

      /*numavg = numcol/numprocs;
        numrem = numcol - numavg*numprocs;
        temp = rand();
        rand_proc = temp/RAND_MAX*numprocs;
        if(procid==rand_proc) numcol = numavg+numrem;
        else numcol = numavg;*/

      /*         #if SIMTYPE!=2
                 randomselect(col_array,numpart,numcol); //..CV
#endif

#if SIMTYPE==2
for(j=0;j<numpart;j++)//..CV
#else
for(j=0;j<col_array.size();j++)
#endif
{
#if SIMTYPE==2
s_pind = j;  //..CV
#else
s_pind = col_array[j]; //..CV
#endif
s_pind3 = 3*s_pind;
s_cell = s_part[id].cell[s_pind];
s_pt = s_part[id].pt[s_pind];
if(s_msh.philoc==1) s_loc = s_pt;
else if(s_msh.philoc==0) s_loc = s_cell;

vtherm = sqrt(s_neutrals.T[s_loc]*kb/s_part[id].mass);

      //std::cout << std::endl << vtherm << "\t" << s_neutrals.T[s_loc] << std::endl;

      for(k=0;k<s_msh.vecdims;k++)   //Find Random Neutral Velocity
      {
      temp=rand();
      temp=temp/RAND_MAX;
      temp=s_mth.erfinv(temp*erf(vupper)+(1.0-temp)*erf(vlower));
      velneut[k] = sqrt(2.0)*vtherm*temp;
      }

      en_inc = 0.0;

      for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

      for(k=0;k<s_msh.vecdims;k++) //Changing frame  
      {  
      vel[k] = vel[k] - velneut[k];
      en_inc += vel[k]*vel[k];
      }
      //if(std::isnan(vel[0])==true || std::isinf(vel[0])==true || vel[0]==0)  std::cout << ".CF.0.." << "\t" << vel[0]  << std::endl ; 
      //if(std::isnan(vel[1])==true || std::isinf(vel[1])==true || vel[1]==0)  std::cout << ".CF.1.." << "\t" << vel[1]  << std::endl ; 
      //if(std::isnan(vel[2])==true || std::isinf(vel[2])==true || vel[2]==0)  std::cout << ".CF.2.." << "\t" << vel[2]  << std::endl ; 

      vmag =   sqrt(en_inc);
      en_inc = en_inc*s_part[id].mass*0.5;           

      //if(std::isnan(en_inc)==true || std::isinf(en_inc)==true)  std::cout << "..en_inc.." << std::endl << en_inc  << std::endl ; 

      nu = s_neutrals.N[s_loc]*vmag*s_part[id].crsssctn[i];

      R1 = rand();
      R1 = R1/RAND_MAX;

      P1 = nu/nu_max;

#if SIMTYPE!=2
if(R1<P1) //CV
{         //CV
#endif        
      //std::cout << "\n........COLLISION.........." << std::endl;

      en_eV = en_inc/qe;

      R3 = rand(); 
      phi_ang = 2.0*pi*R3/RAND_MAX;
      R4 = rand(); 
      R4 = R4/RAND_MAX; 
      chi_ang = acos(sqrt(1.0-R4));
      //delta_en = en_inc*(sin(chi_ang)*sin(chi_ang));
      alpha = sqrt(1.0-sin(chi_ang)*sin(chi_ang));


      //if(std::isnan(alpha)==true || std::isinf(alpha)==true)  std::cout << "alpha" << std::endl << alpha  << std::endl ; 
      //if(std::isnan(chi_ang)==true || std::isinf(chi_ang)==true)  std::cout << "chi" << std::endl << chi_ang  << std::endl ; 
      //if(std::isnan(vmag)==true || std::isinf(vmag)==true)  std::cout << "vmag" << std::endl << vmag  << std::endl ; 

      scatterParticle(vel,vmag,chi_ang,phi_ang);
      vel[0] = alpha*vel[0];
      vel[1] = alpha*vel[1];
      vel[2] = alpha*vel[2];

      //if(std::isnan(vel[0])==true || std::isinf(vel[0])==true)  std::cout << "NF..0" << std::endl << vel[0]  << std::endl ; 
      //if(std::isnan(vel[1])==true || std::isinf(vel[1])==true)  std::cout << "NF..1" << std::endl << vel[1]  << std::endl ; 
      //if(std::isnan(vel[2])==true || std::isinf(vel[2])==true)  std::cout << "NF..2" << std::endl << vel[2]  << std::endl ; 

      for(k=0;k<s_msh.vecdims;k++) vel[k] += velneut[k];  //Changing back frames

      s_part[id].vel[s_pind3] = vel[0];
      s_part[id].vel[s_pind3+1] = vel[1];
      s_part[id].vel[s_pind3+2] = vel[2];

      //if(std::isnan(vel[0])==true)  std::cout << "CB..0" << std::endl << vel[0]  << std::endl ; 
      //if(std::isnan(vel[1])==true)  std::cout << "CB..1" << std::endl << vel[1]  << std::endl ; 
      //if(std::isnan(vel[2])==true)  std::cout << "CB..2" << std::endl << vel[2]  << std::endl ; 

      //if(std::isnan(velneut[0])==true)  std::cout << "VN..0" << std::endl << velneut[0]  << std::endl ; 
      //if(std::isnan(velneut[1])==true)  std::cout << "VN..1" << std::endl << velneut[1]  << std::endl ; 
      //if(std::isnan(velneut[2])==true)  std::cout << "VN..2" << std::endl << velneut[2]  << std::endl ; 

      en_new=0.0;
      for(k=0;k<s_msh.vecdims;k++) en_new += vel[k]*vel[k];  //Calculate new energy
      s_part[id].en[s_pind] = 0.5*en_new*s_part[id].mass;
      //if(std::isnan(en_new)==true)  std::cout << "..en_new.." << std::endl << en_new  << std::endl ; 
      //if(std::isnan(s_part[id].en[s_pind])==true)  std::cout << "..en.." << std::endl << en_new  << std::endl ; 

      /*.........Check CEX...........*//*
                                          for(k=0;k<s_msh.vecdims;k++) vel[k] = velneut[k];  //Changing back frames

                                          s_part[id].vel[s_pind3] = vel[0];
                                          s_part[id].vel[s_pind3+1] = vel[1];
                                          s_part[id].vel[s_pind3+2] = vel[2];

                                          en_new=0.0;
                                          for(k=0;k<s_msh.vecdims;k++) en_new += vel[k]*vel[k];  //Calculate new energy
                                          s_part[id].en[s_pind] = 0.5*en_new*s_part[id].mass;
                                          *//*.........Check CEX...........*/

      /*           #if SIMTYPE!=2
                   } //..CCV
#endif
}
}
else
{
      //................Tabular Cross Section.......................///

      inv_hmass = 2.0/s_part[id].mass;
      nu_max = 2.0*maxvmag*s_part[id].crsssctn_max*maxNneut; //..Calculate maximum collision frequency
      //nu_max = 5.0*s_part[id].invmfp_max*maxNneut; //..Calculate maximum collision frequency

      //std::cout << "\nPmax" << Pmax << "\t" << s_part[id].crsssctn_max << "\t"  << maxNneut  << "t" << nu_max << std::endl;

      Pmax = 1.0 - exp(-nu_max*dt);
      //numcol = Pmax*numpart;
      s_part[id].ncolcount += Pmax*numpart;
      numcol = s_part[id].ncolcount;
      s_part[id].ncolcount -= numcol;

      /*numavg = numcol/numprocs;
      numrem = numcol - numavg*numprocs;
      temp = rand();
      rand_proc = temp/RAND_MAX*numprocs;
      if(procid==rand_proc) numcol = numavg+numrem;
      else numcol = numavg;*/

      //std::cout << "\nneut:\t" << procid << "\t" << s_part[id].ncolcount << "\t" <<  numcol << "\n";

      /*         randomselect(col_array,numpart,numcol); 

                 for(j=0;j<col_array.size();j++)
                 {
                 s_pind = col_array[j]; //..CV
                 s_cell = s_part[id].cell[s_pind];
                 s_pt = s_part[id].pt[s_pind];
                 if(s_msh.philoc==1) s_loc = s_pt;
                 else if(s_msh.philoc==0) s_loc = s_cell;
                 s_pind3 = 3*s_pind;

                 vtherm = sqrt(s_neutrals.T[s_loc]*kb/s_part[id].mass);

                 for(k=0;k<s_msh.vecdims;k++)   //Find Random Neutral Velocity
                 {
                 temp=rand();
                 temp=temp/RAND_MAX;
                 temp=s_mth.erfinv(temp*erf(vupper)+(1-temp)*erf(vlower));
                 velneut[k] = sqrt(2.0)*vtherm*temp;
                 }

                 en_inc = 0.0;

                 for(k=0;k<s_msh.vecdims;k++) vel[k]=s_part[id].vel[s_pind3+k];

                 for(k=0;k<s_msh.vecdims;k++) //Changing frame  
                 {  
                 vel[k] = vel[k] - velneut[k];
                 en_inc += vel[k]*vel[k];
                 }

                 vmag =   sqrt(en_inc);
                 en_inc = en_inc*s_part[id].mass*0.5;           
                 en_eV = en_inc/qe;

                 s_tabind = sqrt(log(en_eV+1.0)*1000)-1;

                 if(s_tabind>81)  //Hard coded size of cross-section table 
                 { 
      //std::cout << "\ns_tabind: \t" <<  s_tabind << std::endl; 
      //std::cout << s_pind << "\t" << s_part[id].en[s_pind] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+1] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+2] << s_part[id].pos[s_pind] <<std::endl; 
      s_tabind = 81;
      weight =1;
      }
      else if(s_tabind<0)
      {
      //std::cout << "\ns_tabind: \t" <<  s_tabind << std::endl; 
      //std::cout << s_pind << "\t" << s_part[id].en[s_pind] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+1] << "\t" << s_part[id].vel[s_pind*s_msh.vecdims+2] << s_part[id].pos[s_pind] << std::endl; 
      s_tabind = 0;
      weight = 1.0;
      }
      else
      {             
      weight = (en_inc-s_part[id].en_crsstable[s_tabind])/(s_part[id].en_crsstable[s_tabind+1]-s_part[id].en_crsstable[s_tabind]);
      }

      nu_array[0] = s_neutrals.N[s_loc]*vmag*(s_part[id].el_crsstable[s_tabind]*(1.0-weight) + s_part[id].el_crsstable[s_tabind+1]*(weight)); 
      nu_array[1] = nu_array[0]+s_neutrals.N[s_loc]*vmag*(s_part[id].inel_crsstable[s_tabind]*(1.0-weight) + s_part[id].inel_crsstable[s_tabind+1]*(weight)); 
      nu_array[2] = nu_array[1]+s_neutrals.N[s_loc]*vmag*(s_part[id].ion_crsstable[s_tabind]*(1.0-weight) + s_part[id].ion_crsstable[s_tabind+1]*(weight)); 

      R1 = rand();
      R1 = R1/RAND_MAX;

      P1 = nu_array[2]/nu_max;

      //std::cout << "\nPind,en_inc,vmag,cell:   " << s_pind << "\t" << en_inc << "\t" << vmag << "\t" <<  s_loc << std::endl;
      //std::cout << "\nNu:   " << nu_max << "\t" << nu_array[0] << "\t" << nu_array[1] << "\t" <<  nu_array[2] << std::endl;
      //std::cout << "\nR1,P1,en_eV,ind: " << R1 << "\t" << P1 << "\t" << en_eV << "\t" << s_tabind << std::endl;

      if(R1<P1)   //...Undergo collision
      {                 
#if SIMTYPE==2
      std::cout << "\n........COLLISION.........." << std::endl;
#endif

      nu_array[0] = nu_array[0]/nu_max; 
      nu_array[1] = nu_array[1]/nu_max; 

      if(R1<nu_array[0]) //...Elastic Collision
      {
         s_part[id].elcount += 1; //CV
#if SIMTYPE==2
         std::cout << "\n........ELASTIC.........." << std::endl;
         std::cout << "\n Number of Ion Elastic collisions:"   << s_part[id].elcount <<"\t" << procid << std::endl;  //CV
#endif
         //std::cout << "\n........ELASTIC.........." << std::endl;

         R3 = rand(); 
         phi_ang = 2.0*pi*R3/RAND_MAX;
         R4 = rand(); 
         R4 = R4/RAND_MAX; 
         chi_ang = acos(sqrt(1.0-R4));
         delta_en = en_inc*(sin(chi_ang)*sin(chi_ang));
         alpha = 1.0-sin(chi_ang)*sin(chi_ang);      //FIX!!

         scatterParticle(vel,vmag,chi_ang,phi_ang);
         vel[0] = alpha*vel[0];
         vel[1] = alpha*vel[1];
         vel[2] = alpha*vel[2];

         for(k=0;k<s_msh.vecdims;k++) vel[k] += velneut[k];  //Changing back frames

         s_part[id].vel[s_pind3] = vel[0];
         s_part[id].vel[s_pind3+1] = vel[1];
         s_part[id].vel[s_pind3+2] = vel[2];

         en_new=0.0;
         for(k=0;k<s_msh.vecdims;k++) en_new += vel[k]*vel[k];  //Calculate new energy
         //s_part[id].en[s_pind] = en_inc-delta_en;//0.5*en_new*s_part[id].mass;
         s_part[id].en[s_pind] =0.5*en_new*s_part[id].mass;
      }   
      else if(R1<nu_array[1]) //...CEX Collision
      {
         s_part[id].exccount += 1; //CV
#if SIMTYPE==2
         std::cout << "\n........CEX.........." << std::endl;
         std::cout << "\n Number of CEX collisions:"   << s_part[id].exccount << "\t" << procid << std::endl;  //CV
#endif
         //std::cout << "\n........CEX.........." << std::endl;

         for(k=0;k<s_msh.vecdims;k++) vel[k] = velneut[k];  //Changing back frames

         s_part[id].vel[s_pind3] = vel[0];
         s_part[id].vel[s_pind3+1] = vel[1];
         s_part[id].vel[s_pind3+2] = vel[2];

         en_new=0.0;
         for(k=0;k<s_msh.vecdims;k++) en_new += vel[k]*vel[k];  //Calculate new energy
         s_part[id].en[s_pind] = 0.5*en_new*s_part[id].mass;
      }
   } //..CCV
   }
   }
}
else
{
   std::cout << "\n\n.............COLLISION TYPE NOT SUPPORTED................\n\n";
}

}
//std::cout << "x";

//   std::cout << "\nELECTRON Elastic: "   << s_part[0].elcount*s_part[0].pwght/s_time << "\t" << "Inelastic: " << s_part[0].exccount*s_part[0].pwght/s_time << "\tIonization: "<< s_part[0].ioncount*s_part[0].pwght/s_time  << std::endl; //CV
// std::cout << "\nION Elastic: "   << s_part[1].elcount*s_part[1].pwght/s_time << "\t" << "Inelastic: "  << s_part[1].exccount*s_part[1].pwght/s_time << std::endl; //CV
}*/


///...................Coulomb Collisions...................///

void solver::coulombCollisions(std::vector<particles> &s_prt, const mesh &s_msh, std::vector<contnm> &s_cont, MPIvars &s_mpiv)
{
   int i,j,k,sp,csp;
   int vd = s_msh.vecdims;
   int pt,ptvd;

   long double tempd,en_temp;
   long double coef1,coef2;
   double mr;
   double lnC = 10.0;
   long double numom,nuE;
   long double f1,f2,f3;
   double vupper = 5.0;
   double vlower = -5.0;
   long double vthA;
   long double tol=69e-12; 

   std::vector<long double> dV;
   std::vector<long double> vth;
   std::vector<long double> vthsq;
   std::vector<long double> velA;

   if(s_cont[0].nsp>2)  std::cout << "\n\nCOULOMB COLLISIONS ONLY CONFIGURED FOR FIRST TWO SPECIES\n\n"; 

   //std::cout << "\nStarting Coulomb....";

   mathFunctions s_mth;

   mr = s_prt[0].mass*s_prt[1].mass/(s_prt[0].mass+s_prt[1].mass);
   for(i=0;i<vd;i++) velA.push_back(0.0);

   coulForce = 1;

   weighContinuumPIC(s_prt,s_cont,s_msh,s_mpiv);

   for(sp=0;sp<s_prt[0].nsp;sp++)  //Intra-species collisions
   {
      coef1 = 4.0/3.0*sqrt(pi)*s_prt[sp].charge*s_prt[sp].charge*s_prt[sp].charge*s_prt[sp].charge*lnC/(s_prt[sp].mass*s_prt[sp].mass)/(4.0*4.0*pi*pi*eps0*eps0);


      for(i=0;i<s_prt[sp].en.size();i++)
      {
         pt = s_prt[sp].pt[i];
         ptvd = 3*pt;

         numom = coef1*s_cont[sp].N[pt]*pow(kb*(s_cont[sp].T[pt]+tol)/s_prt[sp].mass,-1.5);

         vthA = sqrt(2.0*numom*deltaT*kb*s_cont[sp].T[pt]/s_prt[sp].mass);

         //std::cout << "\nvthA:  " << vthA << std::endl;

         for(j=0;j<vd;j++)
         {
            tempd=  (long double) rand();
            tempd= (tempd)/(RAND_MAX);
            tempd= s_mth.erfinv(tempd*erf(vupper)+(1.0-tempd)*erf(vlower));
            tempd = sqrt(2.0)*vthA*tempd;
            velA[j] = tempd;
            //std::cout << "\t " << velA[j];
         }

         for(j=0;j<s_msh.vecdims;j++) s_prt[sp].vel[i*vd+j] += -deltaT*numom*(s_prt[sp].vel[i*vd+j]-s_cont[sp].U[ptvd+j]) + velA[j];

         en_temp = 0.0;

         for(j=0;j<s_msh.vecdims;j++)
         {
            en_temp += s_prt[sp].vel[i*vd+j]*s_prt[sp].vel[i*vd+j];  // Avg before and after
            if(std::isnan(s_prt[sp].vel[i*vd+j]) || std::isinf(s_prt[sp].vel[i*vd+j])) std::cout << "\nCoul1:  " << i << "\t" << j << "\t" << velA[j] << "\t" << vthA << "\t" << numom << "\t" << coef1 << "\t" << s_cont[sp].T[pt]; 
         }
         s_prt[sp].en[i] = fabs(en_temp*(0.5)*s_prt[sp].mass);    
      }
   }

   weighContinuumPIC(s_prt,s_cont,s_msh,s_mpiv); //Suggested by Jones et. al

   if(s_cont[0].nsp>1)  //...Interspecies Collisions
   {
      //std::cout << "\nInterspecies....";

      for(i=0;i<s_cont[0].N.size();i++)
      {
         tempd = 0.0;
         for(j=0;j<vd;j++) tempd += pow((s_cont[0].U[i*vd+j] - s_cont[1].U[i*vd+j]),2);
         dV.push_back(sqrt(tempd));

         vthsq.push_back(2.0*(kb*s_cont[0].T[i]/s_prt[0].mass  +  kb*s_cont[1].T[i]/s_prt[1].mass));
         vth.push_back(sqrt(vthsq[i]));

         //std::cout << "\ni,dV,vth:  " << i << "\t" << dV[i] <<"\t" << vth[i];

      }

      coef1 = 8.0*sqrt(pi)*s_prt[0].charge*s_prt[0].charge*s_prt[1].charge*s_prt[1].charge*lnC/(mr*mr*4*4*pi*pi*eps0*eps0);
      coef2 = 16.0*sqrt(pi)*s_prt[0].charge*s_prt[0].charge*s_prt[1].charge*s_prt[1].charge*lnC/(s_prt[0].mass*s_prt[1].mass*4*4*pi*pi*eps0*eps0);

      for(sp=0;sp<s_cont[0].nsp;sp++)
      {

         if(sp==0) csp=1;
         else if(sp==1) csp=0;

         //std::cout << "\nSpecies:  " << sp << " with " << csp;

         for(i=0;i<s_prt[sp].en.size();i++)
         {
            en_temp = 0.0;
            pt = s_prt[sp].pt[i];
            ptvd = pt*vd;

            numom = coef1*s_cont[csp].N[pt]/(dV[pt]*dV[pt]*dV[pt])*(0.5*sqrt(pi)*erf(dV[pt]/vth[pt]) - dV[pt]/vth[pt]*exp(-dV[pt]*dV[pt]/vthsq[pt]));
            nuE = coef2*s_cont[csp].N[pt]/(vth[pt]*vthsq[pt])*exp(-dV[pt]*dV[pt]/vthsq[pt]);
            if(std::isnan(numom) || std::isnan(nuE)) std::cout << "\nCoul2nu's:  " << numom << "\t" << nuE << "\t" << dV[pt] << "\t" << vth[pt] << "\t" << coef1 << "\t" << coef2; 
            //nuE = coef2*s_cont[csp].N[pt]/(vth[pt]*vthsq[pt]);
            //std::cout << "\nnumom:  " << numom;
            //std::cout << "\nnuE:  " << nuE;

            for(j=0;j<vd;j++)
            {
               velA[j] = s_prt[sp].vel[i*vd+j];

               f1 = numom*mr*(s_cont[csp].U[ptvd+j]-s_cont[sp].U[ptvd+j]);  
               f2 = -numom*mr*mr/s_prt[sp].mass * pow((s_cont[csp].U[ptvd+j]-s_cont[sp].U[ptvd+j]),2)*(s_cont[sp].U[ptvd+j]-s_prt[sp].vel[i*vd+j])/(s_cont[sp].Tvec[ptvd+j]*kb/s_prt[sp].mass - s_cont[sp].U[ptvd+j]*s_cont[sp].U[ptvd+j]+tol);
               f3 = nuE*(kb*s_cont[sp].T[pt]-kb*s_cont[csp].T[pt])*(s_cont[sp].U[ptvd+j]-s_prt[sp].vel[i*vd+j])/(3.0*(s_cont[sp].T[pt]*kb/s_prt[sp].mass)+tol);
               //f3 = nuE*(kb*s_cont[sp].T[pt]-kb*s_cont[csp].T[pt])*(0.0-s_prt[sp].vel[i*vd+j])/(s_cont[sp].Tvec[ptvd+j]*kb/s_prt[sp].mass - 0.0);

               s_prt[sp].vel[i*vd+j] += deltaT/s_prt[sp].mass*(f1+f2+f3);    
               if(std::isnan(s_prt[sp].vel[i*vd+j]) || std::isinf(s_prt[sp].vel[i*vd+j])) std::cout << "\nCoul2:  " << i << "\t" << j << "\t" << velA[j] << "\t" << numom << "\t" << nuE << "\t" << f1 << "\t" << f2 << "\t" << f3 << "\n"; 
               if(std::isnan(s_prt[sp].vel[i*vd+j]) || std::isinf(s_prt[sp].vel[i*vd+j])) std::cout << "\nCoul2:  " << coef1 << "\t" << coef2 << "\t" << s_cont[csp].U[ptvd+j] << "\t" << s_cont[sp].U[ptvd+j] << "\t" << s_cont[sp].Tvec[ptvd+j]; 

               en_temp += velA[j]*velA[j]+s_prt[sp].vel[i*vd+j]*s_prt[sp].vel[i*vd+j];  // Avg before and after
               //en_temp += s_prt[sp].vel[i*vd+j]*s_prt[sp].vel[i*vd+j];  // Avg before and after

               //std::cout << "\nsp,i,j,deltaT,f3,vb,va,numom,nuE:  " << sp << "\t" << i << "\t" << j << "\t" << deltaT << "\t" << f3 << "\t" << velA[j] << "\t" << s_prt[sp].vel[i*vd+j] << "\t" << numom << "\t" << nuE << "\n";
               //std::cout << "\nsp,i,j,deltaT:   " << sp << "\t" << i << "\t" << j << "\t" << deltaT;
               //std::cout << "\nf1,f2,f3,vb,va:  " << f1 << "\t" << f2 << "\t" << f3 << "\t" << velA[j] << "\t" << s_prt[sp].vel[i*vd+j];
               //std::cout << "\nnumom,nuE,dT,dU,dKE:  " << numom << "\t" << nuE << "\t" << (kb*s_cont[sp].T[pt]-kb*s_cont[csp].T[pt]) << "\t" << (s_cont[sp].U[ptvd+j]-s_prt[sp].vel[i*vd+j]) << "\t" << (s_cont[sp].Tvec[ptvd+j]*kb/s_prt[sp].mass - s_cont[sp].U[ptvd+j]*s_cont[sp].U[ptvd+j]) << "\n";
            }
            s_prt[sp].en[i] = fabs(en_temp*(0.25)*s_prt[sp].mass);    
         }
      }
   }

   coulForce = 0;

}


/*void solver::coulombCollisions(std::vector<particles> &s_prt, const mesh &s_msh, std::vector<contnm> &s_cont)
  {
  int i,j,k,sp;
  int ind,ind3,csp;
  long double tempd;
  double dt,numom,nuE;
  double vth,vthsq,vthA,mrd,coef1,coef2;
  double qe4=qe*qe*qe*qe;
  double lnC = 10.0;
  long double vupper = 5.0;
  long double vlower = -5.0;

  long int np;
  long int npg;

  mathFunctions s_mth;
  dt = deltaT;

  std::vector<double> deltav;
  std::vector<double> velA;

  for(i=0;i<s_msh.vecdims;i++) velA.push_back(0.0);

  if(s_prt[0].nsp>1)
  {
  mrd = s_prt[0].mass*s_prt[1].mass/(s_prt[0].mass+s_prt[1].mass);

  coef1 = 8.0*sqrt(pi)*s_prt[0].charge*s_prt[0].charge*s_prt[1].charge*s_prt[1].charge*lnC/(mrd*mrd*16.0*pi*pi*eps0*eps0);
  coef2 = 16.0*sqrt(pi)*s_prt[0].charge*s_prt[0].charge*s_prt[1].charge*s_prt[1].charge*lnC/(s_prt[0].mass*s_prt[1].mass*16.0*pi*pi*eps0*eps0);        

  for(i=0;i<s_cont[0].En.size();i++) 
  {
  tempd = 0.0;
  for(j=0;j<s_msh.vecdims;j++) tempd += (s_cont[0].U[i*s_msh.vecdims+j]-s_cont[1].U[i*s_msh.vecdims+j])*(s_cont[0].U[i*s_msh.vecdims+j]-s_cont[1].U[i*s_msh.vecdims+j]);  
  deltav.push_back(sqrt(tempd)); 
  }

  for(sp=0;sp<s_prt[0].nsp;sp++)
  {
  if(sp==0) csp = 1;
  else if(sp==1) csp = 0;

  for(i=0;i<s_prt[sp].en.size();i++)
  {
  ind = s_prt[sp].pt[i];
  ind3 = 3*ind;

  vthsq =2.0*(kb*s_cont[0].T[ind]/s_prt[0].mass+kb*s_cont[1].T[ind]/s_prt[1].mass);
  vth = sqrt(vthsq);

  numom = coef1*s_cont[csp].N[ind]/(deltav[ind]*deltav[ind]*deltav[ind])*(0.5*sqrt(pi)*erf(deltav[ind]/vth)-(deltav[ind]/vth)*exp(-deltav[ind]*deltav[ind]/vthsq));
  nuE = coef2*s_cont[csp].N[ind]/(vthsq*vth)*exp(-(deltav[ind]*deltav[ind]/vthsq));
  std::cout << "\nnumom1:  " << coef1 << "\t" << numom << std::endl;
  std::cout << "\nnumom2:  " << 1/(deltav[ind]*deltav[ind]*deltav[ind]) << "\t" << (0.5*sqrt(pi)*erf(deltav[ind]/vth)-(deltav[ind]/vth)*exp(-deltav[ind]*deltav[ind]/vthsq));
//std::cout << "\nnuE:  " << coef2 << "\t" << nuE << std::endl;
//std::cout << "\nvth:  " << vth <<  std::endl;

for(j=0;j<s_msh.vecdims;j++)  
{      
velA[j] = s_prt[sp].vel[i*s_msh.vecdims+j];

s_prt[sp].vel[i*s_msh.vecdims+j] += dt/s_prt[sp].mass * (  numom*mrd*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])    -      numom*mrd*mrd/s_prt[sp].mass*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])/(kb*s_cont[sp].Tvec[ind3+j]/s_prt[sp].mass-s_cont[sp].U[ind3+j]*s_cont[sp].U[ind3+j])*(s_cont[sp].U[ind3+j]-s_prt[sp].vel[i*s_msh.vecdims+j])        +          nuE*(kb*s_cont[sp].T[ind]-kb*s_cont[csp].T[ind])/(kb*s_cont[sp].Tvec[ind3+j]/s_prt[sp].mass-s_cont[sp].U[ind3+j]*s_cont[sp].U[ind3+j])*(s_cont[sp].U[ind3+j]-s_prt[sp].vel[i*s_msh.vecdims+j]));       
//if(s_prt[sp].vel[i*s_msh.vecdims+j]>1.0e7 || std::isnan(s_prt[sp].vel[i*s_msh.vecdims+j]) == true ) std::cout << "\n" << s_prt[sp].vel[i*s_msh.vecdims+j] << "\t" << sp << "\t" << i << "\t" << j << "\t" << numom << "\t" << nuE << "\t" << vth << "\t" << dt/s_prt[sp].mass * (  numom*mrd*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j]))   << "\t" << -dt/s_prt[sp].mass*(numom*mrd*mrd/s_prt[sp].mass*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])/(kb*s_cont[sp].Tvec[ind3+j]/s_prt[sp].mass-s_cont[sp].U[ind3+j]*s_cont[sp].U[ind3+j])*(s_cont[sp].U[ind3+j]-s_prt[sp].vel[i*s_msh.vecdims+j])) << "\t" <<  dt/s_prt[sp].mass*(nuE*(kb*s_cont[sp].T[ind]-kb*s_cont[csp].T[ind])/(kb*s_cont[sp].Tvec[ind3+j]/s_prt[sp].mass-s_cont[sp].U[ind3+j]*s_cont[sp].U[ind3+j])*(s_cont[sp].U[ind3+j]-s_prt[sp].vel[i*s_msh.vecdims+j]));       
std::cout << "\n" << velA[j] << "\t" << sp << "\t" << i << "\t" << j << "\t" << numom << "\t" << nuE << "\t" << vth << "\t" << dt/s_prt[sp].mass * (  numom*mrd*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j]))   << "\t" << -dt/s_prt[sp].mass*(numom*mrd*mrd/s_prt[sp].mass*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])*(s_cont[csp].U[ind3+j]-s_cont[sp].U[ind3+j])/(kb*s_cont[sp].Tvec[ind3+j]/s_prt[sp].mass-s_cont[sp].U[ind3+j]*s_cont[sp].U[ind3+j])*(s_cont[sp].U[ind3+j]-velA[j])) << "\t" <<  dt/s_prt[sp].mass*(nuE*(kb*s_cont[sp].T[ind]-kb*s_cont[csp].T[ind])/(kb*s_cont[sp].Tvec[ind3+j]/s_prt[sp].mass-s_cont[sp].U[ind3+j]*s_cont[sp].U[ind3+j])*(s_cont[sp].U[ind3+j]-velA[j]));       
}
}
}
}

/*for(sp=0;sp<s_prt[0].nsp;sp++)
{
coef1 = 4.0/3.0*sqrt(pi)*s_prt[sp].charge*s_prt[sp].charge*s_prt[sp].charge*s_prt[sp].charge*lnC/(s_prt[sp].mass*s_prt[sp].mass)/(4.0*4.0*pi*pi*eps0*eps0);


for(i=0;i<s_prt[sp].en.size();i++)
{
ind = s_prt[sp].pt[i];
ind3 = 3*ind;

numom = coef1*s_cont[sp].N[ind]*pow(kb*s_cont[sp].T[ind]/s_prt[sp].mass,-1.5);

//std::cout << "\nnumom:  " << coef1 << "\t" << numom << std::endl;

vthA = sqrt(2.0*numom*dt*kb*s_cont[sp].T[ind]/s_prt[sp].mass);

//std::cout << "\nvthA:  " << vthA << std::endl;

for(j=0;j<s_msh.vecdims;j++)
{
tempd= rand();
tempd= tempd/RAND_MAX;
tempd= s_mth.erfinv(tempd*erf(vupper)+(1.0-tempd)*erf(vlower));
tempd = sqrt(2.0)*vthA*tempd;
velA[j] = tempd;
//std::cout << "\t " << velA[j];
}

for(j=0;j<s_msh.vecdims;j++) s_prt[sp].vel[i*s_msh.vecdims+j] += -dt*numom*(s_prt[sp].vel[i*s_msh.vecdims+j]-s_cont[sp].U[s_msh.vecdims*ind+j]) + velA[j];

}
}*/



//}

//.....Select random array of particle indices........//

void solver::randomselect(std::vector<int> &randarray, int size, int num)
{
   int i,j,k;
   int ind;
   long double R,P;

   std::vector<int>::iterator check = randarray.begin();

   //srand(time(NULL));

   // std::random_shuffle(i_R.begin(),i_R.end());

   while(randarray.size()<num)
   {
      R = rand(); 
      P = R/RAND_MAX;
      ind = (P*(size-1.0)+0.5);  //..Forced rounding
      //std::cout << std::endl << "ind:  " <<  ind << std::endl; //..PRT

      check = std::find(randarray.begin(),randarray.end(),ind); 

      if(check==randarray.end())
      {
         randarray.push_back(ind);
      }  

   }

}



//.....Seed a single particle, used for ionization.....//

void solver::seedSingleParticle(particles &s_part, const mesh &s_msh, std::vector<double> s_vel, std::vector<double> s_pos, double s_Temp, int s_cell, int s_pt)
{ 

   //.....Maxwellian set based on "Loading and Injection of Maxwellian Distributions in Particle Simulations" by Cartwright, Verboncoeur, and Birdsall....//

   int i,j,k;
   int index;
   double temp;
   double vupper,vlower,vtherm;
   double tempen;

   mathFunctions s_mth;

   vupper = 6.0;
   vlower = -6.0;

   index = s_part.pos.size()/s_msh.meshdims;

   tempen = 0.0;
   vtherm = sqrt(s_Temp*kb/(s_part.mass));

   for(k=0;k<s_msh.meshdims;k++) s_part.pos.push_back(s_pos[k]);

   for(k=0;k<s_msh.vecdims;k++)
   {
      temp=rand();
      temp=temp/RAND_MAX;
      temp=s_mth.erfinv(temp*erf(vupper)+(1-temp)*erf(vlower));
      s_part.vel.push_back(sqrt(2.0)*vtherm*temp+s_vel[k]);
      tempen= tempen + (sqrt(2.0)*vtherm*temp+s_vel[k])*(sqrt(2.0)*vtherm*temp+s_vel[k]);
   }
   s_part.en.push_back((0.5)*tempen*s_part.mass);   //Energy of particle in group, not collection

   s_part.cell.push_back(s_cell);
   s_part.pt.push_back(s_pt);

}

//......Scatter a particle's velocity after collision.......//

void solver::scatterParticle(std::vector<double> &s_vel,double s_vmag,double s_chi, double s_phi)
{
   int i,j,k;
   double inv_sqrt;
   std::vector<double> new_vel;

   for(i=0;i<s_vel.size();i++) new_vel.push_back(0.0);

   inv_sqrt = 1.0/sqrt(s_vel[0]*s_vel[0] + s_vel[1]*s_vel[1]);

   //if(std::isnan(inv_sqrt)==true)  std::cout << std::endl << inv_sqrt  << std::endl ; 

   //if(std::isnan(s_vel[0])==true)  std::cout << "SP Bef..0" << std::endl << s_vel[0]  << std::endl ; 
   //if(std::isnan(s_vel[1])==true)  std::cout << "SP Bef..1" << std::endl << s_vel[1]  << std::endl ; 
   //if(std::isnan(s_vel[2])==true)  std::cout << "SP Bef..2" << std::endl << s_vel[2]  << std::endl ; 
   //if(std::isnan(s_chi)==true)  std::cout << "schi" << std::endl << s_chi  << std::endl ; 
   //if(std::isnan(s_vmag)==true)  std::cout << "svmag" << std::endl << s_vmag  << std::endl ; 
   //if(std::isnan(s_phi)==true)  std::cout << "sphi" << std::endl << s_phi  << std::endl ; 

   new_vel[0] = s_vel[0]*cos(s_chi)    +  s_vel[1]*s_vmag*sin(s_chi)*sin(s_phi)*inv_sqrt   +   s_vel[0]*s_vel[2]*sin(s_chi)*cos(s_phi)*inv_sqrt; 
   new_vel[1] = s_vel[1]*cos(s_chi)  -  s_vel[0]*s_vmag*sin(s_chi)*sin(s_phi)*inv_sqrt     +   s_vel[1]*s_vel[2]*sin(s_chi)*cos(s_phi)*inv_sqrt; 
   new_vel[2] = s_vel[2]*cos(s_chi)  - sin(s_chi)*cos(s_phi)/inv_sqrt; 

   /*if(std::isnan(new_vel[0])==true || std::isnan(new_vel[1])==true || std::isnan(new_vel[2])==true)  
     {
     std::cout << "\nvel" << "\t" << s_vel[0]  << std::endl ; 
     std::cout << "vel" << "\t" << s_vel[1]  << std::endl ; 
     std::cout << "vel" << "\t" << s_vel[2]  << std::endl ; 
     std::cout << "new" << "\t" << new_vel[0]  << std::endl ; 
     std::cout << "new" << "\t" << new_vel[1]  << std::endl ; 
     std::cout << "new" << "\t" << new_vel[2]  << std::endl ; 
     std::cout << "chi" << "\t" << s_chi  << std::endl ; 
     std::cout << "phi" << "\t" << s_phi  << std::endl ; 
     std::cout << "vmag" << "\t" << s_vmag  << std::endl ; 
     }*/

   s_vel[0] = new_vel[0];
   s_vel[1] = new_vel[1];
   s_vel[2] = new_vel[2];

}

//......Find processor with least particles to use as seed processor.......//

void solver::findseedproc(particles &s_part)  //MPI
{ 
   int i,j,k;
   int numprocs,procid; //MPI
   int numpart,s_seed_proc;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   std::vector<int> numpart_array;

   if(procid==0)  for(i=0;i<numprocs;i++) numpart_array.push_back(0);

   numpart = s_part.en.size();

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Gather(&numpart,1,MPI_INT,&numpart_array.front(),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);

   if(procid==0)  s_seed_proc = std::distance(numpart_array.begin(),std::min_element(numpart_array.begin(), numpart_array.end()));

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Bcast(&s_seed_proc,1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);
   s_part.seedproc = s_seed_proc; 
   MPI_Barrier(MPI_COMM_WORLD);

   std::cout << "\nSeed Proc:  " << s_part.name << "\t" << s_part.seedproc << "\t" << procid << "\t" << numpart << std::endl; 
}


//.....Redistribute particles evenly across processors.....//

void solver::redistributeparticles(particles &s_part,int s_vecdims, int s_meshdims)  //MPI
{ 
   int i,j,k;
   int numprocs,procid; //MPI
   int numpart,numpart_total,s_seed_proc,numavg,numrem;

   std::vector<double> s_pos,s_vel,s_en;
   std::vector<int> s_cell,s_pt;

   std::vector<double> s_pos2,s_vel2,s_en2;
   std::vector<int> s_cell2,s_pt2;


   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   std::vector<int> numpart_array;
   std::vector<int> numpos_array,numvel_array,pcnt;
   std::vector<int> disp,meshdisp,vecdisp;
   numpart_total=0;

   for(i=0;i<numprocs;i++) numpart_array.push_back(0);
   for(i=0;i<numprocs;i++) pcnt.push_back(0);

   numpart = s_part.en.size();
   std::cout << "\nStarting particles: " << numpart << std::endl;

   MPI_Barrier(MPI_COMM_WORLD); // block until all processes have reached this routine
   //MPI_Gather(&numpart,1,MPI_INT,&numpart_array.front(),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Allgather(&numpart,1,MPI_INT,&numpart_array.front(),1,MPI_INT,MPI_COMM_WORLD); // Gather data from all taksks and distribute the combined data to all tasks
   MPI_Barrier(MPI_COMM_WORLD); // block until all processes have reached this routine

   MPI_Allreduce(&numpart,&numpart_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD); // Combines values from all processes and distributes the result back to all processes

   for(i=0;i<numpart_total;i++) // iterate over the number of elements in send buffer
   {
      for(j=0;j<s_meshdims;j++) s_pos.push_back(0.0);
      for(j=0;j<s_vecdims;j++) s_vel.push_back(0.0);
      s_en.push_back(0.0);
      s_cell.push_back(0);
      s_pt.push_back(0);

      for(j=0;j<s_meshdims;j++) s_pos2.push_back(0.0);
      for(j=0;j<s_vecdims;j++) s_vel2.push_back(0.0);
      s_en2.push_back(0.0);
      s_cell2.push_back(0);
      s_pt2.push_back(0);
   }


   pcnt[0] = 0;
   for(i=1;i<numprocs;i++) pcnt[i] = pcnt[i-1] + numpart_array[i-1];

   for(i=0;i<numpart;i++)
   {
      for(j=0;j<s_meshdims;j++) s_pos[(i+pcnt[procid])*s_meshdims+j] = s_part.pos[i*s_meshdims+j];
      for(j=0;j<s_vecdims;j++) s_vel[(i+pcnt[procid])*s_vecdims+j] = s_part.vel[i*s_vecdims+j];
      s_en[i+pcnt[procid]] = s_part.en[i];
      s_cell[i+pcnt[procid]] = s_part.cell[i];
      s_pt[i+pcnt[procid]] = s_part.pt[i];
   }

   s_part.pos.clear();
   s_part.vel.clear();
   s_part.en.clear();
   s_part.pt.clear();
   s_part.cell.clear();


   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(&s_pos.front(),&s_pos2.front(),numpart_total*s_meshdims, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_vel.front(),&s_vel2.front(),numpart_total*s_vecdims, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_en.front(),&s_en2.front(),numpart_total, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_cell.front(),&s_cell2.front(),numpart_total, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_pt.front(),&s_pt2.front(),numpart_total, MPI_INT,MPI_SUM,MPI_COMM_WORLD);

   s_pos.clear();
   s_vel.clear();
   s_en.clear();
   s_pt.clear();
   s_cell.clear();

   numavg = numpart_total/numprocs;
   numrem = numpart_total - numavg*numprocs;

   numpart_array[0] = numavg+numrem;
   for(i=1;i<numprocs;i++) numpart_array[i] = numavg;

   pcnt[0] = 0;
   for(i=1;i<numprocs;i++) pcnt[i] = pcnt[i-1] + numpart_array[i-1];

   for(i=0;i<numpart_array[procid];i++)
   {
      for(j=0;j<s_meshdims;j++) s_part.pos.push_back(s_pos2[(pcnt[procid]+i)*s_meshdims+j]);
      for(j=0;j<s_vecdims;j++) s_part.vel.push_back(s_vel2[(pcnt[procid]+i)*s_vecdims+j]);
      s_part.en.push_back(s_en2[pcnt[procid]+i]);
      s_part.cell.push_back(s_cell2[pcnt[procid]+i]);
      s_part.pt.push_back(s_pt2[pcnt[procid]+i]);
   } 

   s_pos2.clear();
   s_vel2.clear();
   s_en2.clear();
   s_pt2.clear();
   s_cell2.clear();

   numpart = s_part.en.size();
   std::cout << "\nProc:  " << procid << "\t" << s_part.name <<"\t" << numpart;

   MPI_Barrier(MPI_COMM_WORLD);
}


void solver::updateNeutralBackground(contnm &s_neut, const mesh &s_msh, double s_dt)
{
   int i,j,k;
   int numprocs,procid; //MPI
   double V=sqrt(5.0/3.0*s_neut.R*s_neut.TL);
   double C=s_dt*V/s_msh.deltax; 
   double invdV = 1/(s_msh.deltax*s_msh.parea[0]);
   std::vector<double> ion_all;
   std::vector<double> N_new;

   //std::cout << "\tUpdating Neutral BG..." << V << "\t" << C;

   for(i=0;i<s_neut.N.size();i++) N_new.push_back(0);

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   for(i=0;i<s_neut.ioncount.size();i++) ion_all.push_back(0.0);
   MPI_Allreduce(&s_neut.ioncount.front(),&ion_all.front(),s_neut.ioncount.size(), MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);

   //...Upwind Scheme;

   N_new[0] = s_neut.N[0] - C*(s_neut.N[0]-s_neut.NL) - ion_all[0]*invdV;  
   for(i=1;i<s_neut.N.size();i++)  
   {
      N_new[i] = s_neut.N[i] - C*(s_neut.N[i]-s_neut.N[i-1]) - ion_all[i]*invdV;  
      //std::cout << " \t" << ion_all[i]*invdV;  
   }

   for(i=0;i<s_neut.ioncount.size();i++) s_neut.ioncount[i] = 0;
   for(i=0;i<s_neut.N.size();i++) s_neut.N[i] = N_new[i];

}






/*










//...........Euler Solver from Class Project(Need to fix)...........//////


double solver::errorCheck(const flow &s_cont,const mesh &s_msh)
{
int i;
double val1,val2,val3,error;

error = 0.0;  

for(i=2;i<(s_msh.cmesh.size()-2);i++)
{
val1 = s_cont.rho[i];
val2 = (1.0+(1.0e-4)*pow(cos(3.14158*s_msh.cmesh[i]/2.0),8.0));

error = error + fabs(val1-val2);

//      std::cout << std::setprecision(10) <<  i << "\t" << "  v1:  " << val1 << "  v2:   " << val2 << "  error:  " << error << std::endl;
}

error = error/(s_msh.cmesh.size()-4);

return error;
}


void solver::updateEulerFlow(flow &s_cont,const mesh &s_msh,double s_dt,int s_tflag,int s_sflag)
{
int i,j,k;
std::vector<double> rhonew, Unew, pnew, Ennew;
std::vector<double> rhoRK, momRK, pRK, EnRK;
std::vector<double> drhoRK, dmomRK, dEnRK;
std::vector<double> fp,fc,fm,Fp,Fm;
std::vector<double> fLp,fRp,fLm,fRm;
std::vector<double> k_RK;
double dx,alpha;
double mfact[4] = {0.5,0.5,1.0,0.0};
double kfact[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};

double rhopp,rhop,rhoc,rhoci,rhom,rhomm;
double mompp,momp,momc,momci,momm,mommm;
double ppp,pp,pc,pci,pm,pmm;
double Enpp,Enp,Enc,Enci,Enm,Enmm;
double rhoqLp,rhoqRp,rhoqLm,rhoqRm;
double momqLp,momqRp,momqLm,momqRm;
double pqLp,pqRp,pqLm,pqRm;
double EnqLp,EnqRp,EnqLm,EnqRm;

for(i=0;i<3;i++) fp.push_back(0.0);   
for(i=0;i<3;i++) fc.push_back(0.0);   
for(i=0;i<3;i++) fm.push_back(0.0);   
for(i=0;i<3;i++) Fp.push_back(0.0);   
for(i=0;i<3;i++) Fm.push_back(0.0);   
for(i=0;i<3;i++) fLp.push_back(0.0);   
for(i=0;i<3;i++) fRp.push_back(0.0);   
for(i=0;i<3;i++) fLm.push_back(0.0);   
for(i=0;i<3;i++) fRm.push_back(0.0);   
for(i=0;i<4;i++) k_RK.push_back(0.0);   


for(i=0;i<s_cont.rho.size();i++) rhonew.push_back(s_cont.rho[i]);
for(i=0;i<s_cont.U.size();i++)  Unew.push_back(s_cont.U[i]);
for(i=0;i<s_cont.p.size();i++)  pnew.push_back(s_cont.p[i]);
for(i=0;i<s_cont.En.size();i++) Ennew.push_back(s_cont.En[i]);
for(i=0;i<s_cont.rho.size();i++) rhoRK.push_back(s_cont.rho[i]);
for(i=0;i<s_cont.U.size();i++)  momRK.push_back(s_cont.rho[i]*s_cont.U[i]);
for(i=0;i<s_cont.p.size();i++)  pRK.push_back(s_cont.p[i]);
for(i=0;i<s_cont.En.size();i++) EnRK.push_back(s_cont.En[i]);
for(i=0;i<s_cont.rho.size();i++) drhoRK.push_back(0.0);
for(i=0;i<s_cont.U.size();i++)  dmomRK.push_back(0.0);
for(i=0;i<s_cont.En.size();i++) dEnRK.push_back(0.0);



//............RK4 with Central Difference................//


if(s_tflag==1 && s_sflag==0) 
{
   dx =  s_msh.cmesh[1]-s_msh.cmesh[0];


   for(j=0;j<4;j++)
   {

      for(i=2;i<(s_msh.cmesh.size()-2);i++)
      {
         rhop = rhoRK[i+1];
         rhoc = rhoRK[i];
         rhom = rhoRK[i-1];

         momp = momRK[i+1];
         momc = momRK[i];
         momm = momRK[i-1];

         pp = pRK[i+1];
         pc = pRK[i];
         pm = pRK[i-1];

         Enp = EnRK[i+1];
         Enc = EnRK[i];
         Enm = EnRK[i-1];

         flux1D(fp,rhop,momp,pp,Enp);
         flux1D(fm,rhom,momm,pm,Enm);

         for(k=0;k<3;k++) k_RK[k] = s_dt*(fm[k]-fp[k])/(2.0*dx) ;

         drhoRK[i] = k_RK[0];
         dmomRK[i] = k_RK[1];
         dEnRK[i]  = k_RK[2];

      }

      for(i=2;i<(s_msh.cmesh.size()-2);i++)
      {
         rhoRK[i] = s_cont.rho[i] + mfact[j]*drhoRK[i];
         momRK[i] =  s_cont.rho[i]*s_cont.U[i] + mfact[j]*dmomRK[i];
         EnRK[i] = s_cont.En[i] + mfact[j]*dEnRK[i];
         pRK[i] = (s_cont.gam-1)*(EnRK[i] - 0.5*momRK[i]*momRK[i]/rhoRK[i]);

         rhonew[i] =  rhonew[i] + kfact[j]*drhoRK[i]; 
         Unew[i] =  Unew[i] + kfact[j]*dmomRK[i]/rhonew[i]; 
         Ennew[i] =  Ennew[i]+ kfact[j]*dEnRK[i]; 
         pnew[i] = (s_cont.gam-1.0)*(Ennew[i] - 0.5*rhonew[i]*Unew[i]*Unew[i]);

         //std::cout <<"\n" << rhonew[i] << "\t" << Unew[i] << "\t" << Ennew[i] << "\t" << pnew[i] << "\n"; 
      }

   }

   for(i=2;i<(s_msh.cmesh.size()-2);i++)
   {
      s_cont.rho[i] = rhonew[i];
      s_cont.U[i] = Unew[i];
      s_cont.En[i] = Ennew[i];
      s_cont.p[i] = pnew[i];
   }
}

//.................RK4 with Lax Friedrich.................//

if(s_tflag==1 && s_sflag==1) 
{
   dx =  s_msh.cmesh[1]-s_msh.cmesh[0];
   alpha = dx/s_dt;

   for(j=0;j<4;j++)
   {

      for(i=2;i<(s_msh.cmesh.size()-2);i++)
      {
         rhop = rhoRK[i+1];
         rhoc = rhoRK[i];
         rhom = rhoRK[i-1];

         momp = momRK[i+1];
         momc = momRK[i];
         momm = momRK[i-1];

         pp = pRK[i+1];
         pc = pRK[i];
         pm = pRK[i-1];

         Enp = EnRK[i+1];
         Enc = EnRK[i];
         Enm = EnRK[i-1];

         flux1D(fp,rhop,momp,pp,Enp);
         flux1D(fc,rhoc,momc,pc,Enc);
         flux1D(fm,rhom,momm,pm,Enm);

         flux1DLF(Fp,fp,fc,rhoc,rhop,momc,momp,Enc,Enp,alpha);
         flux1DLF(Fm,fc,fm,rhom,rhoc,momm,momc,Enm,Enc,alpha);

         for(k=0;k<3;k++) k_RK[k] = s_dt*(Fm[k]-Fp[k])/(dx) ;

         drhoRK[i] = k_RK[0];
         dmomRK[i] = k_RK[1];
         dEnRK[i]  = k_RK[2];

      }

      for(i=2;i<(s_msh.cmesh.size()-2);i++)
      {
         rhoRK[i] = s_cont.rho[i] + mfact[j]*drhoRK[i];
         momRK[i] =  s_cont.rho[i]*s_cont.U[i] + mfact[j]*dmomRK[i];
         EnRK[i] = s_cont.En[i] + mfact[j]*dEnRK[i];
         pRK[i] = (s_cont.gam-1)*(EnRK[i] - 0.5*momRK[i]*momRK[i]/rhoRK[i]);

         rhonew[i] =  rhonew[i] + kfact[j]*drhoRK[i]; 
         Unew[i] =  Unew[i] + kfact[j]*dmomRK[i]/rhonew[i]; 
         Ennew[i] =  Ennew[i]+ kfact[j]*dEnRK[i]; 
         pnew[i] = (s_cont.gam-1.0)*(Ennew[i] - 0.5*rhonew[i]*Unew[i]*Unew[i]);

      }

   }

   for(i=2;i<(s_msh.cmesh.size()-2);i++)
   {
      s_cont.rho[i] = rhonew[i];
      s_cont.U[i] = Unew[i];
      s_cont.En[i] = Ennew[i];
      s_cont.p[i] = pnew[i];
   }
}

//.................RK4 with MUSCL.................//

if(s_tflag==1 && s_sflag==2) 
{
   dx =  s_msh.cmesh[1]-s_msh.cmesh[0];
   alpha = dx/s_dt;

   for(j=0;j<4;j++)
   {

      for(i=2;i<(s_msh.cmesh.size()-2);i++)
      {
         rhopp = rhoRK[i+2];
         rhop = rhoRK[i+1];
         rhoc = rhoRK[i];
         rhom = rhoRK[i-1];
         rhomm = rhoRK[i-2];

         mompp = momRK[i+2];
         momp = momRK[i+1];
         momc = momRK[i];
         momm = momRK[i-1];
         mommm = momRK[i-2];

         ppp = pRK[i+2];
         pp = pRK[i+1];
         pc = pRK[i];
         pm = pRK[i-1];
         pmm = pRK[i-2];

         Enpp = EnRK[i+2];
         Enp = EnRK[i+1];
         Enc = EnRK[i];
         Enm = EnRK[i-1];
         Enmm = EnRK[i-2];

         rhoqLp = queueL(rhop,rhoc,rhom);
         rhoqRp = queueR(rhopp,rhop,rhoc);
         rhoqLm = queueL(rhoc,rhom,rhomm);
         rhoqRm = queueR(rhop,rhoc,rhom);

         momqLp = queueL(momp,momc,momm);
         momqRp = queueR(mompp,momp,momc);
         momqLm = queueL(momc,momm,mommm);
         momqRm = queueR(momp,momc,momm);

         EnqLp = queueL(Enp,Enc,Enm);
         EnqRp = queueR(Enpp,Enp,Enc);
         EnqLm = queueL(Enc,Enm,Enmm);
         EnqRm = queueR(Enp,Enc,Enm);

         pqLp = queueL(pp,pc,pm);
         pqRp = queueR(ppp,pp,pc);
         pqLm = queueL(pc,pm,pmm);
         pqRm = queueR(pp,pc,pm);

         flux1D(fRp,rhoqRp,momqRp,pqRp,EnqRp);
         flux1D(fLp,rhoqLp,momqLp,pqLp,EnqLp);
         flux1D(fRm,rhoqRm,momqRm,pqRm,EnqRm);
         flux1D(fLm,rhoqLm,momqLm,pqLm,EnqLm);

         flux1DLF(Fp,fRp,fLp,rhoqLp,rhoqRp,momqLp,momqRp,EnqLp,EnqRp,alpha);
         flux1DLF(Fm,fRm,fLm,rhoqLm,rhoqRm,momqLm,momqRm,EnqLm,EnqRm,alpha);

         for(k=0;k<3;k++) k_RK[k] = s_dt*(Fm[k]-Fp[k])/(dx) ;

         drhoRK[i] = k_RK[0];
         dmomRK[i] = k_RK[1];
         dEnRK[i]  = k_RK[2];

      }

      for(i=2;i<(s_msh.cmesh.size()-2);i++)
      {
         rhoRK[i] = s_cont.rho[i] + mfact[j]*drhoRK[i];
         momRK[i] =  s_cont.rho[i]*s_cont.U[i] + mfact[j]*dmomRK[i];
         EnRK[i] = s_cont.En[i] + mfact[j]*dEnRK[i];
         pRK[i] = (s_cont.gam-1)*(EnRK[i] - 0.5*momRK[i]*momRK[i]/rhoRK[i]);

         rhonew[i] =  rhonew[i] + kfact[j]*drhoRK[i]; 
         Unew[i] =  Unew[i] + kfact[j]*dmomRK[i]/rhonew[i]; 
         Ennew[i] =  Ennew[i]+ kfact[j]*dEnRK[i]; 
         pnew[i] = (s_cont.gam-1.0)*(Ennew[i] - 0.5*rhonew[i]*Unew[i]*Unew[i]);

      }

   }

   for(i=2;i<(s_msh.cmesh.size()-2);i++)
   {
      s_cont.rho[i] = rhonew[i];
      s_cont.U[i] = Unew[i];
      s_cont.En[i] = Ennew[i];
      s_cont.p[i] = pnew[i];
   }

}

}


void solver::updateEulerFlow2D(flow &s_cont,const mesh &s_msh,double s_dt,int s_tflag,int s_sflag)
{
   int i,j,k,n,ind;
   int numcells;
   std::vector<double> rhonew, Unew, pnew, Ennew;
   std::vector<double> rhoRK, xmomRK, ymomRK,pRK, EnRK;
   std::vector<double> drhoRK, dxmomRK,dymomRK, dEnRK;
   std::vector<double> fr,fcx,fcy,fl,fb,ft,Fr,Fl,Fb,Ft;
   std::vector<double> fLr,fRr,fLl,fRl;
   std::vector<double> fLt,fRt,fLb,fRb;
   std::vector<double> k_RK;
   std::vector<int> neigh,neigh2;   

   double dx,dy,alpha;
   double mfact[4] = {0.5,0.5,1.0,0.0};
   double kfact[4] = {1.0/6.0,1.0/3.0,1.0/3.0,1.0/6.0};

   double rhorr,rhor,rhoc,rhoci,rhol,rholl;
   double rhott,rhot,rhob,rhobb;
   double xmomrr,xmomr,xmomc,xmomci,xmoml,xmomll;
   double xmomtt,xmomt,xmomb,xmombb;
   double ymomrr,ymomr,ymomc,ymomci,ymoml,ymomll;
   double ymomtt,ymomt,ymomb,ymombb;
   double prr,pr,pc,pci,pl,pll;
   double ptt,pt,pb,pbb;
   double Enrr,Enr,Enc,Enci,Enl,Enll;
   double Entt,Ent,Enb,Enbb;
   double rhoqRr,rhoqRl,rhoqLr,rhoqLl;
   double rhoqRt,rhoqRb,rhoqLt,rhoqLb;
   double xmomqRr,xmomqRl,xmomqLr,xmomqLl;
   double xmomqRt,xmomqRb,xmomqLt,xmomqLb;
   double ymomqRr,ymomqRl,ymomqLr,ymomqLl;
   double ymomqRt,ymomqRb,ymomqLt,ymomqLb;
   double pqRr,pqRl,pqLr,pqLl;
   double pqRt,pqRb,pqLt,pqLb;
   double EnqRr,EnqRl,EnqLr,EnqLl;
   double EnqRt,EnqRb,EnqLt,EnqLb;

   for(i=0;i<3;i++) fr.push_back(0.0);   
   for(i=0;i<3;i++) fcy.push_back(0.0);   
   for(i=0;i<3;i++) fcx.push_back(0.0);   
   for(i=0;i<3;i++) fl.push_back(0.0);   
   for(i=0;i<3;i++) ft.push_back(0.0);   
   for(i=0;i<3;i++) fb.push_back(0.0);   
   for(i=0;i<3;i++) Fr.push_back(0.0);   
   for(i=0;i<3;i++) Fl.push_back(0.0);   
   for(i=0;i<3;i++) Ft.push_back(0.0);   
   for(i=0;i<3;i++) Fb.push_back(0.0);   
   for(i=0;i<3;i++) fLr.push_back(0.0);   
   for(i=0;i<3;i++) fRr.push_back(0.0);   
   for(i=0;i<3;i++) fLl.push_back(0.0);   
   for(i=0;i<3;i++) fRl.push_back(0.0);   
   for(i=0;i<3;i++) fLt.push_back(0.0);   
   for(i=0;i<3;i++) fRt.push_back(0.0);   
   for(i=0;i<3;i++) fLb.push_back(0.0);   
   for(i=0;i<3;i++) fRb.push_back(0.0);   


   for(i=0;i<5;i++) k_RK.push_back(0.0);   
   for(i=0;i<4;i++) neigh.push_back(0);   
   for(i=0;i<4;i++) neigh2.push_back(0);   


   for(i=0;i<s_cont.rho.size();i++) rhonew.push_back(s_cont.rho[i]);
   for(i=0;i<s_cont.U.size();i++)  Unew.push_back(s_cont.U[i]);
   for(i=0;i<s_cont.p.size();i++)  pnew.push_back(s_cont.p[i]);
   for(i=0;i<s_cont.En.size();i++) Ennew.push_back(s_cont.En[i]);
   for(i=0;i<s_cont.rho.size();i++) rhoRK.push_back(s_cont.rho[i]);
   for(i=0;i<(s_cont.U.size()/2);i++)  xmomRK.push_back(s_cont.rho[i]*s_cont.U[s_msh.vecdims*i]);
   for(i=0;i<(s_cont.U.size()/2);i++)  ymomRK.push_back(s_cont.rho[i]*s_cont.U[s_msh.vecdims*i+1]);
   for(i=0;i<s_cont.p.size();i++)  pRK.push_back(s_cont.p[i]);
   for(i=0;i<s_cont.En.size();i++) EnRK.push_back(s_cont.En[i]);
   for(i=0;i<s_cont.rho.size();i++) drhoRK.push_back(0.0);
   for(i=0;i<(s_cont.U.size()/2);i++)  dxmomRK.push_back(0.0);
   for(i=0;i<(s_cont.U.size()/2);i++)  dymomRK.push_back(0.0);
   for(i=0;i<s_cont.En.size();i++) dEnRK.push_back(0.0);



   //............RK4 with Central Difference................//


   if(s_tflag==1 && s_sflag==0) 
   {
      dx =  s_msh.cmesh[s_msh.meshdims]-s_msh.cmesh[0];
      dy =  dx;
      numcells = (s_msh.numpoints[0]-1)*(s_msh.numpoints[1]-1);

      //std::cout << std::endl << "dx:   " << dx << "dy:  " << dy << "dt:  " << s_dt << std::endl;

      for(j=0;j<4;j++)
      {

         for(i=0;i<(s_msh.numpoints[1]-1);i++)
         {
            for(n=0;n<(s_msh.numpoints[0]-1);n++)
            {
               ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;
               //std::cout << "Index: \t" << ind << std::endl;

               neigh = s_msh.cneighofc(ind);

               //std::cout << "Neighbor:  " << neigh[0] << "\t" << neigh[1] << "\t" << std::endl;

               rhoc = rhoRK[ind];
               rhol = rhoRK[neigh[0]];
               rhor = rhoRK[neigh[1]];
               rhob = rhoRK[neigh[2]];
               rhot = rhoRK[neigh[3]];

               xmomc = xmomRK[ind];
               xmoml = xmomRK[neigh[0]];
               xmomr = xmomRK[neigh[1]];
               xmomb = xmomRK[neigh[2]];
               xmomt = xmomRK[neigh[3]];

               ymomc = ymomRK[ind];
               ymoml = ymomRK[neigh[0]];
               ymomr = ymomRK[neigh[1]];
               ymomb = ymomRK[neigh[2]];
               ymomt = ymomRK[neigh[3]];


               pc = pRK[ind];
               pl = pRK[neigh[0]];
               pr = pRK[neigh[1]];
               pb = pRK[neigh[2]];
               pt = pRK[neigh[3]];

               Enc = EnRK[ind];
               Enl = EnRK[neigh[0]];
               Enr = EnRK[neigh[1]];
               Enb = EnRK[neigh[2]];
               Ent = EnRK[neigh[3]];

               xflux2D(fl,rhol,xmoml,ymoml,pl,Enl);
               xflux2D(fr,rhor,xmomr,ymomr,pr,Enr);
               yflux2D(ft,rhot,xmomt,ymomt,pt,Ent);
               yflux2D(fb,rhob,xmomb,ymomb,pb,Enb);

               for(k=0;k<4;k++) k_RK[k] = s_dt*(((fl[k]-fr[k])/(2.0*dx)+(fb[k]-ft[k])/(2.0*dy))) ;

               drhoRK[ind]  = k_RK[0];
               dxmomRK[ind] = k_RK[1];
               dymomRK[ind] = k_RK[2];
               dEnRK[ind]   = k_RK[3];

            }
         }

         for(i=0;i<(s_msh.numpoints[1]-1);i++)
         {
            for(n=0;n<(s_msh.numpoints[0]-1);n++)
            {
               ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;

               rhoRK[ind] = s_cont.rho[ind] + mfact[j]*drhoRK[ind];
               xmomRK[ind] =  s_cont.rho[ind]*s_cont.U[ind*s_msh.vecdims] + mfact[j]*dxmomRK[ind];
               ymomRK[ind] =  s_cont.rho[ind]*s_cont.U[ind*s_msh.vecdims+1] + mfact[j]*dymomRK[ind];
               EnRK[ind] = s_cont.En[ind] + mfact[j]*dEnRK[ind];
               pRK[ind] = (s_cont.gam-1.0)*(EnRK[ind] - 0.5*(xmomRK[ind]*xmomRK[ind]+ymomRK[ind]*ymomRK[ind])/rhoRK[ind]);

               rhonew[ind] =  rhonew[ind] + kfact[j]*drhoRK[ind]; 
               Unew[ind*s_msh.vecdims] =  Unew[ind*s_msh.vecdims] + kfact[j]*dxmomRK[ind]/rhonew[ind]; 
               Unew[ind*s_msh.vecdims+1] =  Unew[ind*s_msh.vecdims+1] + kfact[j]*dymomRK[ind]/rhonew[ind]; 
               Ennew[ind] =  Ennew[ind]+ kfact[j]*dEnRK[ind]; 
               pnew[ind] = (s_cont.gam-1.0)*(Ennew[ind] - 0.5*rhonew[ind]*(Unew[ind*s_msh.vecdims]*Unew[ind*s_msh.vecdims]+Unew[ind*s_msh.vecdims+1]*Unew[ind*s_msh.vecdims+1]));
            }
            //std::cout <<"\n" << rhonew[i] << "\t" << Unew[i] << "\t" << Ennew[i] << "\t" << pnew[i] << "\n"; 
         }
      }

      for(i=0;i<(s_msh.numpoints[1]-1);i++)
      {
         for(n=0;n<(s_msh.numpoints[0]-1);n++)
         {
            ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;        
            s_cont.rho[ind] = rhonew[ind];
            s_cont.U[ind*s_msh.vecdims] = Unew[ind*s_msh.vecdims];
            s_cont.U[ind*s_msh.vecdims+1] = Unew[ind*s_msh.vecdims+1];
            s_cont.En[ind] = Ennew[ind];
            s_cont.p[ind] = pnew[ind];
         }
      }
   }

   //.................RK4 with Lax Friedrich.................//

   if(s_tflag==1 && s_sflag==1) 
   {
      dx =  s_msh.cmesh[s_msh.meshdims]-s_msh.cmesh[0];
      dy =  dx;
      numcells = (s_msh.numpoints[0]-1)*(s_msh.numpoints[1]-1);
      alpha = dx/s_dt;

      //std::cout << std::endl << "dx:   " << dx << "dy:  " << dy << "dt:  " << s_dt << std::endl;

      for(j=0;j<4;j++)
      {

         for(i=0;i<(s_msh.numpoints[1]-1);i++)
         {
            for(n=0;n<(s_msh.numpoints[0]-1);n++)
            {
               ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;
               //std::cout << "Index: \t" << ind << std::endl;

               neigh = s_msh.cneighofc(ind);

               //std::cout << "Neighbor:  " << neigh[0] << "\t" << neigh[1] << "\t" << std::endl;

               rhoc = rhoRK[ind];
               rhol = rhoRK[neigh[0]];
               rhor = rhoRK[neigh[1]];
               rhob = rhoRK[neigh[2]];
               rhot = rhoRK[neigh[3]];

               xmomc = xmomRK[ind];
               xmoml = xmomRK[neigh[0]];
               xmomr = xmomRK[neigh[1]];
               xmomb = xmomRK[neigh[2]];
               xmomt = xmomRK[neigh[3]];

               ymomc = ymomRK[ind];
               ymoml = ymomRK[neigh[0]];
               ymomr = ymomRK[neigh[1]];
               ymomb = ymomRK[neigh[2]];
               ymomt = ymomRK[neigh[3]];

               pc = pRK[ind];
               pl = pRK[neigh[0]];
               pr = pRK[neigh[1]];
               pb = pRK[neigh[2]];
               pt = pRK[neigh[3]];

               Enc = EnRK[ind];
               Enl = EnRK[neigh[0]];
               Enr = EnRK[neigh[1]];
               Enb = EnRK[neigh[2]];
               Ent = EnRK[neigh[3]];

               xflux2D(fl,rhol,xmoml,ymoml,pl,Enl);
               xflux2D(fcx,rhoc,xmomc,ymomc,pc,Enc);
               xflux2D(fr,rhor,xmomr,ymomr,pr,Enr);
               yflux2D(ft,rhot,xmomt,ymomt,pt,Ent);
               yflux2D(fb,rhob,xmomb,ymomb,pb,Enb);
               yflux2D(fcy,rhoc,xmomc,ymomc,pc,Enc);

               flux2DLF(Fl,fcx,fl,rhol,rhoc,xmoml,xmomc,ymoml,ymomc,Enl,Enc,alpha);
               flux2DLF(Fr,fr,fcx,rhoc,rhor,xmomc,xmomr,ymomc,ymomr,Enc,Enr,alpha);
               flux2DLF(Fb,fcy,fb,rhob,rhoc,xmomb,xmomc,ymomb,ymomc,Enb,Enc,alpha);
               flux2DLF(Ft,ft,fcy,rhoc,rhot,xmomc,xmomt,ymomc,ymomt,Enc,Ent,alpha);


               for(k=0;k<4;k++) k_RK[k] = s_dt*(((Fl[k]-Fr[k])/(dx)+(Fb[k]-Ft[k])/(dy))) ;

               drhoRK[ind]  = k_RK[0];
               dxmomRK[ind] = k_RK[1];
               dymomRK[ind] = k_RK[2];
               dEnRK[ind]   = k_RK[3];
            }
         }

         for(i=0;i<(s_msh.numpoints[1]-1);i++)
         {
            for(n=0;n<(s_msh.numpoints[0]-1);n++)
            {
               ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;

               rhoRK[ind] = s_cont.rho[ind] + mfact[j]*drhoRK[ind];
               xmomRK[ind] =  s_cont.rho[ind]*s_cont.U[ind*s_msh.vecdims] + mfact[j]*dxmomRK[ind];
               ymomRK[ind] =  s_cont.rho[ind]*s_cont.U[ind*s_msh.vecdims+1] + mfact[j]*dymomRK[ind];
               EnRK[ind] = s_cont.En[ind] + mfact[j]*dEnRK[ind];
               pRK[ind] = (s_cont.gam-1.0)*(EnRK[ind] - 0.5*(xmomRK[ind]*xmomRK[ind]+ymomRK[ind]*ymomRK[ind])/rhoRK[ind]);

               rhonew[ind] =  rhonew[ind] + kfact[j]*drhoRK[ind]; 
               Unew[ind*s_msh.vecdims] =  Unew[ind*s_msh.vecdims] + kfact[j]*dxmomRK[ind]/rhonew[ind]; 
               Unew[ind*s_msh.vecdims+1] =  Unew[ind*s_msh.vecdims+1] + kfact[j]*dymomRK[ind]/rhonew[ind]; 
               Ennew[ind] =  Ennew[ind]+ kfact[j]*dEnRK[ind]; 
               pnew[ind] = (s_cont.gam-1.0)*(Ennew[ind] - 0.5*rhonew[ind]*(Unew[ind*s_msh.vecdims]*Unew[ind*s_msh.vecdims]+Unew[ind*s_msh.vecdims+1]*Unew[ind*s_msh.vecdims+1]));
            }
            //std::cout <<"\n" << rhonew[i] << "\t" << Unew[i] << "\t" << Ennew[i] << "\t" << pnew[i] << "\n"; 
         }
      }

      for(i=0;i<(s_msh.numpoints[1]-1);i++)
      {
         for(n=0;n<(s_msh.numpoints[0]-1);n++)
         {
            ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;        
            s_cont.rho[ind] = rhonew[ind];
            s_cont.U[ind*s_msh.vecdims] = Unew[ind*s_msh.vecdims];
            s_cont.U[ind*s_msh.vecdims+1] = Unew[ind*s_msh.vecdims+1];
            s_cont.En[ind] = Ennew[ind];
            s_cont.p[ind] = pnew[ind];
         }
      }
   }


   //.................RK4 with MUSCL.................//

   if(s_tflag==1 && s_sflag==2) 
   {
      dx =  s_msh.cmesh[s_msh.meshdims]-s_msh.cmesh[0];
      dy =  dx;
      numcells = (s_msh.numpoints[0]-1)*(s_msh.numpoints[1]-1);
      alpha = dx/s_dt;

      //std::cout << std::endl << "dx:   " << dx << "dy:  " << dy << "dt:  " << s_dt << std::endl;

      for(j=0;j<4;j++)
      {

         for(i=0;i<(s_msh.numpoints[1]-1);i++)
         {

            for(n=0;n<(s_msh.numpoints[0]-1);n++)
            {
               ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;

               neigh = s_msh.cneighofc(ind);
               for(k=0;k<4;k++) neigh2[k] = ind+2*(neigh[k]-ind);

               rhoc  = rhoRK[ind];
               rhol  = rhoRK[neigh[0]];
               rholl = rhoRK[neigh2[0]];
               rhor  = rhoRK[neigh[1]];
               rhorr = rhoRK[neigh2[1]];
               rhob  = rhoRK[neigh[2]];
               rhobb = rhoRK[neigh2[2]];
               rhot  = rhoRK[neigh[3]];
               rhott = rhoRK[neigh2[3]];

               xmomc  = xmomRK[ind];
               xmoml  = xmomRK[neigh[0]];
               xmomll = xmomRK[neigh2[0]];
               xmomr  = xmomRK[neigh[1]];
               xmomrr = xmomRK[neigh2[1]];
               xmomb  = xmomRK[neigh[2]];
               xmombb = xmomRK[neigh2[2]];
               xmomt  = xmomRK[neigh[3]];
               xmomtt = xmomRK[neigh2[3]];

               ymomc  = ymomRK[ind];
               ymoml  = ymomRK[neigh[0]];
               ymomll = ymomRK[neigh2[0]];
               ymomr  = ymomRK[neigh[1]];
               ymomrr = ymomRK[neigh2[1]];
               ymomb  = ymomRK[neigh[2]];
               ymombb = ymomRK[neigh2[2]];
               ymomt  = ymomRK[neigh[3]];
               ymomtt = ymomRK[neigh2[3]];

               pc  = pRK[ind];
               pl  = pRK[neigh[0]];
               pll = pRK[neigh2[0]];
               pr  = pRK[neigh[1]];
               prr = pRK[neigh2[1]];
               pb  = pRK[neigh[2]];
               pbb = pRK[neigh2[2]];
               pt  = pRK[neigh[3]];
               ptt = pRK[neigh2[3]];

               Enc  = EnRK[ind];
               Enl  = EnRK[neigh[0]];
               Enll = EnRK[neigh2[0]];
               Enr  = EnRK[neigh[1]];
               Enrr = EnRK[neigh2[1]];
               Enb  = EnRK[neigh[2]];
               Enbb = EnRK[neigh2[2]];
               Ent  = EnRK[neigh[3]];
               Entt = EnRK[neigh2[3]];

               rhoqLr = queueL(rhor,rhoc,rhol);
               rhoqRr = queueR(rhorr,rhor,rhoc);
               rhoqLl = queueL(rhoc,rhol,rholl);
               rhoqRl = queueR(rhor,rhoc,rhol);
               rhoqLt = queueL(rhot,rhoc,rhob);
               rhoqRt = queueR(rhott,rhot,rhoc);
               rhoqLb = queueL(rhoc,rhob,rhobb);
               rhoqRb = queueR(rhot,rhoc,rhob);

               xmomqLr = queueL(xmomr,xmomc,xmoml);
               xmomqRr = queueR(xmomrr,xmomr,xmomc);
               xmomqLl = queueL(xmomc,xmoml,xmomll);
               xmomqRl = queueR(xmomr,xmomc,xmoml);
               xmomqLt = queueL(xmomt,xmomc,xmomb);
               xmomqRt = queueR(xmomtt,xmomt,xmomc);
               xmomqLb = queueL(xmomc,xmomb,xmombb);
               xmomqRb = queueR(xmomt,xmomc,xmomb);

               ymomqLr = queueL(ymomr,ymomc,ymoml);
               ymomqRr = queueR(ymomrr,ymomr,ymomc);
               ymomqLl = queueL(ymomc,ymoml,ymomll);
               ymomqRl = queueR(ymomr,ymomc,ymoml);
               ymomqLt = queueL(ymomt,ymomc,ymomb);
               ymomqRt = queueR(ymomtt,ymomt,ymomc);
               ymomqLb = queueL(ymomc,ymomb,ymombb);
               ymomqRb = queueR(ymomt,ymomc,ymomb);

               EnqLr = queueL(Enr,Enc,Enl);
               EnqRr = queueR(Enrr,Enr,Enc);
               EnqLl = queueL(Enc,Enl,Enll);
               EnqRl = queueR(Enr,Enc,Enl);
               EnqLt = queueL(Ent,Enc,Enb);
               EnqRt = queueR(Entt,Ent,Enc);
               EnqLb = queueL(Enc,Enb,Enbb);
               EnqRb = queueR(Ent,Enc,Enb);

               pqLr = queueL(pr,pc,pl);
               pqRr = queueR(prr,pr,pc);
               pqLl = queueL(pc,pl,pll);
               pqRl = queueR(pr,pc,pl);
               pqLt = queueL(pt,pc,pb);
               pqRt = queueR(ptt,pt,pc);
               pqLb = queueL(pc,pb,pbb);
               pqRb = queueR(pt,pc,pb);


               xflux2D(fl,rhol,xmoml,ymoml,pl,Enl);
               xflux2D(fcx,rhoc,xmomc,ymomc,pc,Enc);
               xflux2D(fr,rhor,xmomr,ymomr,pr,Enr);
               yflux2D(ft,rhot,xmomt,ymomt,pt,Ent);
               yflux2D(fb,rhob,xmomb,ymomb,pb,Enb);
               yflux2D(fcy,rhoc,xmomc,ymomc,pc,Enc);

               flux2DLF(Fl,fcx,fl,rhol,rhoc,xmoml,xmomc,ymoml,ymomc,Enl,Enc,alpha);
               flux2DLF(Fr,fr,fcx,rhoc,rhor,xmomc,xmomr,ymomc,ymomr,Enc,Enr,alpha);
               flux2DLF(Fb,fcy,fb,rhob,rhoc,xmomb,xmomc,ymomb,ymomc,Enb,Enc,alpha);
               flux2DLF(Ft,ft,fcy,rhoc,rhot,xmomc,xmomt,ymomc,ymomt,Enc,Ent,alpha);

               xflux2D(fRr,rhoqRr,xmomqRr,ymomqRr,pqRr,EnqRr);
               xflux2D(fLr,rhoqLr,xmomqLr,ymomqLr,pqLr,EnqLr);
               xflux2D(fRl,rhoqRl,xmomqRl,ymomqRl,pqRl,EnqRl);
               xflux2D(fLl,rhoqLl,xmomqLl,ymomqLl,pqLl,EnqLl);

               yflux2D(fRt,rhoqRt,xmomqRt,ymomqRt,pqRt,EnqRt);
               yflux2D(fLt,rhoqLt,xmomqLt,ymomqLt,pqLt,EnqLt);
               yflux2D(fRb,rhoqRb,xmomqRb,ymomqRb,pqRb,EnqRb);
               yflux2D(fLb,rhoqLb,xmomqLb,ymomqLb,pqLb,EnqLb);

               flux2DLF(Fr,fRr,fLr,rhoqLr,rhoqRr,xmomqLr,xmomqRr,ymomqLr,ymomqRr,EnqLr,EnqRr,alpha);
               flux2DLF(Fl,fRl,fLl,rhoqLl,rhoqRl,xmomqLl,xmomqRl,ymomqLl,ymomqRl,EnqLl,EnqRl,alpha);

               flux2DLF(Ft,fRt,fLt,rhoqLt,rhoqRt,xmomqLt,xmomqRt,ymomqLt,ymomqRt,EnqLt,EnqRt,alpha);
               flux2DLF(Fb,fRb,fLb,rhoqLb,rhoqRb,xmomqLb,xmomqRb,ymomqLb,ymomqRb,EnqLb,EnqRb,alpha);

               for(k=0;k<4;k++) k_RK[k] = s_dt*(((Fl[k]-Fr[k])/(dx)+(Fb[k]-Ft[k])/(dy))) ;

               drhoRK[ind]  = k_RK[0];
               dxmomRK[ind] = k_RK[1];
               dymomRK[ind] = k_RK[2];
               dEnRK[ind]   = k_RK[3];
            }
         }

         for(i=0;i<(s_msh.numpoints[1]-1);i++)
         {
            for(n=0;n<(s_msh.numpoints[0]-1);n++)
            {
               ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;

               rhoRK[ind] = s_cont.rho[ind] + mfact[j]*drhoRK[ind];
               xmomRK[ind] =  s_cont.rho[ind]*s_cont.U[ind*s_msh.vecdims] + mfact[j]*dxmomRK[ind];
               ymomRK[ind] =  s_cont.rho[ind]*s_cont.U[ind*s_msh.vecdims+1] + mfact[j]*dymomRK[ind];
               EnRK[ind] = s_cont.En[ind] + mfact[j]*dEnRK[ind];
               pRK[ind] = (s_cont.gam-1.0)*(EnRK[ind] - 0.5*(xmomRK[ind]*xmomRK[ind]+ymomRK[ind]*ymomRK[ind])/rhoRK[ind]);

               rhonew[ind] =  rhonew[ind] + kfact[j]*drhoRK[ind]; 
               Unew[ind*s_msh.vecdims] =  Unew[ind*s_msh.vecdims] + kfact[j]*dxmomRK[ind]/rhonew[ind]; 
               Unew[ind*s_msh.vecdims+1] =  Unew[ind*s_msh.vecdims+1] + kfact[j]*dymomRK[ind]/rhonew[ind]; 
               Ennew[ind] =  Ennew[ind]+ kfact[j]*dEnRK[ind]; 
               pnew[ind] = (s_cont.gam-1.0)*(Ennew[ind] - 0.5*rhonew[ind]*(Unew[ind*s_msh.vecdims]*Unew[ind*s_msh.vecdims]+Unew[ind*s_msh.vecdims+1]*Unew[ind*s_msh.vecdims+1]));
            }
            //std::cout <<"\n" << rhonew[i] << "\t" << Unew[i] << "\t" << Ennew[i] << "\t" << pnew[i] << "\n"; 
         }
      }

      for(i=0;i<(s_msh.numpoints[1]-1);i++)
      {
         for(n=0;n<(s_msh.numpoints[0]-1);n++)
         {
            ind = (i+2)*(s_msh.numpoints[0]-1+2*s_msh.numghost) + n + 2;        
            s_cont.rho[ind] = rhonew[ind];
            s_cont.U[ind*s_msh.vecdims] = Unew[ind*s_msh.vecdims];
            s_cont.U[ind*s_msh.vecdims+1] = Unew[ind*s_msh.vecdims+1];
            s_cont.En[ind] = Ennew[ind];
            s_cont.p[ind] = pnew[ind];
         }
      }
   }

}


void solver::flux1D(std::vector<double> &s_f, double s_rho,double s_mom, double s_p, double s_En)
{
   s_f[0] = s_mom;
   s_f[1] = s_mom*s_mom/s_rho + s_p;
   s_f[2] = s_mom/s_rho*(s_En+s_p);
}


void solver::xflux2D(std::vector<double> &s_f, double s_rho,double s_xmom,double s_ymom, double s_p, double s_En)
{
   s_f[0] = s_xmom;
   s_f[1] = s_xmom*s_xmom/s_rho + s_p;
   s_f[2] = s_xmom*s_ymom/s_rho;
   s_f[3] = s_xmom/s_rho*(s_En+s_p);
}


void solver::yflux2D(std::vector<double> &s_f, double s_rho,double s_xmom,double s_ymom, double s_p, double s_En)
{
   s_f[0] = s_ymom;
   s_f[1] = s_xmom*s_ymom/s_rho;
   s_f[2] = s_ymom*s_ymom/s_rho + s_p;
   s_f[3] = s_ymom/s_rho*(s_En+s_p);
}


void solver::flux1DLF(std::vector<double> &s_F, const std::vector<double> &s_f2, const std::vector<double> &s_f1, double s_rho1, double s_rho2,double s_mom1, double s_mom2, double s_En1, double s_En2, double alpha)
{
   s_F[0] = (0.5)*(s_f1[0]+s_f2[0]) - (0.5)*alpha*(s_rho2-s_rho1);
   s_F[1] = (0.5)*(s_f1[1]+s_f2[1]) - (0.5)*alpha*(s_mom2-s_mom1);
   s_F[2] = (0.5)*(s_f1[2]+s_f2[2]) - (0.5)*alpha*(s_En2-s_En1);
}

void solver::flux2DLF(std::vector<double> &s_F, const std::vector<double> &s_f2, const std::vector<double> &s_f1, double s_rho1, double s_rho2, double s_xmom1, double s_xmom2,double s_ymom1, double s_ymom2, double s_En1, double s_En2, double alpha)
{
   s_F[0] = (0.5)*(s_f1[0]+s_f2[0]) - (0.5)*alpha*(s_rho2-s_rho1);
   s_F[1] = (0.5)*(s_f1[1]+s_f2[1]) - (0.5)*alpha*(s_xmom2-s_xmom1);
   s_F[2] = (0.5)*(s_f1[2]+s_f2[2]) - (0.5)*alpha*(s_ymom2-s_ymom1);
   s_F[3] = (0.5)*(s_f1[3]+s_f2[3]) - (0.5)*alpha*(s_En2-s_En1);
}



double solver::queueL(double up,double uc, double um)
{
   double qL;
   double theta;

   theta = thetaMUSCL(up,uc,um);
   qL = uc + 0.5*theta*(up-uc);

   return qL;
}

double solver::queueR(double upp,double up,double uc)
{
   double qR;
   double theta;

   theta = thetaMUSCL(upp,up,uc);
   qR = up - 0.5*theta*(upp-up);

   return qR;
}


double solver::thetaMUSCL(double up,double uc,double um)
{
   double ar = (uc-um)/(up-uc);
   double temp,theta;  
   double minval1=std::min(2.0*ar,1.0);
   double minval2=std::min(ar,2.0);

   temp = std::max(0.0,minval1);
   theta  = std::max(temp,minval2);

   return theta;
}

*/

#include "initializesolver.h"
#include "writeoutput.h"
#include "mesh.h"
#include "solver.h"
#include "solvertype.h"
#include "boundary.h"
#include "solverdefs.h"

//..Electrostatic 1D PIC solver..//

void sType::electrostatic1DPIC()
{
   int i,j,k,numthreads,tid;
   int numprocs,procid;


   std::vector<particles> prt;
   std::vector<contnm> cnt;
   MPIvars mpiv;

   //..Initialize code..//

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   srand(time(NULL)+procid);  //Seed Random Number Generator

   initializeSolver pr;
   solverVars svar;
   mesh msh;
   fields EM;
   solver slvr;
   boundary bnd;
   mathFunctions mth;
   std::vector<boundvars> bdv;
   std::vector<spclvars> spclv;

   pr.readdata(svar,msh,bdv,spclv);
   slvr.coulForce = 0;

   //omp_set_num_threads(svar.nsp);

   /*#pragma omp parallel private(tid)  // OMP
     {
     tid = omp_get_thread_num();
     if(tid==0)
     {
     numthreads = omp_get_num_threads();
     std::cout << "\n\nNumber of threads: \t" << numthreads << std::endl;
     }
     }*/

   writeOutput wrt(0,svar.nsp);

   msh.generatePICmesh(svar.lscale); 


   if(procid==0) //MPI
   {
      std::cout << "\nRunning with " << numprocs << " processors...\n"; 

      wrt.writepMesh(msh);
      wrt.writecMesh(msh);
   }

   svar.totalTime = 0.0;

   pr.initializePIC(prt,EM,cnt,msh,bdv,spclv,svar,wrt,mpiv);  //MPI

   slvr.deltaT = svar.dt; 
   slvr.totalTime = svar.totalTime;

   //for(i=0;i<svar.nsp;i++) wrt.writePICField(msh, cnt[i], EM,slvr.totalTime);
   //for(j=0;j<svar.nsp;j++) wrt.writeSingleParticle(prt[j],EM,msh,slvr.totalTime,msh.meshdims,msh.vecdims);
   //wrt.writePICField(msh, cnt, EM,slvr.totalTime);

   //for(j=0;j<svar.nsp;j++) slvr.findseedproc(prt[j]);
   bnd.seedParticles(msh,cnt,prt,svar,bdv);  
   if(svar.nspcl>0) bnd.applySpecialRegionsPIC(msh,cnt,prt,svar,spclv,EM);

   /*for(j=0;j<svar.nsp;j++) 
     {
     for(k=0;k<numprocs;k++) //MPI
     {
     MPI_Barrier(MPI_COMM_WORLD);
     if(k==procid) wrt.writeParticles(prt[j],EM,msh,slvr.totalTime);
     } //MPI
     }*/

#if SIMTYPE==2
   // only for the collision test simulation
   if(svar.outvdf==1)
   {
      for(j=0;j<prt[0].nsp;j++) wrt.findvdf(prt[j],msh,svar.outvdf_flag,slvr.totalTime);  //MPI
      for(j=0;j<prt[0].nsp;j++) wrt.findphasespace(prt[j],msh,svar.outvdf_flag,slvr.totalTime);   //MPI
   }
#endif

   auto t1 = std::chrono::high_resolution_clock::now();

   ////////////...............Main loop..........///////////////////////

   std::cout << "\n\tBeginning Main Loop\n";

   for(i=1; i<(svar.iter+1); i++)
   {
#if SIMTYPE==0 || SIMTYPE==3
      // Only for the full simulation or single particle simulation
      //bnd.applyContinuumBoundaryConditions(msh,EM,cnt,bdv);
      
      // Calculate electric potential (phi)
      slvr.poisson1D(msh,EM.phi,cnt,bdv,prt,slvr.deltaT,slvr.totalTime);   //..MB Comment out for

      // Calculate the electric field from the potential
      slvr.phitoE(EM.phi,EM.E,msh); //..MB Comment out for

      // set the electric field at the boundary
      bnd.applyEfieldBoundary(msh,EM,bdv); //..MB Comment out for
#endif

      if(svar.rst==0 || i>1) //..Avoid double update on restart
      {
         for(j=0;j<svar.nsp;j++) slvr.updatePartVel(prt[j], EM, msh, svar.mvscheme);
      }

      for(j=0;j<svar.nsp;j++) slvr.checkNAN(prt[j], EM, msh); //MPI

      if((i-1)%svar.p_iter == 0)  //........File Outputs
      {
         MPI_Barrier(MPI_COMM_WORLD);  //MPI

         for(j=0;j<svar.nsp;j++) slvr.redistributeparticles(prt[j],msh.vecdims,msh.meshdims);  //MPI
         MPI_Barrier(MPI_COMM_WORLD);  //MPI

         slvr.weighContinuumPIC(prt,cnt,msh,mpiv);  //CHGMPI
         for(j=0;j<svar.nsp;j++) mth.smoothData(cnt[j].N,msh.smthcnt);

         if(procid==0)
         {
            bnd.applyContinuumBoundaryConditions(msh,EM,cnt,bdv);
            if(svar.outcont==1) for(j=0;j<cnt[0].nsp;j++) wrt.writePICField(msh, cnt[j], EM,slvr.totalTime);
            if(svar.outcont==1) wrt.writePICField(msh, cnt, EM,slvr.totalTime);

         }
         MPI_Barrier(MPI_COMM_WORLD); //MPI
         if(svar.outpart==1) 
         {
            for(j=0;j<prt[0].nsp;j++) 
            {
               MPI_Barrier(MPI_COMM_WORLD);
               wrt.writeParticles(prt[j],EM,msh,slvr.totalTime);
            }
         }
         if(svar.outvdf==1)  
         {
            for(j=0;j<prt[0].nsp;j++) wrt.findvdf(prt[j],msh,svar.outvdf_flag,slvr.totalTime);//MPI
            for(j=0;j<prt[0].nsp;j++) wrt.findphasespace(prt[j],msh,svar.outvdf_flag,slvr.totalTime);   //CMPI
         }
      }
      else if(i%svar.outfinalskip == 0 && i > (svar.iter-svar.outfinal))
      {
         MPI_Barrier(MPI_COMM_WORLD);  //MPI

         slvr.weighContinuumPIC(prt,cnt,msh,mpiv);
         for(j=0;j<svar.nsp;j++) mth.smoothData(cnt[j].N,msh.smthcnt);

         if(procid==0)
         {
            bnd.applyContinuumBoundaryConditions(msh,EM,cnt,bdv);
            if(svar.outcont==1) for(j=0;j<cnt[0].nsp;j++) wrt.writePICField(msh, cnt[j], EM,slvr.totalTime);
            if(svar.outcont==1) wrt.writePICField(msh, cnt, EM,slvr.totalTime);

            //for(j=0;j<prt[0].nsp;j++) wrt.writeSingleParticle(prt[j],EM,msh,slvr.totalTime);
         }
         MPI_Barrier(MPI_COMM_WORLD); //MPI
         if(svar.outpart==1) 
         {
            for(j=0;j<prt[0].nsp;j++) 
            {
               MPI_Barrier(MPI_COMM_WORLD);
               wrt.writeParticles(prt[j],EM,msh,slvr.totalTime);
            }
         }
         if(svar.outvdf==1)  
         {
            for(j=0;j<prt[0].nsp;j++) wrt.findvdf(prt[j],msh,svar.outvdf_flag,slvr.totalTime);//MPI
            for(j=0;j<prt[0].nsp;j++) wrt.findphasespace(prt[j],msh,svar.outvdf_flag,slvr.totalTime);   //MPI
         }
      }

      if((i-1)%svar.p_iter == 0 && procid==0 && svar.outrst==1)  wrt.writeRestart(bdv,slvr.totalTime,i,svar.nsp);  
      else if(i%svar.outfinalskip == 0 && i > (svar.iter-svar.outfinal) && svar.outrst==1) wrt.writeRestart(bdv,slvr.totalTime,i,svar.nsp);  

#if SIMTYPE==3
      for(j=0;j<prt[0].nsp;j++) wrt.writeSingleParticle(prt[j],EM,msh,slvr.totalTime);
#endif

      if(svar.outavg==1)  wrt.findglobalenergy(prt,EM,msh,slvr.totalTime);  //MPI

      for(j=0;j<svar.nsp;j++)  slvr.updatePartPos(prt[j], msh, svar.mvscheme);
      bnd.applyParticleBoundaryConditions(msh,cnt,prt,svar,bdv);
      for(j=0;j<svar.nsp;j++) msh.connectpartandmesh(prt[j]);

      //slvr.weighDensityPIC(prt,cnt,msh,mpiv);  //CHGMPI
      //for(j=0;j<svar.nsp;j++) mth.smoothData(cnt[j].N,msh.smthcnt);

      MPI_Barrier(MPI_COMM_WORLD);

#if SIMTYPE!=1
      if(svar.nct>0)
      {   
         for(j=0;j<svar.nsp;j++) slvr.collideParticles(prt,msh,cnt[svar.nsp],slvr.totalTime,j);  //MB
         //slvr.updateNeutralBackground(cnt[svar.nsp],msh,slvr.deltaT);
      }
#endif


      slvr.weighContinuumPIC(prt,cnt,msh,mpiv);   //FIX SUBCYCLE
      for(j=0;j<svar.nsp;j++) mth.smoothData(cnt[j].N,msh.smthcnt);
      if(svar.coul>0)  slvr.coulombCollisions(prt,msh,cnt,mpiv); 
      MPI_Barrier(MPI_COMM_WORLD);
      if(svar.nspcl>0)  bnd.applySpecialRegionsPIC(msh,cnt,prt,svar,spclv,EM);

      for(j=0;j<svar.nsp;j++) msh.connectpartandmesh(prt[j]);

      //slvr.weighContinuumPIC(prt,cnt,msh,mpiv);   //FIX SUBCYCLE
      slvr.weighDensityPIC(prt,cnt,msh,mpiv);   //FIX SUBCYCLE
      for(j=0;j<svar.nsp;j++) mth.smoothData(cnt[j].N,msh.smthcnt);

      //bnd.applyContinuumBoundaryConditions(msh,EM,cnt,bdv);  //MPI

      slvr.totalTime = slvr.totalTime + slvr.deltaT;
      svar.totalTime = slvr.totalTime;
      if(procid==0)  std::cout << "\nIT:  " << i << "\t" << "Time: " << slvr.totalTime;


   }

   auto t2 = std::chrono::high_resolution_clock::now();
   std::cout << "\n\nAnnnnnnd Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() << " ms\n\n";

   //slvr.updateAppliedPotential(pr.init_phi,EM.phi,msh);
   //slvr.timestep(ions,msh.pmesh, msh.meshlength, msh.numpoints, svar.dt, svar.cfl,msh.meshdims,msh.vecdims);

}

/*

//..2D electromagnetic  PIC..//
//Not implemented

void sType::electromagnetic2DPIC()
{
int i;
std::vector<int> test,test2;
for(i=0;i<4;i++)  test.push_back(0);   
for(i=0;i<2;i++)  test2.push_back(0.0);   


std::cout << "Howdy, Welcome to Frans' PIC Solver \n \n";

initializeSolver pr;
solverVars svar;
mesh msh;
particles ions;
fields    EM;
flow      cont;
solver slvr;
boundary bnd;
writeOutput wrt(0,1);

}

//..1D Euler..//
//Not implemented

void sType::euler1D()
{
int i;

std::vector<int> test,test2;
for(i=0;i<4;i++)  test.push_back(0);   
for(i=0;i<2;i++)  test2.push_back(0.0);   

initializeSolver pr;
solverVars svar;
mesh msh;
flow  fluid;
solver slvr;
writeOutput wrt(0,1);
boundary bnd;
std::vector<boundvars> bdv;
std::vector<spclvars> spclv;

pr.readdata(svar,msh,bdv,spclv);

msh.generateEulermesh(svar.lscale); 
wrt.writepMesh(msh);
wrt.writecMesh(msh);

pr.initializeEuler(fluid,msh,bnd,svar);

slvr.deltaT = svar.dt; 

std::clock_t start;
start = std::clock();

wrt.writeFlowEuler(msh,fluid);

for(i=1; i<(svar.iter+1); i++)
{
bnd.applyPeriodicBoundaryConditions(msh,fluid);
slvr.updateEulerFlow(fluid,msh,slvr.deltaT,svar.tdflag,svar.sdflag);
if(i%svar.p_iter == 0) 
{
std::cout << "IT:" << i << std::endl;
wrt.writeFlowEuler(msh,fluid);
}
}

std::cout << "\n\nVolume Average Error:\t" << std::setprecision(10) << slvr.errorCheck(fluid,msh) ;

std::cout << "\n\nTotal Time:" << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << " ms" << std::endl ;

}


//..2D Euler..//
//Not implemented

void sType::euler2D()
{
   int i;

   std::vector<int> test,test2;
   for(i=0;i<4;i++)  test.push_back(0);   
   for(i=0;i<2;i++)  test2.push_back(0.0);   

   initializeSolver pr;
   solverVars svar;
   mesh msh;
   flow  fluid;
   solver slvr;
   writeOutput wrt(0,1);
   boundary bnd;
   std::vector<boundvars> bdv;
   std::vector<spclvars> spclv;

   pr.readdata(svar,msh,bdv,spclv);


   msh.generateEulermesh(svar.lscale); 
   wrt.writepMesh(msh);
   wrt.writecMesh(msh);

   pr.initializeEuler(fluid,msh,bnd,svar);

   slvr.deltaT = svar.dt; 

   std::clock_t start;
   start = std::clock();

   wrt.writeFlowEuler(msh,fluid);

   for(i=1;i<(svar.iter+1); i++)
   {
      bnd.applyBoundaryConditionsEuler(msh,fluid);
      slvr.updateEulerFlow2D(fluid,msh,slvr.deltaT,svar.tdflag,svar.sdflag);

      std::cout << "IT:" << i << std::endl;
      if(i%svar.p_iter == 0) 
      {
         wrt.writeFlowEuler(msh,fluid);
      }
   }

   std::cout << "Total Time:\t" << (std::clock() - start) / (double)(CLOCKS_PER_SEC / 1000) << "\t ms" << std::endl ;

}

*/

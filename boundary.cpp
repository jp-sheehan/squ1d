#include "boundary.h"

//..Initialize boundary functions..//

void boundary::initializeBoundary(mesh &b_msh,std::vector<boundvars> &b_bdv, std::vector<particles> &b_part)
{
   int i,j,k;
   int bcell;

   findBoundaryCells(b_msh,b_bdv);  //Find boundary cells
   initializeParticleCounter(b_bdv,b_part); //Initialize particle counter for leaving domain

   for(j=0;j<b_part[0].nsp;j++) cleanParticles(b_part[j],b_msh.meshbd,b_bdv,b_msh,b_msh.vecdims,b_msh.meshdims,b_msh.perflag); 

   for(i=0;i<b_bdv[0].nbound;i++)
   {
      if(b_msh.meshdims==1)
      {
         bcell=b_bdv[i].cboundcells[0];  
         if(bcell< b_msh.incneighofc(bcell))  b_bdv[i].fdir = 1;
         else if(bcell>b_msh.incneighofc(bcell))  b_bdv[i].fdir  = -1;
      }
   }
}

//..Initialize Special regions..//

void boundary::initializeSpecialRegions(mesh &b_msh,std::vector<spclvars> &b_spclv, std::vector<particles> &b_part)
{
   int i,j,k;

   //std::cout << "\n\tApplying Special Regions....";

   findSpecialRegionCells(b_msh,b_spclv);  //Find special region cells

   //std::cout << "x";
}

//..Apply particle boundary conditions..//

void boundary::applyParticleBoundaryConditions(mesh b_msh, std::vector<contnm> &b_cont, std::vector<particles> &b_part, solverVars b_svar, std::vector<boundvars> &b_bdv)
{
   int i,j,k;

   for(j=0;j<b_svar.nsp;j++) cleanParticles(b_part[j],b_msh.meshbd,b_bdv,b_msh,b_msh.vecdims,b_msh.meshdims,b_msh.perflag);

   if(b_bdv[0].clss=="SOURCE" || b_bdv[1].clss=="SOURCE")   if(b_msh.perflag==0)  seedParticles(b_msh,b_cont,b_part,b_svar,b_bdv);

}

//..Initialize particle counter for those leaving domain..//

void boundary::initializeParticleCounter(std::vector<boundvars> &b_bdv, std::vector<particles> b_part)
{
   int i,j,k;

   for(i=0;i<b_bdv.size();i++)
   {
      //if(b_bdv[i].clss=="WALL" && b_bdv[i].wallfloat==1)
      //{
      for(j=0;j<b_part[0].nsp;j++) b_bdv[i].partcount.push_back(0);
      //}
      b_bdv[i].sigma = 0.0;
      b_bdv[i].Ebd = 0.0;
   }   

}


void boundary::setParticleCounter(std::vector<boundvars> &b_bdv, std::vector<particles> b_part, int i_rst)
{
   int i,j,k;
   int i_iter,temp;
   double tempdouble;
   std::stringstream fname;
   std::string filename;

   std::string name;

   std::cout << "\n\tReading Restart File for particle counter...." ;

   name = "restart";

   fname <<  name.c_str() << i_rst << ".out";

   filename = fname.str();

   std::ifstream rdfile(filename.c_str());

   rdfile >> tempdouble >> tempdouble;

   for(i=0;i<b_part[0].nsp;i++)
   {
      rdfile >> temp;
      rdfile >> temp;
      rdfile >> temp;
      rdfile >> temp;
      rdfile >> temp;
   }

   rdfile >> temp >> temp;

   for(i=0;i<b_bdv[0].nbound;i++)
   {
      //if(b_bdv[i].clss=="WALL" && b_bdv[i].wallfloat==1)
      //{
      for(j=0;j<b_part[0].nsp;j++)
      {
         rdfile >> temp;
         b_bdv[i].partcount[j] = temp;
      }
      rdfile >> tempdouble;
      b_bdv[i].sigma = tempdouble;
      b_bdv[i].Ebd = 0.0;
      //}
   }   

   rdfile.close();

}


//..Seed particles at boundary..//

void boundary::seedParticles(mesh b_msh, std::vector<contnm> &b_cont, std::vector<particles> &b_part, solverVars b_svar,std::vector<boundvars> &b_bdv)
{
   int i,j,k,l,m,nbcs;
   int bcell;
   int numprocs,procid;  //MPI
   double dens,temp;
   std::vector<double> vel;
   std::vector<double> meshpt;
   std::vector<double> meshptin;

   for(i=0;i<b_msh.meshdims;i++) meshpt.push_back(0.0);
   for(i=0;i<b_msh.meshdims;i++) meshptin.push_back(0.0);
   for(i=0;i<b_msh.vecdims;i++) vel.push_back(0.0);

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   solver b_slvr;

   double area = 1.0;   
   nsp = b_bdv[0].nsp;
   nbound = b_bdv[0].nbound;
   nbcs = 3*nsp+2;

   //std::cout << "\n\tSeeding Particles......";

   for(m=0;m<nsp;m++) b_part[m].gp=0;

   for(i=0;i<nbound;i++)
   {
      for(m=0;m<nsp;m++)
      { 

         //..............Density....................//

         if(b_bdv[i].boundtype[3*m] == "DIRICHLET")
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
               dens = beqnparser(b_bdv[i].bddens[m],meshpt);  //UPDATE
            }
         }
         else if(b_bdv[i].boundtype[3*m] == "NEUMANN")  //CHG
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               dens = b_cont[m].N[b_msh.incneighofc(bcell)];
            }
         }
         else if(b_bdv[i].boundtype[3*m] == "PERIODIC")  dens=0.0;  //CHG
         else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";

         //.............Velocity....................//

         if(b_bdv[i].boundtype[3*m+1] == "DIRICHLET")
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
               for(k=0;k<b_msh.vecdims;k++)  vel[k] = beqnparser(b_bdv[i].bdvel[k],meshpt);  //UPDATE
            }
         }
         else if(b_bdv[i].boundtype[3*m+1] == "NEUMANN")  //CHG
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.vecdims;k++)  vel[k] = b_cont[m].U[b_msh.vecdims*b_msh.incneighofc(bcell)+k];
            }
         }
         else if(b_bdv[i].boundtype[3*m+1] == "PERIODIC") 
         {}   //CHG
         else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";

         //..............Temperature....................//

         if(b_bdv[i].boundtype[3*m+2] == "DIRICHLET")
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
               temp = beqnparser(b_bdv[i].bdtemp[m],meshpt); //UPDATE
            }
         }
         else if(b_bdv[i].boundtype[3*m+2] == "NEUMANN")   //CHG
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               temp = b_cont[m].T[b_msh.incneighofc(bcell)];
            }
         }
         else if(b_bdv[i].boundtype[3*m+2] == "PERIODIC") temp = 0.0; //CHG 
         else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";

         //............Call Seeding algorithm.................//

         if(b_bdv[i].clss=="SOURCE")
         {
            if(b_bdv[i].boundtype[3*m] != b_bdv[i].boundtype[3*m+1] || b_bdv[i].boundtype[3*m]  != b_bdv[i].boundtype[3*m+2])
            {
               std::cout <<  "....ERROR:  BC's don't match for particle source......";
               exit(EXIT_FAILURE);
            }

            if(b_msh.meshdims==1)
            {
               //if(procid == b_part[m].seedproc) // MPI
               //{    
               bcell=b_bdv[i].cboundcells[0];  
               if(b_bdv[i].srctype == "FLUX") b_slvr.particlefluxsource1D(b_part[m], b_msh, b_bdv[i].bdddist[m], b_bdv[i].bdthdist[m], vel, dens, temp, b_svar.dt, bcell,b_bdv[i].fdir,b_bdv[i].partcount[m]);
               else if(b_bdv[i].srctype == "CELL") b_slvr.particlecellsource(b_part[m], b_msh, b_bdv[i].bdddist[m], b_bdv[i].bdthdist[m], vel, dens, temp, bcell); 
               //}

               for(j=0;j<b_cont[0].nsp;j++) b_bdv[i].partcount[j] = 0;  //Re-seeding due to particle reflection 
            }
         }
      }
   } 
   //std::cout << "x";
}

//..Clear boundary values..//

void boundary::clearallbdv(boundvars &b_bdv)
{
   b_bdv.boundrange.clear();
   b_bdv.boundtype.clear();
   b_bdv.bdvel.clear();
   b_bdv.bddens.clear();
   b_bdv.bdddist.clear();
   b_bdv.bdtemp.clear();
   b_bdv.bdthdist.clear();
   b_bdv.bdE.clear();
   b_bdv.bdB.clear();
   b_bdv.bdU.clear(); 
}

//..Remove, move, and count particles that are outside the domain..//

void boundary::cleanParticles(particles &b_part, std::vector<double> b_meshbd, std::vector<boundvars> &b_bdv, mesh &b_msh, int b_vec_dims, int b_mesh_dims, int b_perflag)
{
   int i,j;
   int ind,count,num_particles;
   int numprocs, procid;  //MPI
   int total;  //MPI
   std::vector<double> area;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   num_particles = b_part.pos.size()/b_mesh_dims;
   count = 0;
   int lcount = 0;
   int rcount = 0;

   if(b_msh.meshdims==1) //..Interpolating areas for surface charge located at cell interface 
   {
      if(b_msh.philoc == 1)
      {
         area.push_back(b_msh.parea[b_bdv[0].pboundcells[0]]);
         area.push_back(b_msh.parea[b_bdv[1].pboundcells[0]]);
      }
      else
      {
         area.push_back((b_msh.carea[b_bdv[0].cboundcells[0]]+b_msh.carea[b_msh.cintcells[0]])*0.5);
         area.push_back((b_msh.carea[b_bdv[1].cboundcells[0]]+b_msh.carea[b_msh.cintcells[b_msh.cintcells.size()-1]])*0.5);
      }
   }

   //std::cout << "\n\tCleaning and collecting particles...";

   if(b_perflag==0)
   {
      for(i=0;i<num_particles;i++)
      {
         ind = i+count;

         for(j=0;j<b_mesh_dims;j++)   //..1D
         {
            if(b_part.pos[b_mesh_dims*ind+j] > b_meshbd[2*j+1])  
            {
               b_part.pos.erase(b_part.pos.begin() + b_mesh_dims*ind, b_part.pos.begin() + b_mesh_dims*ind+b_mesh_dims);
               b_part.vel.erase(b_part.vel.begin() + b_vec_dims*ind, b_part.vel.begin() + b_vec_dims*ind+b_vec_dims);
               b_part.en.erase(b_part.en.begin() + ind, b_part.en.begin() + ind + 1);
               b_part.cell.erase(b_part.cell.begin() + ind, b_part.cell.begin() + ind + 1);
               b_part.pt.erase(b_part.pt.begin() + ind, b_part.pt.begin() + ind + 1);
               count = count-1;
               rcount += 1;
            }
            else if(b_part.pos[b_mesh_dims*ind] < b_meshbd[2*j])  
            {
               //std::cout << std::endl << ind << "\t" << count << "\t" << b_meshbd[2*j] << "\t"  << b_part.pos[b_mesh_dims*ind] << std::endl;
               b_part.pos.erase(b_part.pos.begin() + b_mesh_dims*ind, b_part.pos.begin() + b_mesh_dims*ind+b_mesh_dims);
               b_part.vel.erase(b_part.vel.begin() + b_vec_dims*ind, b_part.vel.begin() + b_vec_dims*ind+b_vec_dims);
               b_part.en.erase(b_part.en.begin() + ind, b_part.en.begin() + ind + 1);
               b_part.cell.erase(b_part.cell.begin() + ind, b_part.cell.begin() + ind + 1);
               b_part.pt.erase(b_part.pt.begin() + ind, b_part.pt.begin() + ind + 1);
               count = count-1;
               lcount += 1;
            }
         }
      }

      if(numprocs>1)
      {
         total = 0.0;  //MPI
         MPI_Allreduce(&lcount,&total,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD); //MPI
         lcount = total; //MPI

         total = 0.0; //MPI
         MPI_Allreduce(&rcount,&total,1, MPI_INT,MPI_SUM,MPI_COMM_WORLD); //MPI
         rcount = total; //MPI
      }

      //std::cout << "l: " << lcount << "r: " << rcount << "\t" << total << std::endl;

      for(i=0;i<b_bdv[0].nbound;i++)
      {
         //if(b_bdv[i].clss == "WALL" && (b_bdv[i].wallfloat==1 || b_bdv[i].wallcur==1 || b_bdv[i].wallcap==1)) 
         //{
         if(b_bdv[i].fdir == 1) 
         {  
            b_bdv[i].partcount[b_part.pid] += lcount;
            b_bdv[i].sigma += lcount*b_part.wcharge/area[0];
         }
         else if(b_bdv[i].fdir == -1) 
         {
            b_bdv[i].partcount[b_part.pid] += rcount;
            b_bdv[i].sigma += rcount*b_part.wcharge/area[1];
         }
         //std::cout << "At " << b_bdv[i].boundname << " " << b_bdv[i].partcount[b_part.pid] << " " << b_part.name << " particles " << std::endl;
         //}
         /*else    // Sigma set elsewhere
           {
           if(b_bdv[i].fdir == 1) 
           {  
           b_bdv[i].partcount[b_part.pid] += lcount;
           }
           else if(b_bdv[i].fdir == -1) 
           {
           b_bdv[i].partcount[b_part.pid] += rcount;
           }
         //std::cout << "At " << b_bdv[i].boundname << " " << b_bdv[i].partcount[b_part.pid] << " " << b_part.name << " particles " << std::endl;
         }*/
      }
   }
   else if(b_perflag==1)
   {
      for(i=0;i<num_particles;i++)
      {
         for(j=0;j<b_mesh_dims;j++)
         {
            if(b_part.pos[b_mesh_dims*i+j] > b_meshbd[j*2+1])
            {
               b_part.pos[b_mesh_dims*i+j] = b_part.pos[b_mesh_dims*i+j] -(b_meshbd[2*j+1]-b_meshbd[2*j]);
            }
            else if(b_part.pos[b_mesh_dims*i+j] < b_meshbd[j*2])  
            {
               b_part.pos[b_mesh_dims*i+j] = b_part.pos[b_mesh_dims*i+j] + (b_meshbd[2*j+1]-b_meshbd[2*j]);
            }
         }
      }
   }

   //std::cout << "x";
}


//..Find cells which are on boundary..//

void boundary::findBoundaryCells(mesh &b_msh,std::vector<boundvars> &b_bdv)
{
   bool b_chck;
   int i,j,k;

   b_msh.ncbound = 0;
   b_msh.npbound = 0;

   //....Finds cells that are on or outside of the computational domain....//

   if(b_msh.meshdims==1)
   {
      for(i=0;i<b_msh.nbound;i++) b_bdv[i].cboundnum = b_msh.numghost;   //Number of ghost cells per boundary
      for(i=0;i<b_msh.nbound;i++) b_bdv[i].pboundnum = b_msh.numghost;

      for(i=0;i<b_msh.nbound;i++) b_msh.ncbound += b_bdv[i].cboundnum;    //Total boundary cells
      for(i=0;i<b_msh.nbound;i++) b_msh.npbound += b_bdv[i].cboundnum;    //Total boundary points

      std::cout << "\nBoundary Cells:\t" << b_msh.ncbound << "\t" << b_msh.npbound << "\n";

      for(i=0;i<b_msh.numghost;i++) b_bdv[0].cboundcells.push_back(i); //Boundary cell indeces
      for(i=0;i<b_msh.numghost;i++) b_bdv[1].cboundcells.push_back(b_msh.cmesh.size()-b_msh.numghost);
      for(i=b_msh.numghost;i<(b_msh.cmesh.size()-b_msh.numghost);i++) b_msh.cintcells.push_back(i); //Internal cell indices 

      b_bdv[0].pboundcells.push_back(0);
      b_bdv[1].pboundcells.push_back(b_msh.numpoints[0]-1);

      for(i=1;i<(b_msh.pmesh.size()-1);i++) b_msh.pintcells.push_back(i); 

   }
   else if(b_msh.meshdims==2)
   {
      /*  Check all of this.....

          for(i=0;i<b_msh.meshdims;i++) b_msh.cboundnum.push_back(b_msh.numpoints[0]-1+2*b_msh.numghost);
          for(i=0;i<b_msh.meshdims;i++) b_msh.cboundnum.push_back(b_msh.numpoints[0]-1);

          for(i=0;i<b_msh.meshdims;i++) b_msh.pboundnum.push_back(b_msh.numpoints[0]);
          for(i=0;i<b_msh.meshdims;i++) b_msh.pboundnum.push_back(b_msh.numpoints[0]-2);

          b_msh.ncbound = std::accumulate(b_msh.cboundnum.begin(),b_msh.cboundnum.end(),0);
          b_msh.npbound = std::accumulate(b_msh.pboundnum.begin(),b_msh.pboundnum.end(),0);

          std::cout << "\nBoundary Cells:\t" << b_msh.ncbound << "\t" << b_msh.npbound << "\n";

          for(i=0;i<(b_msh.numpoints[0]-1+2*b_msh.numghost);i++) b_msh.cboundcells.push_back(i);
          for(i=0;i<(b_msh.numpoints[0]-1+2*b_msh.numghost);i++) b_msh.cboundcells.push_back(i+b_msh.numpoints[1]*(b_msh.numpoints[0]-1+2*b_msh.numghost));
          for(i=0;i<(b_msh.numpoints[0]-1);i++) b_msh.cboundcells.push_back((i+1)*(b_msh.numpoints[0]-1+2*b_msh.numghost));
          for(i=0;i<(b_msh.numpoints[0]-1);i++) b_msh.cboundcells.push_back((i+2)*(b_msh.numpoints[0]-1+2*b_msh.numghost)-1);

          for(i=0;i<(b_msh.cmesh.size());i++)
          {
          b_chck = false; 
          for(j=0;j<b_msh.cboundcells.size();j++)
          {
          if(i==b_msh.cboundcells[j]) b_chck  = true;  
          }
          if(b_chck==false) b_msh.cintcells.push_back(i);
          }

          for(i=0;i<(b_msh.numpoints[0]);i++) b_msh.pboundcells.push_back(i);
          for(i=0;i<(b_msh.numpoints[0]);i++) b_msh.pboundcells.push_back(i+b_msh.numpoints[0]*(b_msh.numpoints[1]-1));
          for(i=0;i<(b_msh.numpoints[0]-2);i++) b_msh.pboundcells.push_back((i+1)*b_msh.numpoints[0]);
          for(i=0;i<(b_msh.numpoints[0]-2);i++) b_msh.pboundcells.push_back((i+2)*(b_msh.numpoints[0])-1);

          for(i=0;i<(b_msh.pmesh.size());i++)
          {
          b_chck = false ;
          for(j=0;j<b_msh.pboundcells.size();j++)
          {
          if(i==b_msh.pboundcells[j]) b_chck  = true;  
          }
          if(b_chck==false) b_msh.pintcells.push_back(i);
          }*/

   }

   /*for(i=0;i<b_msh.cboundnum.size();i++)  std::cout << b_msh.cboundnum[i] << std::endl;
     std::cout << std::endl << "Cell Boundary"  <<  std::endl;
     for(i=0;i<b_msh.cboundcells.size();i++)  std::cout << b_msh.cboundcells[i] << "\t";
     std::cout << std::endl << "Point Boundary"  <<  std::endl;
     for(i=0;i<b_msh.pboundcells.size();i++)  std::cout << b_msh.pboundcells[i] << "\t";
     std::cout << std::endl << "Cell Internal"  <<  std::endl;
     for(i=0;i<b_msh.cintcells.size();i++)  std::cout << b_msh.cintcells[i] << "\t";
     std::cout << std::endl << "Point Internal"  <<  std::endl;
     for(i=0;i<b_msh.pintcells.size();i++)  std::cout << b_msh.pintcells[i] << "\t";*/

}


//..Equation parser for boundary conditions..//

double boundary::beqnparser(std::string expression_string, std::vector<double> i_point) 
{
   int i_size = i_point.size();
   double i_x, i_y, i_z;
   exprtk::symbol_table<double> symbol_table;

   switch (i_size)
   {
      case 1:
         i_x = i_point[0];
         symbol_table.add_variable("x",i_x);
         break;
      case 2:
         i_x = i_point[0];
         i_y = i_point[1];
         symbol_table.add_variable("x",i_x);
         symbol_table.add_variable("y",i_y);
         symbol_table.add_variable("z",i_z);
         break;
      case 3:
         i_x = i_point[0];
         i_y = i_point[1];
         i_z = i_point[2];
         symbol_table.add_variable("x",i_x);
         symbol_table.add_variable("y",i_y);
         symbol_table.add_variable("z",i_z);
         break;
   }



   exprtk::expression<double> expression;
   expression.register_symbol_table(symbol_table);

   exprtk::parser<double> parser;
   parser.compile(expression_string,expression);

   return expression.value();
}

//..Apply continuum boundary conditions..//

void boundary::applyContinuumBoundaryConditions(mesh b_msh, fields &b_flds, std::vector<contnm> &b_cont,std::vector<boundvars> &b_bdv)
{
   int i,j,k,l,m,nbcs;
   int bcell;

   //CHG NOTES: Add U and T, Change B

   std::vector<double> meshpt;
   std::vector<double> meshptin;

   for(i=0;i<b_msh.meshdims;i++) meshpt.push_back(0.0);
   for(i=0;i<b_msh.meshdims;i++) meshptin.push_back(0.0);

   double area = 1.0;   

   nsp = b_bdv[0].nsp;
   nbound = b_bdv[0].nbound;
   nbcs = 3*nsp+2;

   //std::cout << "\n\tApplying Continuum BC's...";

   if(b_msh.philoc==1)
   {
      for(i=0;i<nbound;i++)
      {
         //..............Density and Energy....................//
         for(m=0;m<nsp;m++)
         { 
            if(b_bdv[i].boundtype[3*m] == "DIRICHLET")
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];
                  for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.pmesh[b_msh.meshdims*bcell+k];
                  //b_cont[m].N[bcell] = beqnparser(b_bdv[i].bddens[m],meshpt); //UPDATE
                  //std::cout << "\nBcell,mshpt "<< bcell << "\t" <<meshpt[0] << "\t";
                  b_bdv[i].bddens_val[m] = beqnparser(b_bdv[i].bddens[m],meshpt); //UPDATE
                  //std::cout << b_bdv[i].bddens_val[m];
                  //b_cont[m].En[bcell] = beqnparser(b_bdv[i].bdEn[m],meshpt);
               }
            }
            else if(b_bdv[i].boundtype[3*m] == "NEUMANN")
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];
                  if(b_msh.meshdims==1) 
                  {
                     meshpt[0] = b_msh.pmesh[bcell];
                     //if(meshpt[0]<b_msh.meshbd[0]) b_cont[m].N[bcell] =  b_cont[m].N[b_msh.incneighofc(bcell)] - (b_msh.meshbd[0]-meshpt[0])*beqnparser(b_bdv[i].bddens[m],meshpt); //UPDATE
                     //if(meshpt[0]<b_msh.meshbd[0]) b_cont[m].En[bcell] =  b_cont[m].En[b_msh.incneighofc(bcell)] - (b_msh.meshbd[0]-meshpt[0])*beqnparser(b_bdv[i].bdEn[m],meshpt);
                     //else if(meshpt[0]>b_msh.meshbd[1]) b_cont[m].N[bcell] =  b_cont[m].N[b_msh.incneighofc(bcell)] + (meshpt[0]-b_msh.meshbd[1])*beqnparser(b_bdv[i].bddens[m],meshpt);  //UPDATE
                     //b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                     //else if(meshpt[0]>b_msh.meshbd[1]) b_cont[m].En[bcell] =  b_cont[m].En[b_msh.incneighofc(bcell)] + (meshpt[0]-b_msh.meshbd[1])*beqnparser(b_bdv[i].bdEn[m],meshpt);
                  }
               }
            }
            else if(b_bdv[i].boundtype[3*m] == "PERIODIC") 
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];

                  if(b_msh.meshdims==1) 
                  {
                     //if(b_msh.intscheme==0)
                     //{
                     meshpt[0] = b_msh.pmesh[bcell];
                     if(meshpt[0]<=b_msh.meshbd[0])
                     {  
                        b_cont[m].N[bcell] =  b_cont[m].N[b_bdv[1].pboundcells[0]];
                        b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                        b_cont[m].En[bcell] =  b_cont[m].En[b_bdv[1].pboundcells[0]];
                     }
                     else if(meshpt[0]>=b_msh.meshbd[1])
                     {   
                        b_cont[m].N[bcell] =  b_cont[m].N[b_bdv[0].pboundcells[0]];
                        b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                        b_cont[m].En[bcell] =  b_cont[m].En[b_bdv[0].pboundcells[0]];
                     }
                     //}
                     /*else if(b_msh.intscheme==1)
                       {
                       meshpt[0] = b_msh.cmesh[bcell];
                       if(meshpt[0]<b_msh.meshbd[0])
                       {  
                       b_cont[m].N[bcell] =  b_cont[m].N[bcell] + b_cont[m].N[b_bdv[1].cboundcells[0]-b_msh.numghost];
                       b_cont[m].En[bcell] =  b_cont[m].En[bcell] + b_cont[m].En[b_bdv[1].cboundcells[0]-b_msh.numghost];
                       b_cont[m].N[b_bdv[1].cboundcells[0]-b_msh.numghost] = b_cont[m].N[bcell];
                       b_cont[m].En[b_bdv[1].cboundcells[0]-b_msh.numghost] = b_cont[m].En[bcell];
                       }
                       else if(meshpt[0]>b_msh.meshbd[1])
                       {   
                       b_cont[m].N[bcell] =  b_cont[m].N[bcell] + b_cont[m].N[b_bdv[0].cboundcells[0]+b_msh.numghost];
                       b_cont[m].En[bcell] =  b_cont[m].En[bcell] + b_cont[m].En[b_bdv[0].cboundcells[0]+b_msh.numghost];
                       b_cont[m].N[b_bdv[0].cboundcells[0]+b_msh.numghost] = b_cont[m].N[bcell];
                       b_cont[m].En[b_bdv[0].cboundcells[0]+b_msh.numghost] = b_cont[m].En[bcell];
                       }
                       }*/
                  }
               }
            }
            else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
         }

         //CHG:  Add Other Continuum Variables (U,T)

         if(b_msh.sflag==2)
         {
            //..............E-field....................//

            if(b_bdv[i].boundtype[3*nsp] == "DIRICHLET")
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];
                  for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.pmesh[b_msh.meshdims*bcell+k];
                  for(k=0;k<b_msh.vecdims;k++) b_flds.E[b_msh.vecdims*bcell+k] = beqnparser(b_bdv[i].bdE[k],meshpt); //UPDATE
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "NEUMANN") //CHG
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];
                  //for(k=0;k<b_msh.vecdims;k++)  b_flds.E[b_msh.vecdims*bcell+k] = b_flds.E[b_msh.vecdims*b_msh.inpneighofp(bcell)+k];
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "PERIODIC")
            {
            }
            else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
         }
         else if(b_msh.sflag==0)
         {
            //..............Potential....................//

            if(b_bdv[i].boundtype[3*nsp] == "DIRICHLET")
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];
                  for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.pmesh[b_msh.meshdims*bcell+k];
                  b_bdv[i].bdphi_val = beqnparser(b_bdv[i].bdphi,meshpt);  //UPDATE
                  if(b_bdv[i].wallcap!=1 && b_bdv[i].wallvolt!=1) b_flds.phi[bcell] = b_bdv[i].bdphi_val;
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "NEUMANN")
            {
               //  for(j=0;j<b_bdv[i].cboundnum;j++)
               //  {
               //     bcell=b_bdv[i].cboundcells[j];
               //     b_flds.phi[bcell] = b_flds.phi[b_msh.incneighofc(bcell)];
               //  }
               if(b_bdv[i].wallcur==1)
               {
                  for(j=0;j<b_bdv[i].pboundnum;j++)
                  {
                     bcell=b_bdv[i].pboundcells[j];
                     for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.pmesh[b_msh.meshdims*bcell+k];
                     //b_flds.phi[bcell] = beqnparser(b_bdv[i].bdphi,meshpt);  //UPDATE
                     b_bdv[i].bdphi_val =  beqnparser(b_bdv[i].bdphi,meshpt);  //UPDATE
                  }
               }
               else
               {
                  for(j=0;j<b_bdv[i].pboundnum;j++)
                  {
                     bcell=b_bdv[i].pboundcells[j];
                     if(b_msh.meshdims==1) 
                     {
                        meshpt[0] = b_msh.pmesh[bcell];
                        //if(meshpt[0]<b_msh.meshbd[0]) b_flds.phi[bcell] =  b_flds.phi[b_msh.incneighofc(bcell)] - (b_msh.meshbd[0]-meshpt[0])*beqnparser(b_bdv[i].bdphi,meshpt); //UPDATE
                        //else if(meshpt[0]>b_msh.meshbd[1]) b_flds.phi[bcell] =  b_flds.phi[b_msh.incneighofc(bcell)] + (meshpt[0]-b_msh.meshbd[1])*beqnparser(b_bdv[i].bdphi,meshpt); //UPDATE
                        b_bdv[i].bdphi_val = b_flds.phi[bcell];
                     }
                  }
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "PERIODIC")
            {
               for(j=0;j<b_bdv[i].pboundnum;j++)
               {
                  bcell=b_bdv[i].pboundcells[j];
                  if(b_msh.meshdims==1)      //  1D
                  {
                     meshpt[0] = b_msh.pmesh[b_msh.meshdims*bcell];
                     //if(meshpt[0]<b_msh.meshbd[j])
                     if(meshpt[0]<=b_msh.meshbd[0])
                     {  
                        meshptin[0] = b_msh.pmesh[b_msh.meshdims*(bcell+b_msh.numghost)];
                        b_flds.phi[bcell] = beqnparser(b_bdv[0].bdphi,meshptin); //UPDATE
                        b_flds.phi[b_bdv[1].pboundcells[0]] = b_flds.phi[bcell];
                     }
                  }
               }
            }
            else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
         }


         //..............Bfield....................//

         if(b_bdv[i].boundtype[3*nsp+1] == "DIRICHLET")
         {
            for(j=0;j<b_bdv[i].pboundnum;j++)
            {
               bcell=b_bdv[i].pboundcells[j];
               for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.pmesh[b_msh.meshdims*bcell+k];
               for(k=0;k<b_msh.vecdims;k++) b_flds.B[b_msh.vecdims*bcell+k] = beqnparser(b_bdv[i].bdB[k],meshpt); //UPDATE
            }
         }
         else if(b_bdv[i].boundtype[3*nsp+1] == "NEUMANN")  //CHG
         {
            for(j=0;j<b_bdv[i].pboundnum;j++)
            {
               bcell=b_bdv[i].pboundcells[j];
               //for(k=0;k<b_msh.vecdims;k++)  b_flds.B[b_msh.vecdims*bcell+k] = b_flds.B[b_msh.vecdims*b_msh.incneighofc(bcell)+k];
            }
         }
         else if(b_bdv[i].boundtype[3*nsp+1] == "PERIODIC")
         {
            for(j=0;j<b_bdv[i].pboundnum;j++)
            {
               bcell=b_bdv[i].pboundcells[j];
               /*for(k=0;k<b_msh.meshdims;k++) 
                 {
                 meshpt[k] = b_msh.pmesh[b_msh.meshdims*bcell+k];
                 if(meshpt[k]<b_msh.meshbd[k*b_msh.meshdims])
                 {  
                 for(l=0;l<b_msh.vecdims;l++)
                 {
                 b_flds.B[b_msh.vecdims*(bcell+k)+l] =  b_flds.B[b_msh.vecdims*(bcell+k)+l] + b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j+b_msh.numghost]-b_msh.numghost+k)+l];
                 b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j+b_msh.numghost]-b_msh.numghost+k)+l] = b_flds.B[b_msh.vecdims*(bcell+k)+l];
                 }
                 }
                 else if(meshpt[k]>b_msh.meshbd[k*b_msh.meshdims+1])
                 {  
                 for(l=0;l<b_msh.vecdims;l++)
                 {
                 b_flds.B[b_msh.vecdims*(bcell+k)+l] =  b_flds.B[b_msh.vecdims*(bcell+k)+l] + b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j-b_msh.numghost]+k+b_msh.numghost)+l];
                 b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j-b_msh.numghost]+k+b_msh.numghost)+l] = b_flds.B[b_msh.vecdims*(bcell+k)+l];
                 }
                 }
                 }*/
            }
         }
         else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
      } 
      //std::cout << "x";
   }
   else
   {
      for(i=0;i<nbound;i++)
      {

         //..............Density and Energy....................//
         for(m=0;m<nsp;m++)
         { 
            if(b_bdv[i].boundtype[3*m] == "DIRICHLET")
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];
                  for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
                  b_cont[m].N[bcell] = beqnparser(b_bdv[i].bddens[m],meshpt); //UPDATE
                  b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                  //b_cont[m].En[bcell] = beqnparser(b_bdv[i].bdEn[m],meshpt);
               }
            }
            else if(b_bdv[i].boundtype[3*m] == "NEUMANN")
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];
                  if(b_msh.meshdims==1) 
                  {
                     meshpt[0] = b_msh.cmesh[bcell];
                     if(meshpt[0]<b_msh.meshbd[0]) b_cont[m].N[bcell] =  b_cont[m].N[b_msh.incneighofc(bcell)] - (b_msh.meshbd[0]-meshpt[0])*beqnparser(b_bdv[i].bddens[m],meshpt); //UPDATE
                     //if(meshpt[0]<b_msh.meshbd[0]) b_cont[m].En[bcell] =  b_cont[m].En[b_msh.incneighofc(bcell)] - (b_msh.meshbd[0]-meshpt[0])*beqnparser(b_bdv[i].bdEn[m],meshpt);
                     else if(meshpt[0]>b_msh.meshbd[1]) b_cont[m].N[bcell] =  b_cont[m].N[b_msh.incneighofc(bcell)] + (meshpt[0]-b_msh.meshbd[1])*beqnparser(b_bdv[i].bddens[m],meshpt);  //UPDATE
                     b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                     //else if(meshpt[0]>b_msh.meshbd[1]) b_cont[m].En[bcell] =  b_cont[m].En[b_msh.incneighofc(bcell)] + (meshpt[0]-b_msh.meshbd[1])*beqnparser(b_bdv[i].bdEn[m],meshpt);
                  }
               }
            }
            else if(b_bdv[i].boundtype[3*m] == "PERIODIC") 
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];

                  if(b_msh.meshdims==1) 
                  {
                     //if(b_msh.intscheme==0)
                     //{
                     meshpt[0] = b_msh.cmesh[bcell];
                     if(meshpt[0]<b_msh.meshbd[0])
                     {  
                        b_cont[m].N[bcell] =  b_cont[m].N[b_bdv[1].cboundcells[0]-b_msh.numghost];
                        b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                        b_cont[m].En[bcell] =  b_cont[m].En[b_bdv[1].cboundcells[0]-b_msh.numghost];
                     }
                     else if(meshpt[0]>b_msh.meshbd[1])
                     {   
                        b_cont[m].N[bcell] =  b_cont[m].N[b_bdv[0].cboundcells[0]+b_msh.numghost];
                        b_bdv[i].bddens_val[m] = b_cont[m].N[bcell];
                        b_cont[m].En[bcell] =  b_cont[m].En[b_bdv[0].cboundcells[0]+b_msh.numghost];
                     }
                     //}
                     /*else if(b_msh.intscheme==1)
                       {
                       meshpt[0] = b_msh.cmesh[bcell];
                       if(meshpt[0]<b_msh.meshbd[0])
                       {  
                       b_cont[m].N[bcell] =  b_cont[m].N[bcell] + b_cont[m].N[b_bdv[1].cboundcells[0]-b_msh.numghost];
                       b_cont[m].En[bcell] =  b_cont[m].En[bcell] + b_cont[m].En[b_bdv[1].cboundcells[0]-b_msh.numghost];
                       b_cont[m].N[b_bdv[1].cboundcells[0]-b_msh.numghost] = b_cont[m].N[bcell];
                       b_cont[m].En[b_bdv[1].cboundcells[0]-b_msh.numghost] = b_cont[m].En[bcell];
                       }
                       else if(meshpt[0]>b_msh.meshbd[1])
                       {   
                       b_cont[m].N[bcell] =  b_cont[m].N[bcell] + b_cont[m].N[b_bdv[0].cboundcells[0]+b_msh.numghost];
                       b_cont[m].En[bcell] =  b_cont[m].En[bcell] + b_cont[m].En[b_bdv[0].cboundcells[0]+b_msh.numghost];
                       b_cont[m].N[b_bdv[0].cboundcells[0]+b_msh.numghost] = b_cont[m].N[bcell];
                       b_cont[m].En[b_bdv[0].cboundcells[0]+b_msh.numghost] = b_cont[m].En[bcell];
                       }
                       }*/
                  }
               }
            }
            else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
         }

         //CHG:  Add Other Continuum Variables (U,T)

         if(b_msh.sflag==2)
         {
            //..............E-field....................//

            if(b_bdv[i].boundtype[3*nsp] == "DIRICHLET")
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];
                  for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
                  for(k=0;k<b_msh.vecdims;k++) b_flds.E[b_msh.vecdims*bcell+k] = beqnparser(b_bdv[i].bdE[k],meshpt); //UPDATE
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "NEUMANN") //CHG
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];
                  for(k=0;k<b_msh.vecdims;k++)  b_flds.E[b_msh.vecdims*bcell+k] = b_flds.E[b_msh.vecdims*b_msh.incneighofc(bcell)+k];
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "PERIODIC")
            {
            }
            else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
         }
         else if(b_msh.sflag==0)
         {
            //..............Potential....................//

            if(b_bdv[i].boundtype[3*nsp] == "DIRICHLET")
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];

                  if(b_msh.meshdims==1) 
                  {
                     meshpt[0] = b_msh.cmesh[bcell];
                     b_bdv[i].bdphi_val = beqnparser(b_bdv[i].bdphi,meshpt);  //UPDATE
                     if(meshpt[0]<b_msh.meshbd[0]) b_flds.phi[bcell] =  2.0*b_bdv[i].bdphi_val-b_flds.phi[b_msh.incneighofc(bcell)]; //UPDATE
                     if(meshpt[0]>b_msh.meshbd[0]) b_flds.phi[bcell] =  2.0*b_bdv[i].bdphi_val-b_flds.phi[b_msh.incneighofc(bcell)]; //UPDATE
                  }
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "NEUMANN")
            {
               //  for(j=0;j<b_bdv[i].cboundnum;j++)
               //  {
               //     bcell=b_bdv[i].cboundcells[j];
               //     b_flds.phi[bcell] = b_flds.phi[b_msh.incneighofc(bcell)];
               //  }
               if(b_bdv[i].wallcur==1)
               {
                  for(j=0;j<b_bdv[i].cboundnum;j++)
                  {
                     bcell=b_bdv[i].cboundcells[j];
                     for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
                     b_flds.phi[bcell] = beqnparser(b_bdv[i].bdphi,meshpt);  //UPDATE
                     b_bdv[i].bdphi_val = b_flds.phi[bcell];
                  }
               }
               else
               {
                  for(j=0;j<b_bdv[i].cboundnum;j++)
                  {
                     bcell=b_bdv[i].cboundcells[j];
                     if(b_msh.meshdims==1) 
                     {
                        meshpt[0] = b_msh.cmesh[bcell];
                        if(meshpt[0]<b_msh.meshbd[0]) b_flds.phi[bcell] =  b_flds.phi[b_msh.incneighofc(bcell)] - (b_msh.meshbd[0]-meshpt[0])*beqnparser(b_bdv[i].bdphi,meshpt); //UPDATE
                        else if(meshpt[0]>b_msh.meshbd[1]) b_flds.phi[bcell] =  b_flds.phi[b_msh.incneighofc(bcell)] + (meshpt[0]-b_msh.meshbd[1])*beqnparser(b_bdv[i].bdphi,meshpt); //UPDATE
                        b_bdv[i].bdphi_val = b_flds.phi[bcell];
                     }
                  }
               }
            }
            else if(b_bdv[i].boundtype[3*nsp] == "PERIODIC")
            {
               for(j=0;j<b_bdv[i].cboundnum;j++)
               {
                  bcell=b_bdv[i].cboundcells[j];
                  if(b_msh.meshdims==1)      //  1D
                  {
                     meshpt[0] = b_msh.cmesh[b_msh.meshdims*bcell];
                     //if(meshpt[0]<b_msh.meshbd[j])
                     if(meshpt[0]<b_msh.meshbd[0])
                     {  
                        meshptin[0] = b_msh.cmesh[b_msh.meshdims*(bcell+b_msh.numghost)];
                        b_flds.phi[bcell+b_msh.numghost] = beqnparser(b_bdv[0].bdphi,meshptin); //UPDATE
                        b_flds.phi[b_bdv[1].cboundcells[0]] = b_flds.phi[bcell+b_msh.numghost];
                     }
                  }
               }
            }
            else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
         }


         //..............Bfield....................//

         if(b_bdv[i].boundtype[3*nsp+1] == "DIRICHLET")
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
               for(k=0;k<b_msh.vecdims;k++) b_flds.B[b_msh.vecdims*bcell+k] = beqnparser(b_bdv[i].bdB[k],meshpt); //UPDATE
            }
         }
         else if(b_bdv[i].boundtype[3*nsp+1] == "NEUMANN")  //CHG
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.vecdims;k++)  b_flds.B[b_msh.vecdims*bcell+k] = b_flds.B[b_msh.vecdims*b_msh.incneighofc(bcell)+k];
            }
         }
         else if(b_bdv[i].boundtype[3*nsp+1] == "PERIODIC")
         {
            for(j=0;j<b_bdv[i].cboundnum;j++)
            {
               bcell=b_bdv[i].cboundcells[j];
               for(k=0;k<b_msh.meshdims;k++) 
               {
                  meshpt[k] = b_msh.cmesh[b_msh.meshdims*bcell+k];
                  if(meshpt[k]<b_msh.meshbd[k*b_msh.meshdims])
                  {  
                     for(l=0;l<b_msh.vecdims;l++)
                     {
                        b_flds.B[b_msh.vecdims*(bcell+k)+l] =  b_flds.B[b_msh.vecdims*(bcell+k)+l] + b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j+b_msh.numghost]-b_msh.numghost+k)+l];
                        b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j+b_msh.numghost]-b_msh.numghost+k)+l] = b_flds.B[b_msh.vecdims*(bcell+k)+l];
                     }
                  }
                  else if(meshpt[k]>b_msh.meshbd[k*b_msh.meshdims+1])
                  {  
                     for(l=0;l<b_msh.vecdims;l++)
                     {
                        b_flds.B[b_msh.vecdims*(bcell+k)+l] =  b_flds.B[b_msh.vecdims*(bcell+k)+l] + b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j-b_msh.numghost]+k+b_msh.numghost)+l];
                        b_flds.B[b_msh.vecdims*(b_bdv[i].cboundcells[j-b_msh.numghost]+k+b_msh.numghost)+l] = b_flds.B[b_msh.vecdims*(bcell+k)+l];
                     }
                  }
               }
            }
         }
         else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";
      } 
   }

   //std::cout << "x";

}


//..Apply periodic fluid boundary conditions (CFD Project)..//

void boundary::applyPeriodicBoundaryConditions(mesh b_mshflow, flow &b_cont) // TEMPORARY
{
   int i,j,k;

   b_cont.rho[0] = b_cont.rho[(b_cont.rho.size()-4)];
   b_cont.rho[1] = b_cont.rho[(b_cont.rho.size()-3)];
   b_cont.rho[b_cont.rho.size()-2] = b_cont.rho[2];
   b_cont.rho[b_cont.rho.size()-1] = b_cont.rho[3];

   b_cont.U[0] = b_cont.U[(b_cont.U.size()-4)];
   b_cont.U[1] = b_cont.U[(b_cont.U.size()-3)];
   b_cont.U[b_cont.U.size()-2] = b_cont.U[2];
   b_cont.U[b_cont.U.size()-1] = b_cont.U[3];

   b_cont.En[0] = b_cont.En[(b_cont.En.size()-4)];
   b_cont.En[1] = b_cont.En[(b_cont.En.size()-3)];
   b_cont.En[b_cont.En.size()-2] = b_cont.En[2];
   b_cont.En[b_cont.En.size()-1] = b_cont.En[3];

   b_cont.p[0] = b_cont.p[(b_cont.p.size()-4)];
   b_cont.p[1] = b_cont.p[(b_cont.p.size()-3)];
   b_cont.p[b_cont.p.size()-2] = b_cont.p[2];
   b_cont.p[b_cont.p.size()-1] = b_cont.p[3];

}



//..Apply boundary conditions for Euler (CFD project)..//

void boundary::applyBoundaryConditionsEuler(mesh b_msh,flow &b_cont)
{

   // Need to Update to New BC implementation
   /*
      int i,j,k;
      int count = 0;
      int densdcount = 0;
      int pdcount = 0;
      int veldcount = 0;

      std::vector<double> meshpt;

      for(i=0;i<b_msh.meshdims;i++) meshpt.push_back(0.0);

      for(i=0;i<nbound;i++)
      {
   //..............Density....................//

   if(boundtype[3*i] == "DIRICHLET")
   {
   for(j=0;j<b_msh.cboundnum[i];j++)
   {
   for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*b_msh.cboundcells[count+j]+k];
   b_cont.rho[b_msh.cboundcells[count+j]] = beqnparser(bddens[densdcount],meshpt);
   }
   densdcount = densdcount + 1; 
   }
   else if(boundtype[3*i] == "NEUMANN")
   {
   for(j=0;j<b_msh.cboundnum[i];j++)
   {
   for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*b_msh.cboundcells[count+j]+k];
   b_cont.rho[b_msh.cboundcells[count+j]] = b_cont.rho[b_msh.incneighofc(b_msh.cboundcells[count+j])];
   }
   }
   else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";

   //..............Velocity....................//

   if(boundtype[3*i+1] == "DIRICHLET")
   {
   for(j=0;j<b_msh.cboundnum[i];j++)
   {
   for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*b_msh.cboundcells[count+j]+k];
   for(k=0;k<b_msh.vecdims;k++) b_cont.U[b_msh.vecdims*b_msh.cboundcells[count+j]+k] = beqnparser(bdU[b_msh.vecdims*veldcount+k],meshpt);
   }
   veldcount = veldcount + 1; 
   }
   else if(boundtype[3*i+1] == "NEUMANN")
   {
   for(j=0;j<b_msh.cboundnum[i];j++)
   {
   for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*b_msh.cboundcells[count+j]+k];
   for(k=0;k<b_msh.vecdims;k++)  b_cont.U[b_msh.vecdims*b_msh.cboundcells[count+j]+k] = b_cont.U[b_msh.vecdims*b_msh.incneighofc(b_msh.cboundcells[count+j])+k];
   }
   }
   else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";

   //..............Pressure and Energy...................//

   if(boundtype[3*i+2] == "DIRICHLET")
   {
   for(j=0;j<b_msh.cboundnum[i];j++)
   {
   for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*b_msh.cboundcells[count+j]+k];
   b_cont.p[b_msh.cboundcells[count+j]] = beqnparser(bdp[pdcount],meshpt);
   }
   pdcount = pdcount + 1; 
   }
   else if(boundtype[3*i+2] == "NEUMANN")
   {
   for(j=0;j<b_msh.cboundnum[i];j++)
   {
   for(k=0;k<b_msh.meshdims;k++) meshpt[k] = b_msh.cmesh[b_msh.meshdims*b_msh.cboundcells[count+j]+k];
   b_cont.p[b_msh.cboundcells[count+j]] = b_cont.p[b_msh.incneighofc(b_msh.cboundcells[count+j])];
   b_cont.En[b_msh.cboundcells[count+j]] = b_cont.En[b_msh.incneighofc(b_msh.cboundcells[count+j])];
}
}
else std::cout << "\n\n ERROR!!  Boundary Condition not supported\n\n";

count = std::accumulate(b_msh.cboundnum.begin(),b_msh.cboundnum.begin()+i+1,0);
}*/

}


//..Find cells for special regions..//

void boundary::findSpecialRegionCells(mesh &b_msh, std::vector<spclvars> &b_spclv)
{
   int i,j,k;
   int check=1;

   //std::cout << "\n\tDefining Special Regions....\n";

   for(i=0;i<b_spclv[0].nspcl;i++)
   {
      for(j=0;j<b_msh.cmesh.size()/b_msh.meshdims;j++)  
      {
         check=1;
         for(k=0;k<b_msh.meshdims;k++)
         {
            if(b_msh.cmesh[j*b_msh.meshdims+k]<b_spclv[i].spclrange[2*k] || b_msh.cmesh[j*b_msh.meshdims+k]>b_spclv[i].spclrange[2*k+1])  check = 0;
         }
         if(check==1)  b_spclv[i].spclcells.push_back(j);
      }

      for(j=0;j<b_msh.pmesh.size()/b_msh.meshdims;j++)  
      {
         check=1;
         for(k=0;k<b_msh.meshdims;k++)
         {
            if(b_msh.pmesh[j*b_msh.meshdims+k]<b_spclv[i].spclrange[2*k] || b_msh.pmesh[j*b_msh.meshdims+k]>b_spclv[i].spclrange[2*k+1])  check = 0;
         }
         if(check==1)  b_spclv[i].spclpoints.push_back(j);
      }

      //for(j=0;j<b_spclv[i].spclcells.size();j++)  std::cout << b_spclv[i].spclcells[j] << std::endl;
      //for(j=0;j<b_spclv[i].spclpoints.size();j++)  std::cout << b_spclv[i].spclpoints[j] << std::endl;
   }
   //std::cout << "x\n";

}

//..Apply special regions (Source, Heating, Etc)..//

void boundary::applySpecialRegionsPIC(mesh b_msh, std::vector<contnm> &b_cont, std::vector<particles> &b_part, solverVars b_svar, std::vector<spclvars> &b_spclv,fields &b_flds)
{
   int i,j,k,m;
   int scell,spnt;
   int b_nspcl = b_spclv[0].nspcl;
   int b_nsp = b_part[0].nsp;
   int numprocs,procid;  //MPI

   double dens,temp,flux;
   double area,dx,volume;
   long double J_cond,J_cond_check,dEperpdt;
   long double total; //MPI

   int spcl_start,spcl_end;

   J_cond = 0.0;
   J_cond_check = 0.0;
   dEperpdt = 0.0;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   std::vector<int> neigh;
   std::vector<int> spcllocs;
   std::vector<double> vel;
   std::vector<double> meshpt;
   std::vector<double> L;
   std::vector<double> J_cond_vec;
   std::vector<double> dEperpdt_vec;

   for(j=0;j<b_msh.vecdims;j++) vel.push_back(0.0);
   for(j=0;j<b_msh.meshdims;j++) meshpt.push_back(0.0);
   for(j=0;j<2*b_msh.meshdims;j++) neigh.push_back(0);
   for(j=0;j<b_msh.meshdims;j++) L.push_back(0.0);

   solver b_slvr;

   area = 1.0;

   for(i=0;i<b_nspcl;i++)
   {
      if(b_msh.philoc==1)
      {
         for(j=0;j<b_spclv[i].spclpoints.size();j++) spcllocs.push_back(b_spclv[i].spclpoints[j]);
         for(j=0;j<b_spclv[i].spclpoints.size();j++) J_cond_vec.push_back(0.0);
         for(j=0;j<b_spclv[i].spclpoints.size();j++) dEperpdt_vec.push_back(0.0);
      }
      else if(b_msh.philoc==0) 
      {
         for(j=0;j<b_spclv[i].spclcells.size();j++) spcllocs.push_back(b_spclv[i].spclcells[j]);
         for(j=0;j<b_spclv[i].spclcells.size();j++) J_cond_vec.push_back(0.0);
         for(j=0;j<b_spclv[i].spclcells.size();j++) dEperpdt_vec.push_back(0.0);
      }

      if(b_spclv[i].spcltype == "SOURCE")  //..Apply source region
      {
         for(j=0;j<b_nsp;j++)
         {
#if SIMTYPE!=1
            for(k=0;k<spcllocs.size();k++)
            {
               scell = spcllocs[k];
#else
               scell = (b_msh.numpoints[0]-1)/2;
#endif
               for(m=0;m<b_msh.meshdims;m++) neigh[m] = b_msh.cofpneigh[2*scell+m];
               for(m=0;m<b_msh.meshdims;m++) meshpt[m] = b_msh.cmesh[b_msh.meshdims*scell + m]; 
               for(m=0;m<b_msh.vecdims;m++)  vel[m] = 0.0;        //...CHG
               flux = beqnparser(b_spclv[i].spclFlux[j],meshpt);  //UPDATE
               //volume = dx*area;
               dens = flux*b_svar.dt;
               temp = beqnparser(b_spclv[i].spclT[j],meshpt);  //UPDATE

               //std::cout << std::endl << scell << "\t" << flux << "\t" << dens << "\t" << temp ;
               //std::cout << std::endl << scell << "\t" << b_spclv[i].spclddist[j] << "\t" << b_spclv[i].spclthdist[j] << "\t" ;

               b_slvr.particlecellsource(b_part[j], b_msh, b_spclv[i].spclddist[j], b_spclv[i].spclthdist[j],vel,dens,temp,scell); 
#if SIMTYPE!=1
            }
#endif
         }
      }
      else if(b_spclv[i].spcltype == "HEATING")  //..Apply heating region
      {
         for(k=0;k<b_msh.meshdims;k++)  L[k] = b_spclv[i].spclrange[2*k+1]- b_spclv[i].spclrange[2*k];
         spcl_start = spcllocs[0];
         spcl_end = spcllocs[spcllocs.size()-1];


         for(j=0;j<b_nsp;j++)
         {
            for(k=0;k<spcllocs.size();k++) //1D
            {
               scell = spcllocs[k];
               J_cond_vec[k] += b_cont[j].charge*b_cont[j].N[scell]*b_cont[j].U[b_msh.vecdims*scell+2]; //1D,CHG: Current, not current density.
               J_cond += b_cont[j].charge*b_cont[j].N[scell]*b_cont[j].U[b_msh.vecdims*scell+2]; //1D,CHG: Current, not current density.
               //if(procid==0) std::cout << scell <<"\t" << J_cond << "\t";
               //if(k==0) spcl_start = scell;
               //if(k==(spcllocs.size()-1)) spcl_end = scell;
            }

            for(k=0;k<b_part[j].pos.size();k++) //1D
            {
               if(b_msh.philoc==1)
               {
                  if(b_part[j].pos[k]>=b_spclv[i].spclrange[0] && b_part[j].pos[k] <= b_spclv[i].spclrange[1])
                  {
                     J_cond_check += b_part[j].wcharge*b_part[j].vel[b_msh.vecdims*k+2];   //1D,CHG: Current, not current density, assumes flux-tube area is 1.
                     //J_cond_check += b_part[j].wcharge*sqrt(b_part[j].vel[b_msh.vecdims*k+2]*b_part[j].vel[b_msh.vecdims*k+2]);
                  }
               }
               else if(b_msh.philoc==0)
               {
                  if(b_part[j].pos[k]>=b_spclv[i].spclrange[0] && b_part[j].pos[k] <= b_spclv[i].spclrange[1])
                  {
                     J_cond_check += b_part[j].wcharge*b_part[j].vel[b_msh.vecdims*k+2];   //1D,CHG: Current, not current density, assumes flux-tube area is 1.
                     //J_cond_check += b_part[j].wcharge*sqrt(b_part[j].vel[b_msh.vecdims*k+2]*b_part[j].vel[b_msh.vecdims*k+2]);
                  }
               }
            }
         }

         J_cond = J_cond/spcllocs.size();
         J_cond_check = J_cond_check/L[0];
         //J_cond = J_cond/spcllocs.size();

         if(numprocs>1)
         {
            total = 0.0; //MPI
            MPI_Allreduce(&J_cond_check,&total,1, MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD); //MPI
            J_cond_check = total; //MPI

            //total = 0.0; //MPI NOT USED!
            //MPI_Allreduce(&J_cond,&total,1, MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD); //MPI
            //J_cond = total; //MPI
         }

         dEperpdt = (b_spclv[i].spclJ*sin(2.0*3.14*b_spclv[i].spclomega*b_svar.totalTime) - J_cond)/eps0;
         //dEperpdt = (b_spclv[i].spclJ - J_cond)/eps0;

         if(procid==0)  std::cout << "\nCurrents:   " << J_cond << "\t" << J_cond_check << "\t" << b_spclv[i].spclJ*cos(2.0*3.14*b_spclv[i].spclomega*b_svar.totalTime) <<std::endl;
         /*if(std::isnan(J_cond)==true) 
           {
           if(procid==0) std::cout << "\n 1n:\n";
           if(procid==0) for(k=0;k<b_cont[0].N.size();k++) std::cout << "\t" << b_cont[0].N[k]; //1D,CHG: Current, not current density.
           if(procid==0) std::cout << "\n 1U:\n";
           if(procid==0) for(k=0;k<b_cont[0].N.size();k++) std::cout << "\t" << b_cont[0].U[b_msh.vecdims*k+2]; //1D,CHG: Current, not current density.
           if(procid==0) std::cout << "\n 2n:\n";
           if(procid==0) for(k=0;k<b_cont[1].N.size();k++) std::cout << "\t" << b_cont[1].N[k]; //1D,CHG: Current, not current density.
           if(procid==0) std::cout << "\n 2U:\n";
           if(procid==0) for(k=0;k<b_cont[1].N.size();k++) std::cout << "\t" << b_cont[1].U[b_msh.vecdims*k+2]; //1D,CHG: Current, not current density.
           exit(EXIT_FAILURE);
           }*/
         //std::cout << "\nst:   " << spcl_start << "\t" << spcl_end << "\t" << std::endl;

         for(k=0;k<b_spclv[i].spclpoints.size();k++) //1D
         {
            //dEperpdt_vec[k] = (b_spclv[i].spclJ*sin(2.0*3.14*b_spclv[i].spclomega*b_svar.totalTime) - J_cond_vec[k])/eps0;
            spnt = b_spclv[i].spclpoints[k];
            b_flds.E[spnt*b_msh.vecdims+2] += dEperpdt*b_svar.dt;
            //b_flds.E[spnt*b_msh.vecdims+2] += dEperpdt_vec[k]*b_svar.dt;
         }
         //std::cout << "\nEfield:  " << dEperpdt << "\t" << b_flds.E[spnt*b_msh.vecdims+2] << "\t";
      }
      else if(b_spclv[i].spcltype == "EFIELD")
      {
         for(k=0;k<b_spclv[i].spclpoints.size();k++) //1D
         {
            spnt = b_spclv[i].spclpoints[k];
            b_flds.E[spnt*b_msh.vecdims+2] = b_spclv[i].spclE*cos(2.0*3.14158*b_spclv[i].spclomega*b_svar.totalTime);
         }
         std::cout << "\nEfield: " << b_flds.E[spnt*b_msh.vecdims+2] << "\t";
      }

      spcllocs.clear();
      //J_cond_vec.clear();
      //dEperpdt_vec.clear();
   }

}



void boundary::applyEfieldBoundary(mesh &b_msh, fields &b_flds, std::vector<boundvars> &b_bdv)
{
   int i,j,k;

   if(b_msh.philoc==1) // when phi is defined at the cell corners
   {
      for(i=0;i<b_bdv[0].nbound;i++)
      {
         if(b_bdv[i].clss == "WALL" || b_bdv[i].clss=="SOURCE")// || b_bdv[i].wallfloat==1 || b_bdv[i].wallcur==1)
         {
            b_flds.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = b_bdv[i].Ebd;
         }
         //  else if(b_bdv[i].clss=="SOURCE") b_flds.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = 0.0  ;
         else if(b_msh.perflag==1)
         {
            b_flds.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = 0.5*(b_flds.phi[b_msh.pmesh.size()-2] - b_flds.phi[1])/b_msh.deltax;
         }
         else     //For other boundary set wall electric field to zero (or maybe half?)
         {
            if(b_bdv[i].fdir == 1 ) b_flds.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = 0.0;  
            if(b_bdv[i].fdir == -1 ) b_flds.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = 0.0; 
            //if(b_bdv[i].fdir == 1 ) b_fields.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = 0.5*b_field.E[b_msh.vecdims*(b_bdv[i].pboundcells[0]+1)];  
            //if(b_bdv[i].fdir == -1 ) b_fields.E[b_msh.vecdims*b_bdv[i].pboundcells[0]] = 0.5*b_field.E[b_msh.vecdims*(b_bdv[i].pboundcells[0]-1)]; 
         }
      }
   } 
}

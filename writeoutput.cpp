#include "writeoutput.h"

//..Initialize writing..//

writeOutput::writeOutput(int temp1,int nsp)
{
  int i;

  flow_fnum = temp1;

  for(i=0;i<(nsp);i++) 
  {
    part_fnum.push_back(temp1);
    cont_fnum.push_back(temp1);
    field_fnum.push_back(temp1);
    vdf_fnum.push_back(temp1);
    ps_fnum.push_back(temp1);
  } 

  totcont_fnum = 0;
  rst_fnum = 0;
}


//..Write corner point mesh..//
 
void writeOutput::writepMesh(const mesh &w_mesh)
{
   int i,j,k;
   int w_total_cells,fnum;

   std::stringstream fname;

   std::string filename;

   std::ofstream wrtfile("pmesh.out");
   w_total_cells = w_mesh.pmesh.size()/w_mesh.meshdims; 

   if(w_mesh.wrtflag==0) //Debug Write
   {   
     if(w_mesh.meshdims==1)
     {  
        for(i = 0; i<w_total_cells; i++)
        {
           wrtfile << w_mesh.pmesh[w_mesh.meshdims*i] << "\t" << 0 << std::endl;
        }
     }
     if(w_mesh.meshdims==2)
     {  
        for(i = 0; i<w_total_cells; i++)
        {
           wrtfile << w_mesh.pmesh[w_mesh.meshdims*i] << "\t"<< w_mesh.pmesh[w_mesh.meshdims*i+1] << std::endl;
        }
     }
   }
   else if(w_mesh.wrtflag==1) //Tecplot Write
   {

     if(w_mesh.meshdims==1)
     {
       wrtfile << "VARIABLES = \"X\", \"Y\"" << std::endl;
       wrtfile << "ZONE T=\"Mesh\", I=" << w_total_cells << ", J=1, DATAPACKING=POINT" << std::endl;
   
       for(i = 0; i<w_total_cells; i++)
       {
         wrtfile << w_mesh.pmesh[w_mesh.meshdims*i] << "\t" << 0 << std::endl;
       }
     }
     else if(w_mesh.meshdims ==2)
     {
       wrtfile << "VARIABLES = \"X\", \"Y\"" << std::endl;
       wrtfile << "ZONE T=\"Mesh\", I=" << w_mesh.numpoints[0] << ", J=" << w_mesh.numpoints[1] <<", DATAPACKING=POINT" << std::endl;

       for(i = 0; i<w_total_cells; i++)
       {
         wrtfile << w_mesh.pmesh[w_mesh.meshdims*i] << "\t"<< w_mesh.pmesh[w_mesh.meshdims*i+1] << std::endl;
       }
     }
   }

   wrtfile.close();   

}


//..Write cell center mesh..//

void writeOutput::writecMesh(const mesh &w_mesh)
{
   int i,j,k;
   int w_total_cells;

   std::ofstream wrtfile("cmesh.out");

   w_total_cells = w_mesh.cmesh.size()/w_mesh.meshdims; 

   if(w_mesh.wrtflag==0) //Debug Write
   {
     if(w_mesh.meshdims==1)
     {
        for(i = 0; i<w_total_cells; i++)
        {
           wrtfile << w_mesh.cmesh[w_mesh.meshdims*i] << "\t" << 0 << std::endl;
        }
     }
     else if(w_mesh.meshdims==2)
     {
        for(i = 0; i<w_total_cells; i++)
        {
           wrtfile << w_mesh.cmesh[w_mesh.meshdims*i] << "\t"<< w_mesh.cmesh[w_mesh.meshdims*i+1] << std::endl;
        }
     }
   }
   else if(w_mesh.wrtflag==1) //Tecplot Write
   {

     if(w_mesh.meshdims==1)
     {
       wrtfile << "VARIABLES = \"X\", \"Y\"" << std::endl;
       wrtfile << "ZONE T=\"MESH\",  I=" << w_total_cells << ", J=1, DATAPACKING=POINT" << std::endl;
   
        for(i = 0; i<w_total_cells; i++)
        {
           wrtfile << w_mesh.cmesh[w_mesh.meshdims*i] << "\t" << 0 << std::endl;
        }

     }
     else if(w_mesh.meshdims ==2)
     {
       wrtfile << "VARIABLES = \"X\", \"Y\"" << std::endl;
       wrtfile << "ZONE  T=\"MESH\", I=" << (w_mesh.numpoints[0]+1) << ", J=" << (w_mesh.numpoints[1]+1) <<", DATAPACKING=POINT" << std::endl;

       for(i = 0; i<w_total_cells; i++)
       {
         wrtfile << w_mesh.cmesh[w_mesh.meshdims*i] << "\t"<< w_mesh.cmesh[w_mesh.meshdims*i+1] << std::endl;
       }
     }
   }
   
   wrtfile.close();   
}

//..Write face mesh..//

void writeOutput::writefMesh(std::vector<double> w_mesh, std::vector<int> w_num_cells, int w_mesh_dims)
{
   int i,j,k;
   int w_total_faces;

   std::ofstream wrtfile("fmesh.out");

   w_total_faces = w_mesh.size()/w_mesh_dims; 

   for(i = 0; i<w_total_faces; i++)
   {
      wrtfile << w_mesh[w_mesh_dims*i] << "\t"<< w_mesh[w_mesh_dims*i+1] << std::endl;
   }
   
   wrtfile.close();   

}

//..Write particle data (2D)..//

void writeOutput::writeParticles(particles w_part, std::vector<int> w_numpoints, int w_mesh_dims, int w_vec_dims, int w_wrtflag,double time)
{
   int i,j,k;
   int numprocs,procid; //MPI
   int w_total_particles = w_part.pos.size()/w_mesh_dims;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   std::stringstream fname;
   std::string filename;

   std::string pname;

   pname = w_part.name; 

   std::cout << "\n\tWriting Particles...";
   fname << pname.c_str() << "Particles" << part_fnum[w_part.spflag] << ".dat";

   filename = fname.str();

   std::ofstream wrtfile(filename.c_str(),std::ios_base::app | std::ios_base::out);  //MPI 
   //std::ofstream wrtfile(filename.c_str());


   if(w_wrtflag == 0)  //Debug Write
   {
      /* only 1D
     if(w_mesh_dims==1)
     {
     */
        for(i = 0; i<w_total_particles; i++)
        {
           wrtfile << w_part.pos[w_mesh_dims*i] << "\t" << 0.0 << "\t";
           for(j=0;j<w_vec_dims; j++)   wrtfile << w_part.vel[w_vec_dims*i+j] << "\t";
           wrtfile << w_part.en[i] << "\t";
           wrtfile << std::endl;
         }
     /* only 1D
     }
     else if(w_mesh_dims==2)
     {
        for(i = 0; i<w_total_particles; i++)
        {
           for(j=0;j<w_mesh_dims; j++)  wrtfile << w_part.pos[w_mesh_dims*i+j] << "\t";
           for(j=0;j<w_vec_dims; j++)   wrtfile << w_part.vel[w_vec_dims*i+j] << "\t";
           wrtfile << w_part.en[w_mesh_dims*i] << "\t";
           wrtfile << std::endl;
        }
     }
     */
   }
   else if(w_wrtflag == 1)  //Tecplot Write
   {
      /* only 1D
     if(w_mesh_dims==1)
     {
     */
        if(procid==0)  //MPI
        {
          wrtfile << "TITLE=\"" << filename.c_str() << "\""<< std::endl; 
          wrtfile << "VARIABLES = \"X\", \"Y\", \"v1\", \"v2\", \"v3\", \"en\"" << std::endl;
          wrtfile << "ZONE T=\"Full Field\", I=" << w_total_particles << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;
        }  //MPI

        for(i = 0; i<w_total_particles; i++)
        {
           wrtfile << w_part.pos[w_mesh_dims*i] << "\t" << 0.0 << "\t";
           for(j=0;j<w_vec_dims; j++)   wrtfile << w_part.vel[w_vec_dims*i+j] << "\t";
           wrtfile << w_part.en[i] << "\t";
           wrtfile << std::endl;
        }

        /* only 1D
     } 
     else if(w_mesh_dims==2)
     {
        if(procid==0)  //MPI
        {
          wrtfile << "TITLE=\"" << filename.c_str() << "\"" << std::endl; 
          wrtfile << "VARIABLES = \"X\", \"Y\", \"v1\", \"v2\", \"v3\", \"en\"" << std::endl;
          wrtfile << "ZONE T=\"Full Field\", I=" << (w_numpoints[0]+1) << ", J=" << (w_numpoints[1]+1) << ", DATAPACKING=POINT, T=" << time << std::endl;
        }   

        for(i = 0; i<w_total_particles; i++)
        {
           for(j=0;j<w_mesh_dims; j++)  wrtfile << w_part.pos[w_mesh_dims*i+j] << "\t";
           for(j=0;j<w_vec_dims; j++)   wrtfile << w_part.vel[w_vec_dims*i+j] << "\t";
           wrtfile << w_part.en[i] << "\t";
           wrtfile << std::endl;
        }

     } 
     */
   }
   
   wrtfile.close();   
   part_fnum[w_part.spflag] = part_fnum[w_part.spflag]+1;

}

//..Write particle data (1D)..//
/*
void writeOutput::writeParticles(particles w_part, fields w_flds, mesh w_msh, double time)
{
   int i,j,k;
   int w_total_particles = w_part.pos.size()/w_msh.meshdims;
   int numprocs,procid; //MPI
   int pindex;
   double entemp;
   std::vector<int> w_neigh;
   std::vector<double> w_weight;
   std::vector<double> enphi;
   std::vector<double> enE,intE;
   std::vector<double> w_pos;

   std::stringstream fname;
   std::string filename;
   std::string pname;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   for(i=0;i<2; i++) w_pos.push_back(0.0);
   for(i=0;i<2; i++) w_neigh.push_back(0);
   for(i=0;i<2; i++) w_weight.push_back(0.0);
   for(i=0;i<w_msh.vecdims;i++)  intE.push_back(0.0);

   solver w_slvr;

   pname = w_part.name; 

   std::cout << "\n\tWriting Particles...";
   fname << pname.c_str() << "Particles" << part_fnum[w_part.spflag] << ".dat";

   filename = fname.str();
   std::ofstream wrtfile(filename.c_str(),std::ios_base::app | std::ios_base::out);  //MPI 
   //std::ofstream wrtfile(filename.c_str());

   if(w_msh.wrtflag == 0) //Debug Write
   {
     if(w_msh.meshdims==1)
     {
        for(i = 0; i<w_total_particles; i++)
        {
           wrtfile << w_part.pos[w_msh.meshdims*i] << "\t" << 0.0 << "\t";
           for(j=0;j<w_msh.vecdims; j++)   wrtfile << w_part.vel[w_msh.vecdims*i+j] << "\t";
           //wrtfile << w_part.en[i] << "\t" << enphi[i] << "\t" << w_part.en[i]+enphi[i];
           wrtfile << w_part.en[i];
           wrtfile << std::endl; 
        }
     } 
     else if(w_msh.meshdims==2)
     {
        for(i = 0; i<w_total_particles; i++)
        {
           for(j=0;j<w_msh.meshdims; j++)  wrtfile << w_part.pos[w_msh.meshdims*i+j] << "\t";
           for(j=0;j<w_msh.vecdims; j++)   wrtfile << w_part.vel[w_msh.vecdims*i+j] << "\t";
           wrtfile << w_part.en[i] << "\t";
           wrtfile << std::endl;  
        }
     }         
   }
   else if(w_msh.wrtflag == 1) //Tecplot Write
   {
     if(w_msh.meshdims ==1)
     {
        if(procid==0) //MPI
        {
          wrtfile << "TITLE=\"" << filename.c_str() << "\"" << std::endl; 
          wrtfile << "VARIABLES = \"X\", \"Y\", \"v1\", \"v2\", \"v3\", \"en\"" << std::endl;
          wrtfile << "ZONE T=\"Full Field\", I=" << w_total_particles << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;
        }   

        for(i = 0; i<w_total_particles; i++)
        {
           wrtfile << w_part.pos[w_msh.meshdims*i] << "\t" << 0.0 << "\t";
           for(j=0;j<w_msh.vecdims; j++)   wrtfile << w_part.vel[w_msh.vecdims*i+j] << "\t";
           //wrtfile << w_part.en[i] << "\t" << enphi[i];
           wrtfile << w_part.en[i];
           wrtfile << std::endl; 
        }

     }
     else if(w_msh.meshdims==2)
     {
        if(procid==0)
        {
          wrtfile << "TITLE=\"" << filename.c_str() << "\"" << std::endl; 
          wrtfile << "VARIABLES = \"X\", \"Y\", \"v1\", \"v2\", \"v3\", \"en\"" << std::endl;
          wrtfile << "ZONE T=\"Full Field\", I=" << (w_msh.numpoints[0]+1) << ", J=" << (w_msh.numpoints[1]+1) << ", DATAPACKING=POINT" << std::endl;
        }

        for(i = 0; i<w_total_particles; i++)
        {
           for(j=0;j<w_msh.meshdims; j++)  wrtfile << w_part.pos[w_msh.meshdims*i+j] << "\t";
           for(j=0;j<w_msh.vecdims; j++)   wrtfile << w_part.vel[w_msh.vecdims*i+j] << "\t";
           wrtfile << w_part.en[i] << "\t";
           wrtfile << std::endl;  
        }
     } 
   }
   wrtfile.close();   
   part_fnum[w_part.spflag] = part_fnum[w_part.spflag]+1;

}
*/

//....Write with single processor by collecting all particles.....//

void writeOutput::writeParticles(particles w_part, fields w_flds, mesh w_msh, double time)
{
   int i,j,k;
   int numpart,numpart_total,s_seed_proc,numavg,numrem;
   int numprocs,procid; //MPI
   int pindex;
   double entemp;
   std::vector<int> w_neigh;
   std::vector<double> w_weight;
   std::vector<double> enphi;
   std::vector<double> enE,intE;
   std::vector<double> w_pos;

   std::stringstream fname;
   std::string filename;
   std::string pname;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   for(i=0;i<2; i++) w_pos.push_back(0.0);
   for(i=0;i<2; i++) w_neigh.push_back(0);
   for(i=0;i<2; i++) w_weight.push_back(0.0);
   for(i=0;i<w_msh.vecdims;i++)  intE.push_back(0.0);

   std::vector<double> s_pos,s_vel,s_en;
   std::vector<int> s_cell,s_pt;

   std::vector<double> s_pos2,s_vel2,s_en2;
   std::vector<int> s_cell2,s_pt2;

   solver w_slvr;

   pname = w_part.name; 

   std::vector<int> numpart_array;
   std::vector<int> numpos_array,numvel_array,pcnt;
   std::vector<int> disp,meshdisp,vecdisp;
   numpart_total=0;

   for(i=0;i<numprocs;i++) numpart_array.push_back(0);
   for(i=0;i<numprocs;i++) pcnt.push_back(0);

   numpart = w_part.en.size();

   MPI_Barrier(MPI_COMM_WORLD);
   //MPI_Gather(&numpart,1,MPI_INT,&numpart_array.front(),1,MPI_INT,0,MPI_COMM_WORLD);
   MPI_Allgather(&numpart,1,MPI_INT,&numpart_array.front(),1,MPI_INT,MPI_COMM_WORLD);
   MPI_Barrier(MPI_COMM_WORLD);

   MPI_Allreduce(&numpart,&numpart_total,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);

   for(i=0;i<numpart_total;i++)
   {
     for(j=0;j<w_msh.meshdims;j++) s_pos.push_back(0.0);
     for(j=0;j<w_msh.vecdims;j++) s_vel.push_back(0.0);
     s_en.push_back(0.0);
     s_cell.push_back(0);
     s_pt.push_back(0);

     for(j=0;j<w_msh.meshdims;j++) s_pos2.push_back(0.0);
     for(j=0;j<w_msh.vecdims;j++) s_vel2.push_back(0.0);
     s_en2.push_back(0.0);
     s_cell2.push_back(0);
     s_pt2.push_back(0);
   }
   
   pcnt[0] = 0;
   for(i=1;i<numprocs;i++) pcnt[i] = pcnt[i-1] + numpart_array[i-1];

   for(i=0;i<numpart;i++)
   {
     for(j=0;j<w_msh.meshdims;j++) s_pos[(i+pcnt[procid])*w_msh.meshdims+j] = w_part.pos[i*w_msh.meshdims+j];
     for(j=0;j<w_msh.vecdims;j++) s_vel[(i+pcnt[procid])*w_msh.vecdims+j] = w_part.vel[i*w_msh.vecdims+j];
     s_en[i+pcnt[procid]] = w_part.en[i];
     s_cell[i+pcnt[procid]] = w_part.cell[i];
     s_pt[i+pcnt[procid]] = w_part.pt[i];
   }

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Allreduce(&s_pos.front(),&s_pos2.front(),numpart_total*w_msh.meshdims, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_vel.front(),&s_vel2.front(),numpart_total*w_msh.vecdims, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_en.front(),&s_en2.front(),numpart_total, MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_cell.front(),&s_cell2.front(),numpart_total, MPI_INT,MPI_SUM,MPI_COMM_WORLD);
   MPI_Allreduce(&s_pt.front(),&s_pt2.front(),numpart_total, MPI_INT,MPI_SUM,MPI_COMM_WORLD);

   s_pos.clear();

   if(procid==0)
   {
     std::cout << "\n\tWriting Particles...";
     fname << pname.c_str() << "Particles" << part_fnum[w_part.spflag] << ".dat";

     filename = fname.str();
     std::ofstream wrtfile(filename.c_str(),std::ios_base::app | std::ios_base::out);  //MPI 
     //std::ofstream wrtfile(filename.c_str());

     if(w_msh.wrtflag == 0) //Debug Write
     {
       if(w_msh.meshdims==1)
       {
          for(i = 0; i<numpart_total; i++)
          {
             wrtfile << s_pos2[w_msh.meshdims*i] << "\t" << 0.0 << "\t";
             for(j=0;j<w_msh.vecdims; j++)   wrtfile << s_vel2[w_msh.vecdims*i+j] << "\t";
             //wrtfile << w_part.en[i] << "\t" << enphi[i] << "\t" << w_part.en[i]+enphi[i];
             wrtfile << s_en2[i];
             wrtfile << std::endl; 
          }
       }  
     } 
     else if(w_msh.wrtflag == 1) //Tecplot Write
     {
       if(w_msh.meshdims ==1)
       {
          if(procid==0) //MPI
          {
            wrtfile << "TITLE=\"" << filename.c_str() << "\"" << std::endl; 
            wrtfile << "VARIABLES = \"X\", \"Y\", \"v1\", \"v2\", \"v3\", \"en\"" << std::endl;
            wrtfile << "ZONE T=\"Full Field\", I=" << numpart_total << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;
          }   

          for(i = 0; i<numpart_total; i++)
          {
             wrtfile << s_pos2[w_msh.meshdims*i] << "\t" << 0.0 << "\t";
             for(j=0;j<w_msh.vecdims; j++)   wrtfile << s_vel2[w_msh.vecdims*i+j] << "\t";
             //wrtfile << w_part.en[i] << "\t" << enphi[i];
             wrtfile << s_en2[i];
             wrtfile << std::endl; 
          }

       }
     } 
     wrtfile.close();   
   }
   part_fnum[w_part.spflag] = part_fnum[w_part.spflag]+1;

}

//..Write single particle..//

void writeOutput::writeSingleParticle(particles w_part, fields w_flds, mesh w_msh, double time)
{
   int i,j,k;
   int w_total_particles = w_part.pos.size()/w_msh.meshdims;
   int pindex;
   double enphi,entot;
   std::vector<int> w_neigh;
   std::vector<double> w_weight;

   std::stringstream fname;
   std::string filename;
   std::string pname;

   pname = w_part.name; 

   for(i=0;i<2*w_msh.meshdims;i++) w_neigh.push_back(0);

   std::cout << "\n\tWriting Particle...";
   fname << pname.c_str() << "Particle" << ".dat";

   filename = fname.str();

   std::ofstream wrtfile(filename.c_str(),std::ios_base::app | std::ios_base::out); 

   if(w_msh.intscheme == 0)
   {
      pindex = w_msh.nearc(w_part,0);
      enphi = w_part.charge*w_flds.phi[pindex];  //DC charge
   }
   else if(w_msh.intscheme == 1)
   {
      pindex = w_msh.nearp(w_part,0);
      w_neigh = w_msh.cneighofp(pindex);

      w_weight = w_msh.cell2part_lineweight(w_part.pos[0],w_neigh);

      enphi = w_part.wcharge*(w_weight[0]*w_flds.phi[w_neigh[0]] + w_weight[1]*w_flds.phi[w_neigh[1]]); //DC charge
   }

   entot = w_part.en[0] + enphi;

   if(w_msh.meshdims==1)
   {
       wrtfile << w_part.pos[0] << "\t" << 0.0 << "\t";
       for(j=0;j<w_msh.vecdims; j++)   wrtfile << w_part.vel[j] << "\t";
       wrtfile << w_part.en[0] << "\t" << enphi << "\t" << entot <<  "\t" << time;
       wrtfile << std::endl; 
   } 
   else if(w_msh.meshdims==2)
   {
       for(j=0;j<w_msh.meshdims; j++)  wrtfile << w_part.pos[j] << "\t";
       for(j=0;j<w_msh.vecdims; j++)   wrtfile << w_part.vel[j] << "\t";
       wrtfile << w_part.en[0] << "\t" << time;
       wrtfile << std::endl;  
   }
   
   wrtfile.close();   
}


//..Write field output..//

void writeOutput::writeFields(fields w_EM, int w_mesh_dims, int w_vec_dims)
{
   int i,j,k;
   int w_total_points;

   std::ofstream wrtEfile("Efield.out");

   w_total_points = w_EM.E.size()/w_mesh_dims/w_vec_dims;

   for(i = 0; i<w_total_points; i++)
   {
      for(j=0;j<w_vec_dims; j++)
      {
         wrtEfile << w_EM.E[w_mesh_dims*i+j] << "\t";
      }
      wrtEfile << std::endl;   
   }
   wrtEfile.close();   

   std::ofstream wrtBfile("Bfield.out");

   w_total_points = w_EM.B.size()/w_mesh_dims/w_vec_dims;

   for(i = 0; i<w_total_points; i++)
   {
      for(j=0;j<w_vec_dims; j++)
      {
         wrtBfile << w_EM.B[w_mesh_dims*i+j] << "\t";
      }
      wrtBfile << std::endl;   
   }
   wrtBfile.close();
}


//..Write Euler output..//

void writeOutput::writeFlowEuler(mesh w_msh, flow w_cont)
{
   int i,j,k;
   int w_total_cells;

   std::stringstream fname;
   std::string filename;

   std::string pname;
 
   std::cout << "\n\tWriting Euler Flow Field..." ;

   pname = "flowfield"; 

   fname << pname.c_str() << flow_fnum << ".dat";

   filename = fname.str();

   std::ofstream wrtfile(filename.c_str());

   wrtfile.precision(15);

   w_total_cells = w_msh.cmesh.size()/w_msh.meshdims; 

   for(i = 0; i<w_total_cells; i++)
   {
      for(j=0;j<w_msh.meshdims;j++) wrtfile << w_msh.cmesh[w_msh.meshdims*i+j] << "\t";
      wrtfile << w_cont.rho[i] << "\t";
      wrtfile << w_cont.p[i] << "\t";
      for(j=0;j<w_msh.vecdims;j++)   wrtfile << w_cont.U[w_msh.vecdims*i+j] << "\t";
      wrtfile << w_cont.En[i] << "\t";
      wrtfile << std::endl;
   }
    
   cont_fnum[w_cont.spflag] = cont_fnum[w_cont.spflag]+1;

   wrtfile.close();
}


//..Write continuum and field output for PIC..//

void writeOutput::writePICField(mesh w_msh, contnm w_cont,fields w_EM,double time)
{
   int i,j,k;
   int w_total_ccells, w_total_pcells;
   std::vector<double>  w_intE,w_intB;

   for(i=0;i<w_msh.vecdims;i++) w_intE.push_back(0.0);
   for(i=0;i<w_msh.vecdims;i++) w_intB.push_back(0.0);

   std::stringstream cfname;
   std::stringstream pfname;
   std::string cfilename;
   std::string pfilename;

   std::string cname;
   std::string pname;

   std::cout << "\n\tWriting PIC Flow Field...." ;

   cname = "Output_cField";
   pname = "Output_pField";

   cfname << w_cont.name << cname.c_str() << cont_fnum[w_cont.spflag] << ".dat";
   pfname << w_cont.name << pname.c_str() << cont_fnum[w_cont.spflag] << ".dat";

   cfilename = cfname.str();
   pfilename = pfname.str();

   std::ofstream cwrtfile(cfilename.c_str());
   std::ofstream pwrtfile(pfilename.c_str());

   cwrtfile.precision(15);
   pwrtfile.precision(15);

   w_total_ccells = w_msh.cmesh.size()/w_msh.meshdims; 
   w_total_pcells = w_msh.pmesh.size()/w_msh.meshdims; 

   if(w_msh.wrtflag==0) //Debug Write
   {
     if(w_msh.philoc==0 && w_msh.Bloc==0 && w_msh.Eloc==0)
     {
        for(i = 0; i<w_total_ccells; i++)
        {
           cwrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           cwrtfile << w_cont.N[i] << "\t";
           cwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.T[i] << "\t";
           cwrtfile << w_cont.En[i] << "\t";
           cwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_EM.phi[i] << "\t";
           cwrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << std::endl;
        }
     }
     else if(w_msh.philoc==0 && w_msh.Bloc==0 && w_msh.Eloc==1)
     {
        for(i = 0; i<w_total_ccells; i++)
        {
           cwrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           cwrtfile << w_cont.N[i] << "\t";
           cwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.T[i] << "\t";
           cwrtfile << w_cont.En[i] << "\t";
           //w_intE = w_msh.extravec_p2c(w_EM.E,i); 
           //cwrtfile << w_intE[0] << "\t"<< w_intE[1] << "\t" << w_intE[2] << "\t";
           cwrtfile << w_EM.phi[i] << "\t";
           cwrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << std::endl;
        }
        for(i = 0; i<w_total_pcells; i++)
        {
           pwrtfile << w_msh.pmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           //pwrtfile << w_cont.N[i] << "\t";
           pwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           //pwrtfile << w_EM.phi[i] << "\t";
           //wrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           pwrtfile << std::endl;
        }
     }
     else if(w_msh.philoc==0 && w_msh.Bloc==1)
     {
        for(i = 0; i<w_total_ccells; i++)
        {
           cwrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           cwrtfile << w_cont.N[i] << "\t";
           cwrtfile << w_cont.En[i] << "\t";
           cwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.T[i] << "\t";
           //cwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_EM.phi[i] << "\t";
           //w_intB = w_msh.extravec_p2c(w_EM.E,i); 
           //wrtfile << w_intB[0] << "\t"<< w_intB[1] << "\t" << w_intB[2] << "\t";
           cwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << std::endl;
        }
      }
     else if(w_msh.philoc==1 && w_msh.Bloc==1 && w_msh.Eloc==1)
     {
        for(i = 0; i<w_total_pcells; i++)
        {
           pwrtfile << w_msh.pmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           pwrtfile << w_cont.N[i] << "\t";
           pwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_cont.T[i] << "\t";
           pwrtfile << w_cont.En[i] << "\t";
           pwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_EM.phi[i] << "\t";
           pwrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           pwrtfile << std::endl;
        }
     }
   }
   else if(w_msh.wrtflag==1) //Tecplot Write
   {
     if(w_msh.philoc==0 && w_msh.Bloc==0 && w_msh.Eloc==0)
     {
        cwrtfile << "TITLE=\"" << cfilename.c_str() << "\"" << std::endl; 
        cwrtfile << "VARIABLES = \"X\", \"Y\", \"N\", \"U\", \"V\", \"W\", \"Tx\", \"Ty\", \"Tz\", \"T\", \"En\", \"Ex\", \"Ey\", \"Ez\", \"phi\", \"Bx\", \"By\", \"Bz\", \"Qx\", \"Qy\", \"Qz\"" << std::endl;
        cwrtfile << "ZONE T=\"Full Field\", I=" << w_total_ccells << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;
   
        for(i = 0; i<w_total_ccells; i++)
        {
           cwrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           cwrtfile << w_cont.N[i] << "\t";
           cwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.T[i] << "\t";
           cwrtfile << w_cont.En[i] << "\t";
           cwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_EM.phi[i] << "\t";
           cwrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << std::endl;
        }

     }
     else if(w_msh.philoc==0 && w_msh.Bloc==0 && w_msh.Eloc==1)
     {
        cwrtfile << "TITLE=\"" << cfilename.c_str() << "\"" << std::endl; 
        cwrtfile << "VARIABLES = \"X\", \"Y\", \"N\", \"U\", \"V\", \"W\", \"Tx\", \"Ty\", \"Tz\", \"T\", \"En\", \"phi\", \"Bx\", \"By\", \"Bz\", \"Qx\", \"Qy\", \"Qz\"" << std::endl;
        cwrtfile << "ZONE T=\"Full Field\", I=" << w_total_ccells << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;
   
        for(i = 0; i<w_total_ccells; i++)
        {
           cwrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           cwrtfile << w_cont.N[i] << "\t";
           cwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.T[i] << "\t";
           cwrtfile << w_cont.En[i] << "\t";
           //w_intE = w_msh.extravec_p2c(w_EM.E,i); 
           //cwrtfile << w_intE[0] << "\t"<< w_intE[1] << "\t" << w_intE[2] << "\t";
           cwrtfile << w_EM.phi[i] << "\t";
           cwrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           cwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           cwrtfile << std::endl;
        }
               
        pwrtfile << "TITLE=\"" << pfilename.c_str() << "\"" << std::endl; 
        pwrtfile << "VARIABLES = \"X\", \"Y\", \"Ex\", \"Ey\", \"Ez\"" << std::endl;
        pwrtfile << "ZONE T=\"Full Field\", I=" << w_total_pcells << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time<< std::endl;
   
        for(i = 0; i<w_total_pcells; i++)
        {
           pwrtfile << w_msh.pmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           //pwrtfile << w_cont.N[i] << "\t";
           pwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           //pwrtfile << w_EM.phi[i] << "\t";
           //wrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           pwrtfile << std::endl;
        }

     }
     else if(w_msh.philoc==1 && w_msh.Bloc==1 && w_msh.Eloc==1)
     {
        pwrtfile << "TITLE=\"" << pfilename.c_str() << "\"" << std::endl; 
        pwrtfile << "VARIABLES = \"X\", \"Y\", \"N\", \"U\", \"V\", \"W\", \"Tx\", \"Ty\", \"Tz\", \"T\", \"En\", \"Ex\", \"Ey\", \"Ez\", \"phi\", \"Bx\", \"By\", \"Bz\", \"Qx\", \"Qy\", \"Qz\"" << std::endl;
        pwrtfile << "ZONE T=\"Full Field\", I=" << w_total_pcells << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;
   
        for(i = 0; i<w_total_pcells; i++)
        {
           pwrtfile << w_msh.pmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
           pwrtfile << w_cont.N[i] << "\t";
           pwrtfile << w_cont.U[w_msh.vecdims*i] << "\t"<< w_cont.U[w_msh.vecdims*i+1] << "\t" << w_cont.U[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_cont.Tvec[w_msh.vecdims*i] << "\t"<< w_cont.Tvec[w_msh.vecdims*i+1] << "\t" << w_cont.Tvec[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_cont.T[i] << "\t";
           pwrtfile << w_cont.En[i] << "\t";
           pwrtfile << w_EM.E[w_msh.vecdims*i] << "\t"<< w_EM.E[w_msh.vecdims*i+1] << "\t" << w_EM.E[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_EM.phi[i] << "\t";
           pwrtfile << w_EM.B[w_msh.vecdims*i] << "\t"<< w_EM.B[w_msh.vecdims*i+1] << "\t" << w_EM.B[w_msh.vecdims*i+2] << "\t";
           pwrtfile << w_cont.Qvec[w_msh.vecdims*i] << "\t"<< w_cont.Qvec[w_msh.vecdims*i+1] << "\t" << w_cont.Qvec[w_msh.vecdims*i+2] << "\t";
           pwrtfile << std::endl;
        }

     }
   }
  
   cont_fnum[w_cont.spflag] = cont_fnum[w_cont.spflag]+1;

   cwrtfile.close();
   pwrtfile.close();
}

//..Write energies output for PIC..//

void writeOutput::writePICField(mesh w_msh, std::vector<contnm> w_cont,fields w_EM, double time)  //FIX,GS
{
   int i,j,k;
   int w_total_ccells, w_total_pcells, w_total_cells;
   std::vector<double>  w_intE,w_intB,rho,en_ke,en_phi,en_tot;
   double dx,en_E;

   for(i=0;i<w_msh.vecdims;i++) w_intE.push_back(0.0);
   for(i=0;i<w_msh.vecdims;i++) w_intB.push_back(0.0);

   std::stringstream cfname;
   std::stringstream pfname;
   std::stringstream fname;
   std::string cfilename;
   std::string pfilename;
   std::string filename;

   std::string cname;
   std::string pname;
   std::string name;

   std::cout << "\n\tWriting PIC Flow Field...." ;

   //cname = "Output_cField";
   //pname = "Output_pField";
   name = "Output_Field";

   //cfname <<  cname.c_str() << totcont_fnum << ".dat";
   //pfname <<  pname.c_str() << totcont_fnum << ".dat";
   fname <<  name.c_str() << totcont_fnum << ".dat";

   //cfilename = cfname.str();
   //pfilename = pfname.str();
   filename = fname.str();

   //std::ofstream cwrtfile(cfilename.c_str());
   //std::ofstream pwrtfile(pfilename.c_str());
   std::ofstream wrtfile(filename.c_str());

   //cwrtfile.precision(15);
   //pwrtfile.precision(15);
   wrtfile.precision(15);

   w_total_ccells = w_msh.cmesh.size()/w_msh.meshdims; 
   w_total_pcells = w_msh.pmesh.size()/w_msh.meshdims; 

   dx = w_msh.cmesh[1]-w_msh.cmesh[0];

   if(w_msh.philoc==1) w_total_cells=w_total_pcells;
   else if(w_msh.philoc==0) w_total_cells=w_total_ccells;

   for(i=0;i<w_total_cells;i++)
   {
     rho.push_back(0.0);
     en_phi.push_back(0.0);
     en_ke.push_back(0.0);
     en_tot.push_back(0.0);
   }

   for(j=0;j<w_cont[0].nsp;j++)
   {
     for(i = 0; i<w_total_cells; i++)
     {
       rho[i] = rho[i] + w_cont[j].charge*w_cont[j].N[i];
       en_phi[i] = en_phi[i] + dx*0.5*w_cont[j].charge*w_cont[j].N[i]*w_EM.phi[i]; //DC;
       en_ke[i] = en_ke[i] + w_cont[j].En[i];
       en_tot[i] = en_ke[i] + en_phi[i];
     }
   }

   if(w_msh.wrtflag==0)
   {
     for(i = 0; i<w_total_cells; i++)
     {
       if(w_msh.philoc==1) wrtfile << w_msh.pmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
       else if(w_msh.philoc==0) wrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
       if(w_cont.size()>w_cont[0].nsp) wrtfile << w_cont[w_cont[0].nsp].N[i] << "\t";
       else wrtfile << 0.0 << "\t";
       wrtfile << rho[i] << "\t";
       wrtfile << en_phi[i] << "\t";
       wrtfile << en_ke[i] << "\t";
       wrtfile << en_tot[i] << "\t";
       wrtfile << std::endl;
     }
   }
   else if(w_msh.wrtflag==1)
   {
     wrtfile << "TITLE=\"" << cfilename.c_str() << "\"" << std::endl; 
     wrtfile << "VARIABLES = \"X\", \"Y\", \"N\", \"rho\", \"EnPhi\", \"EnKe\", \"EnTot\"" << std::endl;
     wrtfile << "ZONE T=\"Full Field\", I=" << w_total_ccells << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;

     for(i = 0; i<w_total_cells; i++)
     {
       if(w_msh.philoc==1) wrtfile << w_msh.pmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
       else if(w_msh.philoc==0) wrtfile << w_msh.cmesh[w_msh.meshdims*i] <<"\t"<< 0.0 << "\t";
       if(w_cont.size()>w_cont[0].nsp) wrtfile << w_cont[w_cont[0].nsp].N[i] << "\t";
       else wrtfile << 0.0 << "\t";
       wrtfile << rho[i] << "\t";
       wrtfile << en_phi[i] << "\t";
       wrtfile << en_ke[i] << "\t";
       wrtfile << en_tot[i] << "\t";
       wrtfile << std::endl;
     }

     /*pwrtfile << "TITLE=\"" << pfilename.c_str() << "\"" << std::endl; 
     pwrtfile << "VARIABLES = \"X\", \"Y\", \"EnE\"" << std::endl;
     pwrtfile << "ZONE T=\"Full Field\", I=" << w_total_pcells << ", J=1, DATAPACKING=POINT, STRANDID=1, SOLUTIONTIME=" << time << std::endl;

     for(i = 0; i<w_total_pcells; i++)
     {
       en_E = 0.0;
       pwrtfile << w_msh.pmesh[w_msh.meshdims*i] << "\t" << 0.0 << "\t";
       for(j=0;j<w_msh.vecdims;j++) en_E = en_E + (0.5)*eps0*w_EM.E[w_msh.vecdims*i+j]*w_EM.E[w_msh.vecdims*i+j]; 
       pwrtfile << en_E;
       pwrtfile << std::endl;
     }*/
    
     /*for(i = 0; i<w_total_ccells; i++)
     {
        cwrtfile <<  0.0 << "\t";
     }
     cwrtfile << std::endl;

     for(i = 0; i<w_total_ccells; i++)
     {
        cwrtfile << rho[i] << "\t";
     } 
     cwrtfile << std::endl;

     for(i = 0; i<w_total_ccells; i++)
     {
        cwrtfile << en_phi[i] << "\t";
     }
     cwrtfile << std::endl;

     for(i = 0; i<w_total_ccells; i++)
     {
        cwrtfile << en_ke[i] << "\t";
     }
     cwrtfile << std::endl;

     for(i = 0; i<w_total_ccells; i++)
     {
        cwrtfile << en_tot[i] << "\t";
     }
     cwrtfile << std::endl;*/

   }
    
   totcont_fnum = totcont_fnum+1;

   wrtfile.close();
   //cwrtfile.close();
   //pwrtfile.close();
}


//..Find global energies..//

void writeOutput::findglobalenergy(const std::vector<particles> &s_part, const fields &s_EM, const mesh &s_msh, double time)
{
   int i,j,k,l,s_index;
   double s_cellvolume,glob_EMEn,avg_E,dx;
   std::vector<int> s_neighbors;
   std::vector<double> s_arearatio;
   std::vector<double> glob_En,glob_meanKEn,glob_thermEn,glob_U;
   double glob_En_sum, glob_meanKEn_sum,glob_thermEn_sum,glob_EnTot;

   int numprocs,procid; //MPI

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI


   for(i=0;i<2; i++) s_neighbors.push_back(0);
   for(i=0;i<2; i++) s_arearatio.push_back(0.0);
   for(i=0;i<s_part[0].nsp*s_msh.vecdims; i++) glob_U.push_back(0.0);
   for(i=0;i<s_part[0].nsp; i++) glob_En.push_back(0.0);
   for(i=0;i<s_part[0].nsp; i++) glob_meanKEn.push_back(0.0);
   for(i=0;i<s_part[0].nsp; i++) glob_thermEn.push_back(0.0);

   //std::cout << "\n\tFinding Global Energies...";


   for(k=0; k<s_part[0].nsp;k++)  //CHGMPI
   {
      for(i=0;i<(s_part[k].pos.size()/s_msh.meshdims);i++)
      {
        glob_En[k] = glob_En[k] + s_part[k].pwght*s_part[k].en[i];
        for(j=0; j<s_msh.vecdims; j++)  glob_U[s_msh.vecdims*k+j] = glob_U[s_msh.vecdims*k+j] + s_part[k].vel[s_msh.vecdims*i+j];
      }
      for(j=0; j<s_msh.vecdims; j++)  glob_U[s_msh.vecdims*k+j] = glob_U[s_msh.vecdims*k+j]/(s_part[k].pos.size()/s_msh.meshdims);
      for(j=0; j<s_msh.vecdims; j++)  glob_meanKEn[k] = glob_meanKEn[k] + 0.5*s_part[k].wmass*s_part[k].pos.size()*glob_U[s_msh.vecdims*k+j]*glob_U[s_msh.vecdims*k+j];
      glob_thermEn[k] = glob_En[k] - glob_meanKEn[k]; 
      //wrtfile << glob_En[k] << "\t" << glob_meanKEn[k] << "\t" << glob_thermEn[k] << "\t";
   }
 


   glob_EMEn = 0.0;
   avg_E = 0.0;
   dx = s_msh.cmesh[1]-s_msh.cmesh[0];

   if(s_msh.Eloc==0)
   {
     for(i=0;i<(s_msh.cmesh.size()/s_msh.meshdims-1);i++)
     {
       for(j=0;j<s_msh.vecdims;j++) glob_EMEn = glob_EMEn + s_EM.E[i*s_msh.vecdims+j]*s_EM.E[i*s_msh.vecdims+j];
       for(j=0;j<s_msh.vecdims;j++) avg_E = avg_E + fabs(s_EM.E[i*s_msh.vecdims+j]);
     }
     avg_E = avg_E/(s_msh.cmesh.size()/s_msh.meshdims-1);
     glob_EMEn = 0.5*glob_EMEn*dx*eps0;
   } 
   else if(s_msh.Eloc==1)
   {
     for(i=0;i<(s_msh.pmesh.size()/s_msh.meshdims-1);i++)
     {
       for(j=0;j<s_msh.vecdims;j++) glob_EMEn = glob_EMEn + s_EM.E[i*s_msh.vecdims+j]*s_EM.E[i*s_msh.vecdims+j];
       for(j=0;j<s_msh.vecdims;j++) avg_E = avg_E + fabs(s_EM.E[i*s_msh.vecdims+j]);
     }
     avg_E = avg_E/(s_msh.pmesh.size()/s_msh.meshdims-1);
     glob_EMEn = 0.5*glob_EMEn*dx*eps0;
   } 

   glob_En_sum = std::accumulate(glob_En.begin(),glob_En.end(),0.0);
   glob_meanKEn_sum = std::accumulate(glob_meanKEn.begin(),glob_meanKEn.end(),0.0);
   glob_thermEn_sum =  std::accumulate(glob_thermEn.begin(),glob_thermEn.end(),0.0); 

   glob_EnTot = glob_En_sum+glob_EMEn;

   double glob_En_sum_MPI,glob_meanKEn_sum_MPI,glob_thermEn_sum_MPI,glob_EnTot_MPI;

   MPI_Reduce(&glob_En_sum,&glob_En_sum_MPI,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&glob_meanKEn_sum,&glob_meanKEn_sum_MPI,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   MPI_Reduce(&glob_thermEn_sum,&glob_thermEn_sum_MPI,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   //MPI_Reduce(&glob_EnTot,&glob_EnTot_MPI,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);


   if(procid==0)
   {
     glob_EnTot_MPI = glob_En_sum_MPI+glob_EMEn;
     std::ofstream wrtfile("energyhistory.dat",std::ios_base::app | std::ios_base::out); 
     wrtfile << time << "\t";
     wrtfile << glob_En_sum_MPI << "\t" << glob_meanKEn_sum_MPI << "\t" << glob_thermEn_sum_MPI << "\t" << glob_EMEn << "\t" << avg_E << "\t" << glob_EnTot_MPI;
     wrtfile << std::endl;

     wrtfile.close();
   }
   //std::cout << "x";


}


//..Find Velocity distribution..//

void writeOutput::findvdf(particles s_part, mesh s_msh, int s_index, double time)
{
   int i,j,k;
   int ind;
   double maxen,maxvel,maxvelcheck,minvelcheck,maxvmag;
   std::vector<double> pdf, velocity;
   std::vector<double> pdfx, velocityx;
   std::vector<double> pdfy, velocityy;
   std::vector<double> pdfz, velocityz;
   double vel,numpart,numpartx,delvel,delvelx;
   double numparty,delvely;
   double numpartz,delvelz;
   double velnd;  //CV
   int numvels;
   std::stringstream fname;
   std::string filename,pname;
   long int np = s_part.pos.size()/s_msh.meshdims;
   long int npg  = np - s_part.gp;

   int numprocs,procid; //MPI

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI


   //s_msh.connectpartandmesh(s_part);
   maxen = *std::max_element(s_part.en.begin(),s_part.en.end());
   //maxvel = sqrt(maxen*2.0/s_part.mass);
   maxvelcheck = *std::max_element(s_part.vel.begin(),s_part.vel.end());
   minvelcheck = fabs(*std::min_element(s_part.vel.begin(),s_part.vel.end()));

   if(maxen==0) maxen = 1e-30;

   maxvmag = sqrt(maxen*2.0/s_part.mass);
   //std::cout << std::endl <<  maxen << "\t" << maxvelcheck << "\t" << minvelcheck << "\t" << maxvmag << std::endl;

   maxvel = maxvmag;
   if(maxvelcheck > maxvmag) maxvel = maxvelcheck;
   if(minvelcheck > maxvel) maxvel = minvelcheck;
   maxvmag = maxvel;
   //velnd = 8.38e5; //CV
   //maxvel = maxvel/velnd; //CV

   MPI_Allreduce(&maxvmag,&maxvel,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);

   pname = s_part.name;
   numvels = 100;//(0.5)*s_part.pos.size()/s_msh.cmesh.size();
   delvel = 2.1*maxvel/numvels;
   delvelx = 2.1*maxvel/numvels;
   delvely = 2.1*maxvel/numvels;
   delvelz = 2.1*maxvel/numvels;


   for(i=0;i<numvels;i++) pdf.push_back(0);
   for(i=0;i<numvels;i++) velocity.push_back(delvel*(i+0.5));
   for(i=0;i<numvels;i++) pdfx.push_back(0);
   for(i=0;i<numvels;i++) pdfy.push_back(0);
   for(i=0;i<numvels;i++) pdfz.push_back(0);
   for(i=0;i<numvels;i++) velocityx.push_back(delvelx*((i-numvels/2)+0.5));
   for(i=0;i<numvels;i++) velocityy.push_back(delvely*((i-numvels/2)+0.5));
   for(i=0;i<numvels;i++) velocityz.push_back(delvelz*((i-numvels/2)+0.5));


   if(s_index<0)
   {
     fname << pname.c_str() << "totalvdf" <<  vdf_fnum[s_part.spflag] << ".dat";
     for(i=0;i<npg;i++)
     {
       vel = 0.0;
       for(j=0;j<s_msh.vecdims;j++) vel  = vel + (s_part.vel[s_msh.vecdims*i+j]*s_part.vel[s_msh.vecdims*i+j]);//(velnd*velnd);//CV
       vel = sqrt(vel);
       ind = (vel)/delvel;
       pdf[ind] +=1;
     }
     
     for(i=0;i<npg;i++)
     {
       vel  = s_part.vel[s_msh.vecdims*i];//velnd; //CV
       ind = (vel+0.5*numvels*delvelx)/delvelx;
       pdfx[ind] += 1;

       vel  = s_part.vel[s_msh.vecdims*i+1];//velnd; //CV
       ind = (vel+0.5*numvels*delvely)/delvely;
       pdfy[ind] += 1;

       vel  = s_part.vel[s_msh.vecdims*i+2];//velnd; //CV
       ind = (vel+0.5*numvels*delvelz)/delvelz;
       pdfz[ind] += 1;
     }

     filename = fname.str();
     std::ofstream wrtfile(filename.c_str());

     numpart = std::accumulate(pdf.begin(),pdf.end(),0);
     numpartx = std::accumulate(pdfx.begin(),pdfx.end(),0);
     numparty = std::accumulate(pdfy.begin(),pdfy.end(),0);
     numpartz = std::accumulate(pdfz.begin(),pdfz.end(),0);
   }
   else
   {
     fname << pname.c_str() << "cell" << s_index << "vdf" <<  vdf_fnum[s_part.spflag] << ".dat";
     for(i=0;i<npg;i++)
     {
       if(s_part.cell[i] == s_index)
       {
         vel = 0.0;
         for(j=0;j<s_msh.vecdims;j++) vel  = vel + s_part.vel[s_msh.vecdims*i+j]*s_part.vel[s_msh.vecdims*i+j];
         vel = sqrt(vel);
         ind = (vel)/delvel;
         pdf[ind] += 1;
       }
     }
   
     for(i=0;i<npg;i++)
     { 
       if(s_part.cell[i] == s_index)
       {
         vel  = s_part.vel[s_msh.vecdims*i];
         ind = (vel+0.5*numvels*delvelx)/delvelx;
         pdfx[ind] += 1;

         vel  = s_part.vel[s_msh.vecdims*i+1];
         ind = (vel+0.5*numvels*delvely)/delvelx;
         pdfy[ind] += 1;

         vel  = s_part.vel[s_msh.vecdims*i+2];
         ind = (vel+0.5*numvels*delvelz)/delvelz;
         pdfz[ind] += 1;
       }
     }

     filename = fname.str();

     numpart = std::accumulate(pdf.begin(),pdf.end(),0);
     numpartx = std::accumulate(pdfx.begin(),pdfx.end(),0);
     numparty = std::accumulate(pdfy.begin(),pdfy.end(),0);
     numpartz = std::accumulate(pdfz.begin(),pdfz.end(),0);
   }

   if(numprocs>1) //MPI
   {
     double numpartsum,numpartxsum, numpartysum, numpartzsum;
     std::vector<double> pdfsum;
     std::vector<double> pdfxsum;
     std::vector<double> pdfysum;
     std::vector<double> pdfzsum;
     for(i=0;i<numvels;i++) pdfsum.push_back(0);
     for(i=0;i<numvels;i++) pdfxsum.push_back(0);
     for(i=0;i<numvels;i++) pdfysum.push_back(0);
     for(i=0;i<numvels;i++) pdfzsum.push_back(0);
     
     MPI_Barrier(MPI_COMM_WORLD); //MPI  
     MPI_Reduce(&numpart,&numpartsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&numpartx,&numpartxsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&numparty,&numpartysum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&numpartz,&numpartzsum,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

     MPI_Reduce(&pdf.front(),&pdfsum.front(),pdf.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&pdfx.front(),&pdfxsum.front(),pdfx.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&pdfy.front(),&pdfysum.front(),pdfy.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&pdfz.front(),&pdfzsum.front(),pdfz.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

     if(procid==0)
     {
       std::ofstream wrtfile(filename.c_str());
 
       if(s_msh.wrtflag==0) //Debug Write
       {
         for(i=0;i<numvels;i++) wrtfile << velocity[i] << "\t" << pdfsum[i]/(numpartsum*delvel) << "\t" << velocityx[i] << "\t" << pdfxsum[i]/(numpartxsum*delvelx) << std::endl;
       }
       else if(s_msh.wrtflag==1) //Tecplot Write
       {
         wrtfile << "VARIABLES = \"V\", \"f\", \"Vx\", \"fx\", \"Vy\", \"fy\", \"Vz\", \"fz\"" << std::endl;
         wrtfile << "ZONE T=\"VDF\", I=" << numvels << ", DATAPACKING=POINT, SOLUTIONTIME=" << time << std::endl;
         for(i=0;i<numvels;i++) wrtfile << velocity[i] << "\t" << pdfsum[i]/(numpartsum*delvel) << "\t" << velocityx[i] << "\t" << pdfxsum[i]/(numpartxsum*delvelx) << "\t" << velocityy[i] << "\t" << pdfysum[i]/(numpartysum*delvely) << "\t" << velocityz[i] << "\t" << pdfzsum[i]/(numpartzsum*delvelz) << std::endl;
       }
       wrtfile.close();
     }
   }
   else
   {
     std::ofstream wrtfile(filename.c_str());

     if(s_msh.wrtflag==0) //Debug Write
     {
       for(i=0;i<numvels;i++) wrtfile << velocity[i] << "\t" << pdf[i]/(numpart*delvel) << "\t" << velocityx[i] << "\t" << pdfx[i]/(numpartx*delvelx) << std::endl;
     }
     else if(s_msh.wrtflag==1) //Tecplot Write
     {
       wrtfile << "VARIABLES = \"V\", \"f\", \"Vx\", \"fx\", \"Vy\", \"fy\", \"Vz\", \"fz\"" << std::endl;
       wrtfile << "ZONE T=\"VDF\", I=" << numvels << ", DATAPACKING=POINT, SOLUTIONTIME=" << time << std::endl;
       for(i=0;i<numvels;i++) wrtfile << velocity[i] << "\t" << pdf[i]/(numpart*delvel) << "\t" << velocityx[i] << "\t" << pdfx[i]/(numpartx*delvelx) << "\t" << velocityy[i] << "\t" << pdfy[i]/(numparty*delvely) << "\t" << velocityz[i] << "\t" << pdfz[i]/(numpartz*delvelz) << std::endl;
     }
     wrtfile.close();
   }

   vdf_fnum[s_part.spflag] = vdf_fnum[s_part.spflag]+1;

}


//..Find Phase-Space distribution..//

void writeOutput::findphasespace(particles s_part, mesh s_msh, int s_index, double time)
{
   int i,j,k;
   double maxvel = *std::max_element(s_part.vel.begin(),s_part.vel.end());
   double maxen = *std::max_element(s_part.en.begin(),s_part.en.end());
   std::vector<double> ps, velocity;
   std::vector<double> psx, velocityx;
   std::vector<double> psy, velocityy;
   std::vector<double> psz, velocityz;
   std::vector<double> numpartx,numparty,numpartz,numpart;
   double vel,delvel,delvelx,delvely,delvelz;
   double maxvelcheck1,maxvelcheck2,vmag,maxvmag;
   double maxvelcheck,minvelcheck;
   double velnd;  //TCH
   int numvels,ncells,npts,ind;
   std::stringstream fname;
   std::string filename,pname;
   long int np = s_part.pos.size()/s_msh.meshdims;
   long int npg  = np - s_part.gp;


   ncells = s_msh.cmesh.size();  //CHK

   int numprocs,procid; //MPI

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI


   maxen = *std::max_element(s_part.en.begin(),s_part.en.end());
   maxvelcheck = *std::max_element(s_part.vel.begin(),s_part.vel.end());
   minvelcheck = fabs(*std::min_element(s_part.vel.begin(),s_part.vel.end()));

   if(maxen==0) maxen = 1e-30;

   maxvmag = sqrt(maxen*2.0/s_part.mass);
   //std::cout << std::endl <<  maxen << "\t" << maxvelcheck << "\t" << minvelcheck << "\t" << maxvmag << std::endl;

   maxvel = maxvmag;
   if(maxvelcheck > maxvmag) maxvel = maxvelcheck;
   if(minvelcheck > maxvel) maxvel = minvelcheck;
   maxvmag = maxvel;
   //velnd = 8.38e5; //CV
   //maxvel = maxvel/velnd; //CV

   MPI_Allreduce(&maxvmag,&maxvel,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);


   pname = s_part.name;
   numvels = 100;
   delvel = 2.1*maxvel/numvels;
   delvelx = 2.1*maxvel/numvels;
   delvely = 2.1*maxvel/numvels;
   delvelz = 2.1*maxvel/numvels;

   npts = ncells*numvels;

   for(i=0;i<npts;i++) ps.push_back(0);
   for(i=0;i<numvels;i++) velocity.push_back(delvel*(i+0.5));
   //for(i=0;i<numvels;i++) velocity.push_back(delvel*((i-numvels/2)+0.5));
   for(i=0;i<npts;i++) psx.push_back(0);
   for(i=0;i<npts;i++) psy.push_back(0);
   for(i=0;i<npts;i++) psz.push_back(0);
   for(i=0;i<numvels;i++) velocityx.push_back(delvelx*((i-numvels/2)+0.5));
   for(i=0;i<numvels;i++) velocityy.push_back(delvely*((i-numvels/2)+0.5));
   for(i=0;i<numvels;i++) velocityz.push_back(delvelz*((i-numvels/2)+0.5));
   for(i=0;i<ncells;i++) numpart.push_back(1.0);
   for(i=0;i<ncells;i++) numpartx.push_back(1.0);
   for(i=0;i<ncells;i++) numparty.push_back(1.0);
   for(i=0;i<ncells;i++) numpartz.push_back(1.0);


   fname << pname.c_str() << "phasespace" <<  ps_fnum[s_part.spflag] << ".dat";

   for(i=0;i<npg;i++)
   {
     vel = 0.0;
     for(j=0;j<s_msh.vecdims;j++) vel  = vel + (s_part.vel[s_msh.vecdims*i+j]*s_part.vel[s_msh.vecdims*i+j]);//(velnd*velnd);//CV
     vel = sqrt(vel);
     ind = (vel)/delvel;
     ps[s_part.cell[i]*numvels+ind] +=1;
     
   }
     
   for(i=0;i<npg;i++)
   {
     vel  = s_part.vel[s_msh.vecdims*i];//velnd; //CV
     ind = (vel+0.5*numvels*delvelx)/delvelx;
     psx[s_part.cell[i]*numvels+ind] += 1;

     vel  = s_part.vel[s_msh.vecdims*i+1];//velnd; //CV
     ind = (vel+0.5*numvels*delvely)/delvely;
     psy[s_part.cell[i]*numvels+ind] += 1;

     vel  = s_part.vel[s_msh.vecdims*i+2];//velnd; //CV
     ind = (vel+0.5*numvels*delvelz)/delvelz;
     psz[s_part.cell[i]*numvels+ind] += 1;

   }

   filename = fname.str();
   std::ofstream wrtfile(filename.c_str());

   for(i=1;i<(ncells-1);i++)
   {
     numpart[i] = std::accumulate(ps.begin()+i*numvels,ps.begin()+(i+1)*numvels,0);
     numpartx[i] = std::accumulate(psx.begin()+i*numvels,psx.begin()+(i+1)*numvels,0);
     numparty[i] = std::accumulate(psy.begin()+i*numvels,psy.begin()+(i+1)*numvels,0);
     numpartz[i] = std::accumulate(psz.begin()+i*numvels,psz.begin()+(i+1)*numvels,0);
   }   

   if(numprocs>1) //MPI
   {
     std::vector<double> numpartxsum,numpartysum,numpartzsum,numpartsum;
     std::vector<double> pssum;
     std::vector<double> psxsum;
     std::vector<double> psysum;
     std::vector<double> pszsum;
     for(i=0;i<npts;i++) pssum.push_back(0);
     for(i=0;i<npts;i++) psxsum.push_back(0);
     for(i=0;i<npts;i++) psysum.push_back(0);
     for(i=0;i<npts;i++) pszsum.push_back(0);
     for(i=0;i<ncells;i++) numpartsum.push_back(1.0);
     for(i=0;i<ncells;i++) numpartxsum.push_back(1.0);
     for(i=0;i<ncells;i++) numpartysum.push_back(1.0);
     for(i=0;i<ncells;i++) numpartzsum.push_back(1.0);
     
     MPI_Barrier(MPI_COMM_WORLD); //MPI  
     MPI_Reduce(&numpart.front(),&numpartsum.front(),ncells,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&numpartx.front(),&numpartxsum.front(),ncells,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&numparty.front(),&numpartysum.front(),ncells,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&numpartz.front(),&numpartzsum.front(),ncells,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

     MPI_Reduce(&ps.front(),&pssum.front(),ps.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&psx.front(),&psxsum.front(),psx.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&psy.front(),&psysum.front(),psy.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
     MPI_Reduce(&psz.front(),&pszsum.front(),psz.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
   
     if(procid==0)
     {
       if(s_msh.wrtflag==0) //Debug Write
       {
         for(j=0;j<ncells;j++) 
         {
           for(i=0;i<numvels;i++)
           {
             wrtfile << s_msh.cmesh[j] << "\t" << velocity[i] << "\t" << pssum[j*numvels+i]/(numpartsum[j]*delvel) << "\t" << velocityx[i] << "\t" << psxsum[j*numvels+i]/(numpartxsum[j]*delvelx) << "\t" << psysum[j*numvels+i]/(numpartysum[j]*delvely) << "\t" << pszsum[j*numvels+i]/(numpartzsum[j]*delvelz) << std::endl;
           }
         }
       }
       else if(s_msh.wrtflag==1) //Tecplot Write
       {
         wrtfile << "VARIABLES = \"X\", \"V\", \"f\", \"Vi\", \"fx\", \"fy\", \"fz\"" << std::endl;
         wrtfile << "ZONE T=\"PHASESPACE\", I=" << ncells << ",J=" << numvels << ", DATAPACKING=POINT, SOLUTIONTIME=" << time << std::endl;
     
         for(i=0;i<numvels;i++)
         {
           for(j=0;j<ncells;j++) 
           {
             wrtfile << s_msh.cmesh[j] << "\t" << velocity[i] << "\t" << pssum[j*numvels+i]/(numpartsum[j]*delvel) << "\t" << velocityx[i] << "\t" << psxsum[j*numvels+i]/(numpartxsum[j]*delvelx) << "\t" << psysum[j*numvels+i]/(numpartysum[j]*delvely) << "\t" << pszsum[j*numvels+i]/(numpartzsum[j]*delvelz) << std::endl;
           }
         }
       }
     }
   }
   else
   {
     if(s_msh.wrtflag==0) //Debug Write
     {
       for(j=0;j<ncells;j++) 
       {
         for(i=0;i<numvels;i++)
         {
           wrtfile << s_msh.cmesh[j] << "\t" << velocity[i] << "\t" << ps[j*numvels+i]/(numpart[j]*delvel) << "\t" << psx[j*numvels+i]/(numpartx[j]*delvelx) << "\t" << psy[j*numvels+i]/(numparty[j]*delvely) << "\t" << psz[j*numvels+i]/(numpartz[j]*delvelz) << std::endl;
         }
       }
     }
     else if(s_msh.wrtflag==1) //Tecplot Write
     {
       wrtfile << "VARIABLES = \"X\", \"V\", \"f\", \"Vi\", \"fx\", \"fy\", \"fz\"" << std::endl;
       wrtfile << "ZONE T=\"PHASESPACE\", I=" << ncells << ",J=" << numvels << ", DATAPACKING=POINT, SOLUTIONTIME=" << time << std::endl;
     
       for(i=0;i<numvels;i++)
       {
         for(j=0;j<ncells;j++) 
         {
           wrtfile << s_msh.cmesh[j] << "\t" << velocity[i] << "\t" << ps[j*numvels+i]/(numpart[j]*delvel) << "\t" << velocityx[i] << "\t"<< psx[j*numvels+i]/(numpartx[j]*delvelx) << "\t" << psy[j*numvels+i]/(numparty[j]*delvely) << "\t" << psz[j*numvels+i]/(numpartz[j]*delvelz) << std::endl;
         }
       }
     }
   }

   ps_fnum[s_part.spflag] = ps_fnum[s_part.spflag]+1;

   wrtfile.close();
}

//...Write information used to restart simulation...//

void writeOutput::writeRestart(const std::vector<boundvars> &w_boundvars, double w_time, int w_iter, int w_nsp)
{
   int i,j;
   std::stringstream fname;
   std::string filename;

   std::string name;

   std::cout << "\n\tWriting Restart File...." ;

   name = "restart";

   fname <<  name.c_str() << rst_fnum << ".out";

   filename = fname.str();

   std::ofstream wrtfile(filename.c_str());
  
   wrtfile << w_iter << "\t" << w_time << std::endl << std::endl;


   for(i=0;i<w_nsp;i++)
   {
     wrtfile << part_fnum[i]-1 << "\t";
     wrtfile << cont_fnum[i]-1 << "\t";
     wrtfile << field_fnum[i]-1 << "\t";
     wrtfile << vdf_fnum[i]-1 << "\t";
     wrtfile << ps_fnum[i]-1 << std::endl;
   } 

   wrtfile << totcont_fnum-1 << std::endl;
   wrtfile << rst_fnum << std::endl;
   rst_fnum = rst_fnum + 1;
 
   
   for(i=0;i<w_boundvars[0].nbound;i++)
   {
     for(j=0;j<w_nsp;j++)
     {
       wrtfile << w_boundvars[i].partcount[j] << "\t";
     }
     wrtfile << std::endl;
     wrtfile << w_boundvars[i].sigma << "\n";
   }

   wrtfile.close();   
}

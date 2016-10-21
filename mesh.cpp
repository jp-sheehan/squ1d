#include "mesh.h"

//..Generate Euler mesh (CFD class)..//

void mesh::generateEulermesh(double m_lscale)
{
   int i;
   std::vector<double> offset;

   for(i=0;i<meshdims; i++) offset.push_back(0.0); 

   // Offset not fully coded in but left here....

   meshPoint(pmesh,offset,m_lscale); 
   meshCenter(cmesh,offset,m_lscale); 
}

//..Generate PIC mesh..//

void mesh::generatePICmesh(double m_lscale)
{
   int i;
   std::vector<double> offset;

   for(i=0;i<meshdims; i++) offset.push_back(0.0); 

   // Offset not fully coded in but left here....

   meshPoint(pmesh,offset,m_lscale); 
   meshCenter(cmesh,offset,m_lscale); 

   //..Find neighbors of cells and points..//  
   meshPointNeighbors(pofpneigh,cofpneigh);
   //meshPointNeighbors();
   meshCenterNeighbors(cofcneigh,pofcneigh);

   // Not setup for 2D

   deltax = (meshlength[0])/(numpoints[0]-1);
   hdeltax  = 0.5*deltax;// half delta x?
   //std::cout << "\nDelta:\t" << deltax << "\t" << hdeltax << std::endl;

}


//..Generate mesh for points(cell corners)..//

void mesh::meshPoint(std::vector<double> & m_mesh, std::vector<double> m_offset, double m_lscale)
{
   int i,j,k;
   int m_totalcells = 1;
   double temp1,temp2,temp3;
   double tol = 0.0;

   std::cout << "\n\nGenerating Corner Mesh....... \n";

   for (i=0; i<meshdims; i++) m_totalcells = m_totalcells*numpoints[i];

   std::cout << "Total Points:\t"  << m_totalcells << std::endl;

   if(meshdims == 1)
   {
      for(i=0; i<numpoints[0] ;i++)
      {
         temp1 = i*m_lscale*meshlength[0]/(numpoints[0]-1)+m_offset[0]+m_lscale*meshstart[0];
         m_mesh.push_back(temp1+tol); // when called from generatePICmesh, m_mesh = pmesh
      }
   }
   /*
   else if (meshdims == 2)
   {
      for(j=0; j<numpoints[1]; j++)
      {
         temp2 = j*m_lscale*meshlength[1]/(numpoints[1]-1) + m_offset[1]+m_lscale*meshstart[1];

         for(i=0; i<numpoints[0]; i++) 
         {
            temp1 = i*m_lscale*meshlength[0]/(numpoints[0]-1) + m_offset[0]+m_lscale*meshstart[0];
            m_mesh.push_back(temp1+tol);
            m_mesh.push_back(temp2+tol);
         }
      }
   }
   else if (meshdims == 3)
   {
      for(k=0; k<numpoints[2]; k++)
      {
         temp3 = k*m_lscale*meshlength[2]/(numpoints[2]-1) + m_offset[2]+m_lscale*meshstart[2];

         for(j=0; j<numpoints[1]; j++) 
         {
            temp2 = j*m_lscale*meshlength[1]/(numpoints[1]-1) + m_offset[1] + m_lscale*meshstart[1];

            for(i=0; i<numpoints[0]; i++)
            {  
               temp1 = i*m_lscale*meshlength[0]/(numpoints[0]-1) + m_offset[0] + m_lscale*meshstart[0];
               m_mesh.push_back(temp1+tol);
               m_mesh.push_back(temp2+tol);
               m_mesh.push_back(temp3+tol);
            }
         } 
      }
   } 
   */
}

//..Generate mesh for cell centers..//

void mesh::meshCenter(std::vector<double> & m_mesh, std::vector<double> m_offset, double m_lscale)
{
   int i,j,k;
   int m_totalcells = 1;
   double temp1,temp2,temp3;
   double tol = 0.0;

   std::cout << "\nGenerating Center Mesh....... \n";

   for (i=0; i<meshdims; i++) m_totalcells = m_totalcells*(numpoints[i]+2*numghost-1);

   std::cout << "Total Cells:\t"  << m_totalcells << std::endl;

   if(meshdims == 1)
   {
      for(i=0; i<(numpoints[0]+2*numghost-1) ;i++)
      {
         temp1 = (i-numghost+0.5)*m_lscale*meshlength[0]/(numpoints[0]-1) + m_offset[0] + m_lscale*meshstart[0];
         m_mesh.push_back(temp1+tol); // when called from generatePICmesh, m_mesh = cmesh
      }
   }
   /*
   else if (meshdims == 2)
   {
      for(j=0; j<(numpoints[1]+2*numghost-1); j++)
      {
         temp2 = (j-numghost+0.5)*m_lscale*meshlength[1]/(numpoints[1]-1) + m_lscale*meshstart[1];

         for(i=0; i<(numpoints[0]+2*numghost-1); i++) 
         {
            temp1 = (i-numghost+0.5)*m_lscale*meshlength[0]/(numpoints[0]-1) + m_lscale*meshstart[0];
            m_mesh.push_back(temp1+tol);
            m_mesh.push_back(temp2+tol);
         }
      }
   }
   else if (meshdims == 3)
   {
      for(k=0; k<(numpoints[2]+2*numghost-1); k++)
      {
         temp1 = (0.5)*(2*k-numghost-1)*m_lscale*meshlength[2]/(numpoints[2]-1) + m_offset[2] + m_lscale*meshstart[2];

         for(j=0; j<(numpoints[1]+2*numghost-1); j++) 
         {
            temp2 = (0.5)*(2*j-numghost-1)*m_lscale*meshlength[1]/(numpoints[1]-1) + m_offset[1] + m_lscale*meshstart[1];

            for(i=0; i<(numpoints[0]+2*numghost-1); i++)
            {  
               temp3 = (0.5)*(2*i-numghost-1)*m_lscale*meshlength[0]/(numpoints[0]-1) + m_offset[0] + m_lscale*meshstart[0];
               m_mesh.push_back(temp1+tol);
               m_mesh.push_back(temp2+tol);
               m_mesh.push_back(temp3+tol);
            }
         } 
      }
   } 
   */
}

//..Generate mesh for cell faces..//

void mesh::meshFace(double m_lscale)
{
   int i,j,k;
   int m_totalcells = 1;
   double temp1,temp2,temp3;

   std::cout << "\nGenerating Wall Mesh....... \n";

   for (i=0; i<meshdims; i++) m_totalcells = m_totalcells*(numpoints[i]-1);

   std::cout << "Total Cells:\t"  << m_totalcells << std::endl;

   if(meshdims == 1)
   {
      for(i=0; i<(numpoints[0]) ;i++)
      {
         temp1 = i*m_lscale*meshlength[0]/(numpoints[0]);
         fmesh.push_back(temp1);
      }
   }
   /*
   else if (meshdims == 2)
   {
      for(j=0; j<(numpoints[1]); j++)
      {
         temp2 = (j)*m_lscale*meshlength[1]/(numpoints[1]-1);

         for(i=0; i<(numpoints[0]-1); i++) 
         {
            temp1 = (0.5)*(2*i+1)*m_lscale*meshlength[0]/(numpoints[0]-1);
            fmesh.push_back(temp1);
            fmesh.push_back(temp2);
         }
      }

      for(j=0; j<(numpoints[1]-1); j++)
      {
         temp2 = (0.5)*(2*j+1)*m_lscale*meshlength[1]/(numpoints[1]-1);

         for(i=0; i<(numpoints[0]); i++) 
         {
            temp1 = i*m_lscale*meshlength[0]/(numpoints[0]-1);
            fmesh.push_back(temp1);
            fmesh.push_back(temp2);
         }
      }
   }
   else if (meshdims == 3)
   {
      for(k=0; k<(numpoints[2]-1); k++)
      {
         temp1 = (0.5)*(2*k+1)*m_lscale*meshlength[2]/(numpoints[2]-1);

         for(j=0; j<(numpoints[1]-1); j++) 
         {
            temp2 = (0.5)*(2*j+1)*m_lscale*meshlength[1]/(numpoints[1]-1);

            for(i=0; i<(numpoints[0]-1); i++)
            {  
               temp3 = (0.5)*(2*i+1)*m_lscale*meshlength[0]/(numpoints[0]-1);
               fmesh.push_back(temp1);
               fmesh.push_back(temp2);
               fmesh.push_back(temp3);
            }
         } 
      }
   } 
   */
}

//..Find size of mesh..//

int mesh::mesh_size(std::vector<double> m_mesh)
{
   return m_mesh.size()/meshdims;
}


//..Find internal cell center neighbor of cell center..//

int mesh::incneighofc(int m_cindex) const   //Find Internal Neighbor
{
   int i,inneigh;
   std::vector<int> neighbors;

   for(i=0;i<meshdims*meshdims;i++) neighbors.push_back(0);
   neighbors = cneighofc(m_cindex);

   if(neighbors[0] == -1) inneigh = neighbors[1];
   else if(neighbors[1] == -1) inneigh = neighbors[0];
   else if(neighbors[2] == -1) inneigh = neighbors[3];
   else if(neighbors[3] == -1) inneigh = neighbors[2];
   else inneigh = -1;

   return inneigh;
}

//..Find cell center neighbors of point..//

std::vector<int> mesh::cneighofp(int m_cindex) const
{
   std::vector<int> neighbors;
   //int counter = m_cindex/(numpoints[0]);

   if(meshdims ==1 )
   {
      neighbors.push_back(m_cindex);
      neighbors.push_back(m_cindex+1);
   }
   /*
   else if(meshdims == 2)
   {
      neighbors.push_back(m_cindex+counter);
      neighbors.push_back(m_cindex+counter+1);
      neighbors.push_back(m_cindex+numpoints[0]-1+2*numghost+counter);
      neighbors.push_back(m_cindex+numpoints[0]+2*numghost+counter);
   }
   */

   return neighbors;
}

//..Find point neighbors of cell center..//

std::vector<int> mesh::pneighofc(int m_cindex) const
{
   int i;
   std::vector<int> neighbors;
   //int counter = m_cindex/(numpoints[0]-1+2*numghost);

   if(meshdims==1)
   {
      neighbors.push_back(m_cindex-1);
      neighbors.push_back(m_cindex);
      if(m_cindex == (numpoints[0]-1+2*numghost))     neighbors[1] = -1; 
   }
   /*
   else if(meshdims==2)
   {
      neighbors.push_back(m_cindex-(numpoints[0]-1+2*numghost)-counter);
      neighbors.push_back(m_cindex-(numpoints[0]-1+2*numghost)-counter+1);
      neighbors.push_back(m_cindex-counter-1);
      neighbors.push_back(m_cindex-counter);

      if(m_cindex == (counter*(numpoints[0]-1+2*numghost)))
      {
         neighbors[0] = -1; 
         neighbors[2] = -1;
      }

      if(m_cindex == ((counter+1)*(numpoints[0]-1+2*numghost)-1)) 
      {
         neighbors[1] = -1; 
         neighbors[3] = -1;
      }

      if(m_cindex < (numpoints[0]-1+2*numghost))
      {
         neighbors[0] = -1;
         neighbors[1] = -1;
      }

      if(m_cindex >= (numpoints[0]-1+2*numghost)*(numpoints[1]-2+2*numghost))
      {
         neighbors[2] = -1;
         neighbors[3] = -1;
      }   
   }
   */

   return neighbors;
}

//..Find cell center neighbors of cell center..//

std::vector<int> mesh::cneighofc(int m_cindex) const
{
   std::vector<int> neighbors;

   if(meshdims==1)
   {
      neighbors.push_back(m_cindex-1);
      neighbors.push_back(m_cindex+1);

      if( (m_cindex+1) == (numpoints[0]-1+2*numghost)) neighbors[1] =-1;
   }
   /*
   else if(meshdims==2)
   {
      neighbors.push_back(m_cindex-1);
      neighbors.push_back(m_cindex+1);
      neighbors.push_back(m_cindex-(numpoints[0]-1+2*numghost));
      neighbors.push_back(m_cindex+(numpoints[0]-1+2*numghost));

      if(m_cindex != 0 && m_cindex % (numpoints[0]-1+2*numghost) == 0)  neighbors[0] = -1;
      if((m_cindex+1) % (numpoints[0]-1+2*numghost) == 0) neighbors[1] =-1;
      if(m_cindex < (numpoints[0])-1+2*numghost)  neighbors[2] = -1;
      if(m_cindex >= ((numpoints[0]-1+2*numghost)*(numpoints[1]-2+2*numghost))) neighbors[3] = -1;
   }
   */

   return neighbors;
}


//..Find point neighbors of point..//

std::vector<int> mesh::pneighofp(int m_cindex) const
{
   std::vector<int> neighbors;

   if(meshdims==1)
   {
      neighbors.push_back(m_cindex-1);
      neighbors.push_back(m_cindex+1);
      if((m_cindex+1) == numpoints[0]) neighbors[1] =-1;

      if(perflag==1 && neighbors[0]==-1 && philoc==1) neighbors[0] = numpoints[0]-2; 
      else if(perflag==1 && neighbors[1]==-1 && philoc==1) neighbors[1] = 1; 
   }
   /*
   else if(meshdims==2)   //FIX,GS
   { 
      neighbors.push_back(m_cindex-1);
      neighbors.push_back(m_cindex+1);
      neighbors.push_back(m_cindex-(numpoints[0]));
      neighbors.push_back(m_cindex+(numpoints[0]));

      if(m_cindex < numpoints[0])  neighbors[2] = -1;
      if(m_cindex >= ((numpoints[0])*(numpoints[1]-1))) neighbors[3] = -1;
      if(m_cindex != 0 && m_cindex % (numpoints[0]) == 0)  neighbors[0] = -1;
      if((m_cindex+1) % (numpoints[0]) == 0) neighbors[1] =-1;
   }  
   */

   return neighbors;
}

//......Unfinished....Initially for mesh shifting neighbors.....///

std::vector<int> mesh::pneighofp(std::vector<double> m_mesh_shift, int m_neigh,int m_meshpt,int m_cindex) const
{
   int i;
   std::vector<double>  m_correct;
   std::vector<int> neighbors;

   for(i=0;i<meshdims;i++)  
   {
      if(m_mesh_shift[m_meshpt*meshdims+i] == m_mesh_shift[m_neigh*meshdims+i])  m_correct.push_back(0);
      else if(m_mesh_shift[m_meshpt*meshdims+i] < m_mesh_shift[m_neigh*meshdims+i])  m_correct.push_back(-1);
      else if(m_mesh_shift[m_meshpt*meshdims+i] > m_mesh_shift[m_neigh*meshdims+i])  m_correct.push_back(0);
   }  

   neighbors.push_back(m_cindex-(numpoints[0]));
   neighbors.push_back(m_cindex-1);
   neighbors.push_back(m_cindex+1);
   neighbors.push_back(m_cindex+(numpoints[0]));

   if((m_cindex) % (numpoints[0]) == 0)  neighbors[1] = -1;
   if(m_cindex >= ((numpoints[0])*(numpoints[1]-1))) neighbors[3] = -1;
   if((m_cindex+1) % (numpoints[0]) == 0) neighbors[2] =-1;
   if(m_cindex < numpoints[0])  neighbors[0] = -1;

   return neighbors;
}

//..Find face neighbors of cell center..//

std::vector<int> mesh::fneighofc(int m_cindex) const
{
   std::vector<int> neighbors;
   int counter = m_cindex/(numpoints[0]-1);

   neighbors.push_back(m_cindex);
   neighbors.push_back(m_cindex+numpoints[0]);
   neighbors.push_back(m_cindex+numpoints[0]*(numpoints[1]+1)+counter);
   neighbors.push_back(m_cindex+numpoints[0]*(numpoints[1]+1)+1+counter);

   return neighbors;
}

//..Find cell center neighbors of face..//

std::vector<int> mesh::cneighoff(int m_cindex) const
{
   std::vector<int> neighbors;

   if(m_cindex < numpoints[0])
   {  
      neighbors.push_back(-1);
      neighbors.push_back(m_cindex);
   }
   else if(m_cindex >= numpoints[0] && m_cindex < (numpoints[0]*(numpoints[1])))
   {
      neighbors.push_back(m_cindex-(numpoints[0]));
      neighbors.push_back(m_cindex);
   }
   else if(m_cindex >= (numpoints[0]*(numpoints[1])) && m_cindex < (numpoints[0]*(numpoints[1]+1)))
   {
      neighbors.push_back(m_cindex-(numpoints[0]));
      neighbors.push_back(-1);
   }
   else if(m_cindex >= (numpoints[0]*(numpoints[1]+1)))
   {
      if((m_cindex % numpoints[0])==0)  neighbors.push_back(-1);
      else  neighbors.push_back(m_cindex-numpoints[0]*(numpoints[1]+1));

      if(((m_cindex+1) % numpoints[0])==0)  neighbors.push_back(-1);
      else  neighbors.push_back(m_cindex-numpoints[0]*(numpoints[1]+1)-1);
   }

   return neighbors;
}


//..Find point neighbors of face..//

std::vector<int> mesh::pneighoff(int m_cindex) const
{
   std::vector<int> neighbors;
   int counter = m_cindex/(numpoints[0]-1);

   if(m_cindex < (numpoints[0]*(numpoints[1]+1)))
   { 
      neighbors.push_back(m_cindex+counter);
      neighbors.push_back(m_cindex+counter+1);
   }
   else 
   {  
      neighbors.push_back((m_cindex-numpoints[0]*(numpoints[1]+1)));
      neighbors.push_back((m_cindex-numpoints[0]*(numpoints[1]+1))+(numpoints[0]+1));
   }

   return neighbors;
}

//..Find nearest point (corner point) to particle..//

int mesh::nearp(particles const&  m_part, int m_pindex) const
{
   int i,j,nearestpoint;
   double tempdist1, tempdist2;
   std::vector<double> temp1,temp2;

   nearestpoint = 8675309;
   for(i=0;i<meshdims;i++)
   {  
      temp1.push_back(m_part.pos[m_pindex*meshdims+i]);
      temp2.push_back(0.0);
   }

   for(i=0;i<(pmesh.size()/meshdims); i++)
   {
      for(j=0;j<meshdims;j++) temp2[j] = pmesh[i*meshdims+j];

      tempdist1 = dist(temp1, temp2);
      //std::cout << tempdist1 << "\t" << tempdist1-hdeltax << "\t";
      if((tempdist1-hdeltax) < 1.0e-14)
      {
         //std::cout << ".....FOUND!.......";
         nearestpoint = i;
         break;
      }
   }
   return nearestpoint; 
}



//..Find nearest point (corner point) to particle 1D..//

int mesh::nearp1D(particles const&  m_part, int m_pindex) const
{
   int i,j,nearestpoint;
   double tempdist1, tempdist2;
   double temp1,temp2;

   temp1 = m_part.pos[m_pindex];
   nearestpoint = 8675309;

   for(i=0;i<pmesh.size(); i++)
   {
      tempdist1 = dist1D(temp1, pmesh[i]);
      //std::cout << tempdist1 << "\t" << hdeltax << "\t" << tempdist1-hdeltax << "\t";

      if((tempdist1-hdeltax) < 1.0e-14)
      {
         //std::cout << ".....FOUND!.......";
         nearestpoint = i;
         break;
      }
   }
   return nearestpoint; 
}


//..Find nearest point (corner point) to particle 1D..//   //DC

int mesh::nearp1Dnew(particles const&  m_part, int m_pindex) const
{
   int i,j,nearestpoint;
   double temp;

   temp = (m_part.pos[m_pindex*meshdims]+hdeltax-meshstart[0])/deltax;    //Forced rounding using truncation
   nearestpoint = temp;

   return nearestpoint; 
}



//..Find nearest point (corner point) to data..//

int mesh::nearp(std::vector<double> temp1) const
{
   int i,j;
   int nearestpoint = -1;
   double tempdist1, tempdist2;
   std::vector<double> temp2;

   for(i=0;i<meshdims;i++)  temp2.push_back(0.0);

   for(i=0;i<(pmesh.size()/meshdims); i++)
   {
      for(j=0;j<meshdims;j++)  temp2[j] = pmesh[i*meshdims+j];

      tempdist1 = dist(temp1, temp2);
      if((tempdist1-hdeltax) < 1.0e-14)
      {
         nearestpoint = i;
         break;
      }
   }
   return nearestpoint; 
}

//..Find nearest cell center to particle..//

int mesh::nearc(particles const& m_part, int m_pindex) const
{
   int i,j,nearestcell;
   double tempdist1, tempdist2;
   std::vector<double> temp1,temp2;
   nearestcell = 8675309;

   for(i=0;i<meshdims;i++)
   {  
      temp1.push_back(m_part.pos[m_pindex*meshdims+i]);
      temp2.push_back(0.0);
   }

   for(i=0;i<(cmesh.size()/meshdims); i++)
   {
      for(j=0;j<meshdims;j++) temp2[j] = cmesh[i*meshdims+j];
      tempdist1 = dist(temp1, temp2);

      if((tempdist1-hdeltax) < 1.0e-14)
      {
         nearestcell = i;
         break;
      }
   }
   return nearestcell; 
}

//..Find nearest cell center to particle (1D)..//

int mesh::nearc1D(particles const& m_part, int m_pindex) const
{
   int i,j,nearestcell;
   double tempdist1, tempdist2;
   double temp1,temp2;

   temp1 = (m_part.pos[m_pindex]);
   nearestcell = 8675309;

   for(i=0;i<cmesh.size(); i++)
   {
      tempdist1 = dist1D(temp1,cmesh[i]);
      //std::cout << tempdist1 << "\t" << hdeltax << "\t" << tempdist1-hdeltax << "\t";
      if((tempdist1-hdeltax) < 1.0e-14)
      {
         //std::cout << ".....FOUND!.......";
         nearestcell = i;
         break;
      }
   }
   return nearestcell; 
}


//..Find nearest cell center to particle (1D)..//

int mesh::nearc1Dnew(particles const& m_part, int m_pindex) const         //DC
{
   int i,j,nearestcell;
   double temp;

   temp = (m_part.pos[m_pindex*meshdims]+deltax-meshstart[0])/deltax;    //Forced rounding using truncation
   nearestcell = temp;

   return nearestcell; 
}




//..Find nearest cell center to data..//

int mesh::nearc(std::vector<double> temp1) const
{
   int i,j;
   int nearestcell = -1;
   double tempdist1;
   std::vector<double> temp2;

   for(i=0;i<meshdims;i++)  temp2.push_back(0.0);

   for(i=0;i<(cmesh.size()/meshdims); i++)
   {
      for(j=0;j<meshdims;j++) temp2[j] = cmesh[i*meshdims+j];

      tempdist1 = dist(temp1, temp2);
      if((tempdist1-hdeltax) < 1.0e-14)
      {
         nearestcell = i;
         break;
      }
   }
   return nearestcell; 
}


//..Find nearest face point to particle..//

int mesh::nearf(particles const& m_part, int m_pindex) const
{
   int i,j;
   int nearestcell = -1;
   double tempdist1, tempdist2;
   std::vector<double> temp1,temp2;

   //tempdist1 = 1E10;
   tempdist2 = 1E10;

   for(i=0;i<meshdims;i++)
   {  
      temp1.push_back(m_part.pos[m_pindex*meshdims+i]);
      temp2.push_back(0);
   }

   for(i=0;i<(fmesh.size()/meshdims); i++)
   {
      for(j=0;j<meshdims;j++)
      {
         temp2[j] = fmesh[i*meshdims+j];
      }

      tempdist1 = dist(temp1, temp2);

      if(tempdist1<tempdist2)
      {
         tempdist2 = tempdist1;
         nearestcell = i;
      }
   }
   return nearestcell; 
}

//..Distance between two points..//

double mesh::dist(std::vector<double> m_pt1, std::vector<double> m_pt2) const
{
   int i;
   double temp = 0.0;

   for(i=0; i<meshdims; i++) temp = temp + (m_pt1[i] - m_pt2[i])*(m_pt1[i] - m_pt2[i]);

   return sqrt(temp);
}

//..Distance between two points (1D)..//

double mesh::dist1D(double m_pt1,double m_pt2) const
{
   return fabs(m_pt1 - m_pt2);
}

//..Cell volume..//

/*double mesh::cellvolume(int m_cindex) const 
  {
  int i;
  std::vector<int> neighbors;
  double volume; 
  volume = 1.0;    

  for(i=0;i<2*meshdims; i++) neighbors.push_back(0);

//neighbors = pneighofc(m_cindex);
for(i=0;i<2*meshdims;i++) neighbors[i] = pofcneigh[2*m_cindex+i];   

if (meshdims==1) volume = (pmesh[neighbors[0]]-pmesh[neighbors[1]])*area[m_cindex];
else if(meshdims==2) volume = squarearea(pmesh[meshdims*neighbors[0]],pmesh[meshdims*neighbors[0]+1],pmesh[meshdims*neighbors[3]],pmesh[meshdims*neighbors[3]+1]);

return fabs(volume);
}*/


//..Cell volume..//      //EFF

double mesh::cellvolume(int m_cindex) const 
{
   return deltax*carea[m_cindex];
}


//..Volume around point(corner point)..//

/*double mesh::pointvolume(int m_pindex) const 
  {
  int i;
  std::vector<int> neighbors;
  double volume; 
  volume = 1.0;    


  for(i=0;i<2*meshdims; i++) neighbors.push_back(0);

//neighbors = cneighofp(m_pindex);
for(i=0;i<2*meshdims;i++) neighbors[i] = cofpneigh[2*m_pindex+i];   

if (meshdims==1) volume = (cmesh[neighbors[0]]-cmesh[neighbors[1]])*0.5*(area[neighbors[0]]+area[neighbors[1]]);
else if (meshdims == 2) volume = squarearea(cmesh[meshdims*neighbors[0]],cmesh[meshdims*neighbors[0]+1],cmesh[meshdims*neighbors[3]],cmesh[meshdims*neighbors[3]+1]);

return fabs(volume);
}*/


//..Volume around point(corner point)..//   //EFF

double mesh::pointvolume(int m_pindex) const 
{
   return deltax*(parea[m_pindex]);
}


//..Area of square..//

double mesh::squarearea(double x1,double y1,double x2,double y2) const
{
   double area; 

   area = (x2-x1)*(y2-y1);

   return fabs(area);
}

//..Area of face..//

double mesh::facearea(int m_cindex) const 
{
   int i;
   std::vector<int> neighbors;
   double area;

   for(i=0;i<meshdims; i++) neighbors.push_back(0);

   neighbors = pneighoff(m_cindex);

   area  = pmesh[meshdims*neighbors[1]]-pmesh[meshdims*neighbors[0]];

   return area;
}

//..Area weighting for point based on particle location..//

std::vector<double> mesh::pt_areaweight(double px,double py, int pindex) const
{
   double dx,dy,temp;
   std::vector<double> areaweight;

   dx = (0.5)*meshlength[0]/(numpoints[0]-1);
   dy = (0.5)*meshlength[1]/(numpoints[1]-1);

   areaweight.push_back(squarearea((px-dx),(py-dy),pmesh[meshdims*pindex],pmesh[meshdims*pindex+1]));
   areaweight.push_back(squarearea((px+dx),(py-dy),pmesh[meshdims*pindex],pmesh[meshdims*pindex+1]));
   areaweight.push_back(squarearea((px-dx),(py+dy),pmesh[meshdims*pindex],pmesh[meshdims*pindex+1]));
   areaweight.push_back(squarearea((px+dx),(py+dy),pmesh[meshdims*pindex],pmesh[meshdims*pindex+1]));

   return areaweight;
}

//..Line weighting for point based on particle location..//

std::vector<double> mesh::pt_lineweight(double px,int pindex) const
{
   double dx,dy,temp;
   std::vector<double> lineweight;

   lineweight.push_back(fabs((px-hdeltax)-pmesh[pindex]));
   lineweight.push_back(fabs((px+hdeltax)-pmesh[pindex]));

   return lineweight;
}

//..Line weighting for point based on particle location..//

void mesh::pt_lw(double lineweight[],double px,int pindex) const
{
   lineweight[0] = (fabs((px-hdeltax)-pmesh[pindex]));
   lineweight[1] = (fabs((px+hdeltax)-pmesh[pindex]));
}

//..Line weighting for cell based on particle location..//

void mesh::cell_lw(double lineweight[],double px,int pindex) const
{
   lineweight[0] = (fabs((px-hdeltax)-cmesh[pindex]));
   lineweight[1] = (fabs((px+hdeltax)-cmesh[pindex]));
}



//..Line weighting from point to particle..//

std::vector<double> mesh::pt2part_lineweight(double px,std::vector<int> pindex) const
{
   //double dx,dy,temp;             //EFF
   std::vector<double> lineweight;

   lineweight.push_back((deltax-fabs(px-pmesh[pindex[0]]))/deltax);
   //lineweight.push_back((deltax-fabs(px-pmesh[pindex[1]]))/dx);
   lineweight.push_back(1.0-lineweight[0]);

   return lineweight;
}


//..Line weighting from point to particle..//

void mesh::p2p_lineweight(double px,double lineweight[],int pindex[]) const
{
   lineweight[0] = ((deltax-fabs(px-pmesh[pindex[0]]))/deltax);
   lineweight[1] = (1.0-lineweight[0]);
}


//..Area weighting for cell center based on particle location..//

std::vector<double> mesh::cell_areaweight(double px,double py, int pindex) const
{
   double dx,dy,temp;
   std::vector<double> areaweight;

   dx = (0.5)*meshlength[0]/(numpoints[0]-1);
   dy = (0.5)*meshlength[1]/(numpoints[1]-1);

   areaweight.push_back(squarearea((px-dx),(py-dy),cmesh[meshdims*pindex],cmesh[meshdims*pindex+1]));
   areaweight.push_back(squarearea((px+dx),(py-dy),cmesh[meshdims*pindex],cmesh[meshdims*pindex+1]));
   areaweight.push_back(squarearea((px-dx),(py+dy),cmesh[meshdims*pindex],cmesh[meshdims*pindex+1]));
   areaweight.push_back(squarearea((px+dx),(py+dy),cmesh[meshdims*pindex],cmesh[meshdims*pindex+1]));

   return areaweight;
}

//..Line weighting for cell based on particle location..//

std::vector<double> mesh::cell_lineweight(double px,int pindex) const
{
   double dx,dy,temp;
   std::vector<double> lineweight;

   dx = (0.5)*meshlength[0]/(numpoints[0]-1);

   lineweight.push_back(fabs((px-dx)-cmesh[pindex]));
   lineweight.push_back(fabs((px+dx)-cmesh[pindex]));

   return lineweight;
}


//..Line weighting from cell to particle..//

std::vector<double> mesh::cell2part_lineweight(double px,std::vector<int> pindex)  const
{
   double dx,dy,temp;
   std::vector<double> lineweight;

   dx = meshlength[0]/(numpoints[0]-1);

   lineweight.push_back((dx-fabs(px-cmesh[pindex[0]]))/dx);
   lineweight.push_back((dx-fabs(px-cmesh[pindex[1]]))/dx);

   return lineweight;
}


//..Line weighting from cell to particle..//

void mesh::c2p_lineweight(double px,double lineweight[],int pindex[]) const 
{
   lineweight[0] = ((deltax-fabs(px-cmesh[pindex[0]]))/deltax);
   lineweight[1] = 1-lineweight[0];
}

//..Extrapolate vector from cell to point..//

std::vector<double> mesh::extravec_c2p(std::vector<double> s_vector, int s_index)
{
   int i,j;
   std::vector<int> s_neighbors;
   std::vector<double> s_vecreturn;

   for(i=0;i<(meshdims*meshdims);i++) s_neighbors.push_back(0.0);
   for(i=0;i<(vecdims);i++) s_neighbors.push_back(0.0);

   s_neighbors = cneighofp(s_index);

   for(i=0;i<vecdims;i++) 
   {
      for(j=0;j<s_neighbors.size();j++) s_vecreturn[i] = s_vecreturn[i]+s_vector[s_neighbors[j]*vecdims+i]/(meshdims*meshdims); 
   }

   return s_vecreturn;
}


//..Extrapolate vector from point to cell..//

std::vector<double> mesh::extravec_p2c(std::vector<double> s_vector, int s_index)
{
   int i,j;
   std::vector<int> s_neighbors;
   std::vector<double> s_vecreturn;

   for(i=0;i<(meshdims*meshdims);i++) s_neighbors.push_back(0.0);
   for(i=0;i<(vecdims);i++) s_vecreturn.push_back(0.0);

   s_neighbors = pneighofc(s_index);

   for(i=0;i<vecdims;i++) 
   {
      for(j=0;j<s_neighbors.size();j++) s_vecreturn[i] = s_vecreturn[i]+s_vector[s_neighbors[j]*vecdims+i]/(meshdims*meshdims); 
   }

   return s_vecreturn;
}


//..Extrapolate scalar from cell to point..//

double mesh::extrascal_c2p(std::vector<double> s_scal, int s_index)
{
   int i,j;
   std::vector<int> s_neighbors;
   double s_scalreturn = 0.0;

   for(i=0;i<(meshdims*meshdims);i++) s_neighbors.push_back(0.0);

   s_neighbors = cneighofp(s_index);

   for(j=0;j<s_neighbors.size();j++) s_scalreturn = s_scalreturn + s_scal[s_neighbors[j]*vecdims+i]/vecdims; 

   return s_scalreturn;
}


//..Extrapolate scalar from point to cell..//

double mesh::extrascal_p2c(std::vector<double> s_scal, int s_index)
{
   int i,j;
   std::vector<int> s_neighbors;
   double s_scalreturn = 0;

   for(i=0;i<(meshdims*meshdims);i++) s_neighbors.push_back(0.0);

   s_neighbors = pneighofc(s_index);

   for(j=0;j<s_neighbors.size();j++) s_scalreturn = s_scalreturn + s_scal[s_neighbors[j]*vecdims+i]/vecdims; 

   return s_scalreturn;
}

//..Initially connect particles and mesh..//

void mesh::initconnectpartandmesh(particles & m_part)
{
   int i,j,k;

   //std::cout << "\n\tConnecting Particles and Mesh...";

   std::ofstream wrtfile("connect.out");

   if(meshdims==1)
   {
      for(i=0;i<m_part.pos.size()/meshdims;i++)//  m_part.cell.push_back(nearc(m_part,i)); 
      {  
         m_part.cell.push_back(nearc1Dnew(m_part,i)); 
         //std::cout << i << "Cell: \t" << m_part.pos.size()/meshdims << "\t" << m_part.pos[i] << "\t" << m_part.cell[i-1] << "\t" << nearc(m_part,i) << std::endl;
         //}
         //for(i=0;i<m_part.pos.size()/meshdims;i++) m_part.pt.push_back(nearp(m_part,i));
         m_part.pt.push_back(nearp1Dnew(m_part,i));
         //std::cout << i << "Point: \t" << m_part.pos[i] << "\t" << m_part.pt[i] << "\t" << nearp(m_part,i) << std::endl;
         //wrtfile << i << "\t" << m_part.cell[i] << "\t" << m_part.pt[i] << std::endl;
   }
}
else
{
   for(i=0;i<m_part.pos.size()/meshdims;i++) m_part.cell.push_back(nearc(m_part,i)); 
   for(i=0;i<m_part.pos.size()/meshdims;i++) m_part.pt.push_back(nearp(m_part,i));
}

//for(i=0;i<m_part.cell.size();i++) std::cout << m_part.cell[i] << std::endl; 
//for(i=0;i<m_part.pt.size();i++) std::cout << m_part.pt[i] << std::endl;
//std::cout << "x";
wrtfile.close();
}

//..Connect particles and mesh..//

void mesh::connectpartandmesh(particles & m_part)
{
   int i,j,k;

   //std::cout << "\n\tConnecting Particles and Mesh...";

   if(meshdims==1)
   {
      //#pragma omp parallel  //OMP
      //{
      //#pragma omp for schedule(static)  //OMP
      for(i=0;i<m_part.pos.size()/meshdims;i++)  m_part.cell[i] = nearc1Dnew(m_part,i); 
      //#pragma omp for schedule(static)  //OMP
      for(i=0;i<m_part.pos.size()/meshdims;i++)  m_part.pt[i] = nearp1Dnew(m_part,i);
      //}
   }
   else
   {
      //#pragma omp parallel  //OMP
      //{
      //#pragma omp for schedule(static)  //OMP
      for(i=0;i<m_part.pos.size()/meshdims;i++) m_part.cell[i] = nearc(m_part,i); 
      //#pragma omp for schedule(static)  //OMP
      for(i=0;i<m_part.pos.size()/meshdims;i++) m_part.pt[i] = nearp(m_part,i);
      //}
   }

   //for(i=0;i<m_part.cell.size();i++) std::cout << m_part.cell[i] << std::endl; 
   //for(i=0;i<m_part.pt.size();i++) std::cout << m_part.pt[i] << std::endl;
   //std::cout << "x";
}

//..Disconnect particles and mesh..//

void mesh::disconnectpartandmesh(particles & m_part)
{
   int i,j,k;

   m_part.cell.erase(m_part.cell.begin(),m_part.cell.end());
   m_part.pt.erase(m_part.pt.begin(),m_part.pt.end());
}

//..Find mesh area..//

void mesh::mesharea(fields m_EM)
{
   int i,j,k;
   double r_area;

   if(Bloc==0)
   {
      if(q1dflag==0)
      {
         for(i=0;i<pmesh.size();i++)  parea.push_back(areain);
         for(i=0;i<cmesh.size();i++)  carea.push_back(areain);
         //for(i=0;i<pmesh.size();i++)  prarea.push_back(areain);
         //for(i=0;i<cmesh.size();i++)  crarea.push_back(areain);
         //for(i=0;i<cmesh.size();i++)  area.push_back(fabs(m_EM.B[0]*areain/m_EM.B[vecdims*i]));
         //for(i=0;i<cmesh.size();i++)  area.push_back(fabs(m_EM.B[0]*areain/m_EM.B[vecdims*i]));
      }
      else if(q1dflag==1)
      {
#if SIMTYPE==1
         for(i=0;i<cmesh.size();i++)  carea.push_back(areain);  //..MB check
         for(i=0;i<pmesh.size();i++)  parea.push_back(areain);  //..MB check
#else
         for(i=0;i<cmesh.size();i++)  carea.push_back(fabs(m_EM.B[0]*areain/m_EM.B[vecdims*i]));
         for(i=0;i<pmesh.size();i++)  parea.push_back(0.5*(carea[i]+carea[i+1]));
#endif
         /*for(i=0;i<cmesh.size();i++) 
           {
           r_area = 2.0*pi*sqrt(area[i]/pi)*deltax;
           rarea.push_back(r_area);
           }*/
      }
   }
   else if(Bloc==1)
   {
      if(q1dflag==0)
      {
         for(i=0;i<cmesh.size();i++)  carea.push_back(areain);
         for(i=0;i<pmesh.size();i++)  parea.push_back(areain);
         //for(i=0;i<pmesh.size();i++)  rarea.push_back(areain);
      }
      else if(q1dflag==1)
      {
         for(i=0;i<(pmesh.size());i++)  parea.push_back(fabs(m_EM.B[0]*areain/(m_EM.B[vecdims*i])));

         carea.push_back(areain);
         for(i=0;i<(cmesh.size()-2);i++)  carea.push_back(0.5*(parea[i]+parea[i+1]));
         carea.push_back(parea[pmesh.size()-2]);
         /*for(i=0;i<pmesh.size();i++) 
           {
           r_area = 2.0*pi*sqrt(area[i]/pi)*deltax;
           rarea.push_back(r_area);
           }*/
      }
   }

   //std::cout << std::endl << std::endl << "Area:   ";
   //for(i=0;i<parea.size();i++) std::cout << pmesh[i*meshdims] << "\t" << parea[i] << "\n";
   //std::cout << std::endl << std::endl;  

}

//....Create array of neighbors to cell corners......//

void mesh::meshPointNeighbors(std::vector<int> & m_pofpneigh, std::vector<int> & m_cofpneigh)
{
   int i,j;
   std::vector<int> neighbors;

   for(i=0;i<2*meshdims;i++) neighbors.push_back(0);

   for(i=0;i<pmesh.size();i++)  
   {
      neighbors = pneighofp(i);
      for(j=0;j<neighbors.size();j++)   pofpneigh.push_back(neighbors[j]);

      neighbors = cneighofp(i);
      for(j=0;j<neighbors.size();j++)   cofpneigh.push_back(neighbors[j]);
   }

   //for(i=0;i<pofpneigh.size();i++) std::cout << pofpneigh[i] << "\t";

}

//....Create array of neighbors to cell centers......//

void mesh::meshCenterNeighbors(std::vector<int> & m_cofcneigh, std::vector<int> & m_pofcneigh)
{
   int i,j;
   std::vector<int> neighbors;

   for(i=0;i<2*meshdims;i++) neighbors.push_back(0);

   for(i=0;i<cmesh.size();i++)  
   {
      neighbors = cneighofc(i);
      for(j=0;j<neighbors.size();j++)   m_cofcneigh.push_back(neighbors[j]);

      neighbors = pneighofc(i);
      for(j=0;j<neighbors.size();j++)   m_pofcneigh.push_back(neighbors[j]);
   }

}

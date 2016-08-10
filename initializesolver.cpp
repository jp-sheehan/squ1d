#include "initializesolver.h"
#include <ctime>

//....Reads in type of solver....//

int initializeSolver::readsolver()             
{
   std::string temp;

   //....Solver type file....//

   rd_solverfilename = "SolverType.inp"; 

   std::cout << "Reading solver type...";

   std::ifstream rdfile(rd_solverfilename.c_str());

   rdfile >> temp >> setype >> sptype;

   std::cout << "\n\n" << setype << " " << sptype <<  " solver engaging...\n\n";

   rdfile.close();

   if(setype == "ELECTROSTATIC" && sptype == "PIC") sflag = 0;
   //else if(setype == "ELECTROMAGNETIC" && sptype == "HYBRID") sflag = 2;
   //else if(setype == "EULER" && sptype == "1D") sflag = 10;
   //else if(setype == "EULER" && sptype == "2D") sflag = 11;
   else
   {
      std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
      std::cout << "Only ELECTROSTATIC PIC is allowed." << std::endl;
      exit(EXIT_FAILURE); 
   }

   //....Returns solver flag....//
   return sflag;            
}



//....Read data for specific run and solver type from "SolverType.inp" and "SolverInput.inp" 

void initializeSolver::readdata(solverVars & i_svar, mesh & i_msh, std::vector<boundvars> & i_bdv, std::vector<spclvars> &i_spclv) 
{
   std::string temp, tempstring; // variables for capturing junk
   int i,j,k,l,tempint,nbcs;
   double tempdouble;
   double tempbool;
   std::vector<double> tempvec;
   boundvars tempbdv;
   spclvars tempspclv;
   std::string strang;

   boundary i_bnd;

   rd_solverfilename = "SolverType.inp";
   rd_inputfilename = "SolverInput.inp";
   rd_solverlogfilename = rd_solverfilename +".log";

   //....Reading Solver Type....//

   std::ifstream rdfile(rd_solverfilename.c_str());
   std::ifstream infile(rd_inputfilename.c_str());

   rdfile >> temp >> setype >> sptype;
   rdfile >> temp >> rd_meshread;
   rdfile >> temp >> rd_meshwrite;
   //rdfile >> temp >> rd_meshfilename;
   //rdfile >> temp >> rd_meshfileformat;
   rdfile >> temp >> i_svar.rst;
   rdfile >> temp >> i_msh.meshdims >> i_msh.sym;
   rdfile >> temp >> i_msh.vecdims;

   if(rd_meshwrite == "SIMPLE") i_msh.wrtflag = 0; 
   else if(rd_meshwrite == "TECPLOT") i_msh.wrtflag =1; 
   else
   {
      std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
      std::cout << "rd_meshwrite: " << rd_meshwrite << std::endl;
      exit(EXIT_FAILURE); 
   }

   if(i_msh.sym == "PLANAR")
   {
      i_svar.q1dflag = 0;
      i_msh.q1dflag = 0;
      i_msh.axisflag = 0;
   }
   else if(i_msh.sym == "AXISYMMETRIC")
   { 
      i_svar.q1dflag = 1;
      i_msh.q1dflag = 1;
      i_msh.axisflag = 1;
   }
   else
   {
      std::cout << "\nExiting...ERROR in SYMMETRY\n";
      std::cout << "Dimensions: " << i_msh.meshdims << std::endl;
      std::cout << "Symmetry: " << i_msh.sym << std::endl;
      exit(EXIT_FAILURE); 
   }



   //....Different types of solvers....//
   if(setype == "ELECTROSTATIC" && sptype == "PIC") sflag = 0;
   /* not implemented
      else if(setype == "ELECTROMAGNETIC" && sptype == "HYBRID") sflag = 2;
      else if(setype == "EULER" && sptype == "1D") sflag = 10;
      else if(setype == "EULER" && sptype == "2D") sflag = 11;
      */
   else
   {
      std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
      std::cout << "setype: " << setype << std::endl;
      std::cout << "sptype: " << sptype << std::endl;
      exit(EXIT_FAILURE); 
   }

   i_msh.sflag=sflag;

   if (i_msh.meshdims == 1)
   { 
      std::cout << "Mesh Type:   1D " << i_msh.sym << " with   " << i_msh.vecdims << " Vector Dimensions \n\n"; 
   }
   /* not implemented
      else if (i_msh.meshdims == 2)
      {
      std::cout << "Mesh Type:   2D  " << i_msh.sym << " with " << i_msh.vecdims << " Vector Dimensions \n\n"; 
      }
      else if (i_msh.meshdims == 3)
      {
      std::cout << "Mesh Type:   3D " << i_msh.sym << " with   " << i_msh.vecdims << " Vector Dimensions \n\n"; 
      }
      */
   else
   {
      std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
      std::cout << "Number of mesh dimension must be 1.  Input meshdims: " << i_msh.meshdims << std::endl;
      exit(EXIT_FAILURE); 
   }

   if (rd_meshread == 0) // Create mesh
   {
      // Set mesh length
      rdfile >> temp;
      std::cout << "Length:\t \t";

      for (i=0; i<i_msh.meshdims; i++)
      {
         rdfile >> tempdouble;
         i_msh.meshlength.push_back(tempdouble);  
         std::cout << i_msh.meshlength[i] << "\t";
      }

      // Set mesh starting point
      rdfile >> temp;
      std::cout << "\nStart Point:\t \t";

      for (i=0; i<i_msh.meshdims; i++)
      {
         rdfile >> tempdouble;
         i_msh.meshstart.push_back(tempdouble);  
         std::cout << i_msh.meshstart[i] << "\t";
      }

      // Set mesh boundaries (I think)
      for(i=0;i<i_msh.meshdims;i++)
      {
         i_msh.meshbd.push_back(i_msh.meshstart[i]);
         i_msh.meshbd.push_back(i_msh.meshlength[i]+i_msh.meshstart[i]);
      }

      // set number of points on mesh
      rdfile >> temp;
      std::cout << "\nGrid:\t \t";

      for (i=0; i<i_msh.meshdims; i++)
      {
         rdfile >> tempint;
         i_msh.numpoints.push_back(tempint);  
         std::cout << i_msh.numpoints[i] << "\t";
      }

      // Set cell weight
      rdfile >> temp;
      std::cout << std::endl << "Cell Weight:\t";

      for (i=0; i<i_msh.meshdims; i++) 
      { 
         rdfile >> tempdouble;
         i_msh.cellweight.push_back(tempdouble);  
         std::cout << i_msh.cellweight[i] << "\t";
      }

      // Set cell area
      rdfile >> temp;
      std::cout << std::endl << "Cell Area:\t";

      rdfile >> tempdouble;
      areain = tempdouble;  
      i_msh.areain = areain;
      std::cout << areain;

      //......Mesh Shifting Algorithm....Not currently included

      /*rdfile >> temp >>  mesh_num;

        std::cout << "\n\nNumber of Meshes:  \t" << mesh_num; 
        std::cout << "\nShifts: ";

        rdfile >> temp;

        for(i=0;i<mesh_num;i++)
        {
        for(j=0;j<i_msh.meshdims;j++)
        {
        rdfile >> tempdouble;
        mesh_shift.push_back(tempdouble);
        std::cout  << mesh_shift[i*imsh.meshdims+j] << "\t";
        }
        std::cout << std::endl << "\t";
        } */

      // Set number of ghost cells
      rdfile >> temp >> i_msh.numghost;

      // Set magnetic field location
      rdfile >> temp >> i_msh.Bloc;

      std::cout << std::endl;     

      if(sflag == 0)
      {
         // Set electric potential location
         rdfile >> i_msh.philoc;

         if (i_msh.philoc == 0)   i_msh.Eloc = 1;
         else if(i_msh.philoc ==1)  i_msh.Eloc = 1;
         if (i_msh.philoc == 0)   std::cout << "Potential at Cell Centers\n";
         else if(i_msh.philoc ==1)  std::cout << "Potential at Cell Corners\n";
         if (i_msh.Bloc == 0)   std::cout << "Magnetic Fields at Cell Centers\n";
         else if(i_msh.Bloc ==1)  std::cout << "Magnetic Fields at Cell Corners\n";
      }
      /*
      else if(sflag==2)
      {
         rdfile >> i_msh.Eloc;
         if (i_msh.Eloc == 0)   std::cout << "Electric Fields at Cell Centers\n";
         else if(i_msh.Eloc ==1)  std::cout << "Electric Fields at Cell Corners\n";
         if (i_msh.Bloc == 0)   std::cout << "Magnetic Fields at Cell Centers\n";
         else if(i_msh.Bloc ==1)  std::cout << "Magnetic Fields at Cell Corners\n";
      }
      else if(sflag==10)
      {
         rdfile >> i_msh.Eloc;
         if (i_msh.Eloc == 0)   std::cout << "No Electric or Magnetic Fields\n";
         else if(i_msh.Eloc ==1)  std::cout << "No Electric or Magnetic Fields\n";
      }
      else if(sflag==11)
      {
         rdfile >> i_msh.Eloc;
         if (i_msh.Eloc == 0)   std::cout << "No Electric or Magnetic Fields\n";
         else if(i_msh.Eloc ==1)  std::cout << "No Electric or Magnetic Fields\n";
      }
      */

      // Set particle distribution
      rdfile >> temp >> i_msh.pdist;

      if(i_msh.pdist == 0) std::cout << "\nParticles distributed over entire domain\n";
      else if(i_msh.pdist==1) std::cout << "\nParticles distributed cell by cell\n";

      // Set smoothing
      rdfile >> temp >> i_msh.smthcnt;
      if(i_msh.pdist == 0) std::cout << "\nSmoothed "<< i_msh.smthcnt << " times\n";
   }


   //.....Read information for outputs.......//

   rdfile >> temp >> i_svar.outavg;  
   if(i_svar.outavg==1) std::cout << "\nPrinting out volume averages\n";
   rdfile >> temp >> i_svar.outcont;  
   if(i_svar.outcont==1) std::cout << "\nPrinting out continuum\n";
   rdfile >> temp >> i_svar.outpart;  
   if(i_svar.outpart==1) std::cout << "\nPrinting out particles\n";
   rdfile >> temp >> i_svar.outvdf >> i_svar.outvdf_flag;  
   if(i_svar.outvdf==1) std::cout << "\nPrinting out vdfs\n";
   rdfile >> temp >> i_svar.outrst;  
   if(i_svar.outrst==1) std::cout << "\nPrinting out restart files\n";
   rdfile >> temp >> i_svar.outfinal >> i_svar.outfinalskip;  
   if(i_svar.outfinal>=1) std::cout << "\nPrinting out " << i_svar.outfinal << " time steps skipping every " << i_svar.outfinalskip;


   //......Read solver specific information......//

   rdfile >> temp >> i_svar.lscale;
   rdfile >> temp >> i_svar.tscale;
   rdfile >> temp >> i_svar.cfl;
   rdfile >> temp >> i_svar.dt;
   rdfile >> temp >> i_svar.iter;
   rdfile >> temp >> i_svar.p_iter;

   std::cout << "Length scale:\t" << i_svar.lscale << "\t Time scale:\t " << i_svar.tscale << "\t CFL:\t" << i_svar.cfl << "\t dt:\t" << i_svar.dt << "\t Iterations:\t" << i_svar.iter << std::endl;

   if(sflag<10) // if solver is electrostatic (or electromagnetic (not implemented))
   {
      // Set particle weight
      rdfile >> temp >> i_svar.pwght;
      std::cout <<  "Particle Weight:\t" << i_svar.pwght <<  std::endl; 

      // Set interpolation scheme
      rdfile >> temp >> temp;
      std::cout <<  "Interpolation Scheme:\t" << temp <<  std::endl; 

      if(temp=="ZEROTH") i_msh.intscheme = 0; // nearest
      else if(temp=="FIRST") i_msh.intscheme = 1; // linear
      else
      {
         std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
         std::cout << "Interpolation scheme must be ZEROTH or FIRST." << std::endl
            << "Inerpolation scheme: " << temp << std::endl;
         exit(EXIT_FAILURE); 
      }

      // Set particle mover algorithm
      rdfile >> temp >> temp;
      std::cout <<  "Particle Mover:\t" << temp <<  std::endl; 

      if(temp=="SIMPLE") i_msh.mvscheme = 1;
      else if(temp=="BORIS") i_msh.mvscheme = 2;
      else if(temp=="Q1D") i_msh.mvscheme = 3;
      else
      {
         std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
         std::cout << "Particle mover must be SIMPLE, BORIS, or Q1D." << std::endl
            << "Particle mover: " << temp << std::endl;
         exit(EXIT_FAILURE); 
      }

      i_svar.mvscheme=i_msh.mvscheme;

      // Set background density charge
      rdfile >>  temp >> rhoback;
      std::cout << "Background Density:\t"<< rhoback << std::endl;

      // End of input file
      rdfile >> temp;

      if(temp !="EOF") 
      {
         std::cout << "\nExiting...EOF ERROR in SOLVER TYPE\n";
         std::cout << "No inputs should be included after the background density." << std::endl;
         exit(EXIT_FAILURE); 
      }

      std::cout << std::endl << std::endl;

      rdfile.close();

      //....Read solver input file....//

      std::cout << "Reading input file...\n";

      infile  >>  temp >>  i_svar.nsp; 
      infile  >> temp >> temp >> temp >> temp; //>> temp;

      nsp = i_svar.nsp; // Number of species

      for(i=0;i<i_svar.nsp;i++) // iterate of number of species
      {
         infile >> temp;
         sname.push_back(temp);
         infile >> tempdouble;
         scharge.push_back(tempdouble);
         infile >> tempdouble;
         smass.push_back(tempdouble);
         infile >> tempbool;
         smag.push_back(tempbool);
         //infile >> tempint;
         //scycIT.push_back(tempint);
      }

      //....Read in Coulomb Collisionality Information.....//

      infile  >>  temp >> i_svar.coul;

      //....Read in Collisionality Information.....//

      infile  >>  temp >> i_svar.nct;
      nct = i_svar.nct; // Number of collision types

      if(i_svar.nct>0)
      {
         infile >> temp >> temp >> temp;

         for(i=0;i<i_svar.nct;i++)
         {
            infile >> temp;
            scllsnname.push_back(temp);
            infile >> temp;
            scllsntype.push_back(temp);
            infile >> temp;
            scrsstype.push_back(temp);

            if(scrsstype[i] == "CONSTANT")
            {
               infile >> tempdouble;
               scrsssctn.push_back(tempdouble);
            }
            else scrsssctn.push_back(0.0);

            //if(scllsntype[i] == "NEUTRAL")
            //{
            //}
         }
         infile >> temp >> neutdens>> temp >> neuttemp;

         std::cout << "\n\nCollisions ON \n";
         for(i=0;i<i_svar.nct;i++) std::cout << "\t" << scllsntype[i] << std::endl;
      }

      //....Read in boundary conditions....//

      infile >> temp >> i_svar.nbound;
      i_msh.nbound = i_svar.nbound;

      std::cout << "\n\nNumber of Boundaries: \t \n"  <<  i_svar.nbound << std::endl ;

      nbcs = 3*nsp+2;

      for(i=0; i<i_svar.nbound; i++)
      {
         tempbdv.nsp = i_svar.nsp; // Number of species
         tempbdv.nDN = 0; // Number of Dirichlet and Neumann boundary conditions
         tempbdv.nbound = i_svar.nbound;

         infile >> temp >> temp;
         tempbdv.boundname = temp;

         std::cout << std::endl << tempbdv.boundname << ": \n";  //Name of boundary

         for(k=0; k<i_msh.meshdims; k++)
         {
            infile >> temp;
            for(j=0; j<i_msh.meshdims; j++)
            {
               infile >> tempdouble;
               tempbdv.boundrange.push_back(tempdouble); //Range of boundary
            }
         }

         infile >> temp >> temp;
         tempbdv.clss = temp;       //Read Boundary classification

         if(tempbdv.clss=="PERIODIC") 
         {
            i_msh.perflag=1;
            tempbdv.wallflag = 0;
            std::cout << std::endl << "Type:  " << tempbdv.clss << std::endl; 
         }
         else if(tempbdv.clss=="WALL") 
         {
            i_msh.perflag=0;
            tempbdv.wallflag = 1;
            infile >> temp;
            tempbdv.walltype = temp;
            if(temp == "FLOATING") 
            {
               tempbdv.wallfloat = 1;
               tempbdv.wallcur = 0;
               tempbdv.wallcap = 0;
               tempbdv.wallvolt = 0;
            }
            else if(temp == "CURRENT") 
            {
               tempbdv.wallcur = 1;
               tempbdv.wallfloat = 0;
               tempbdv.wallcap = 0;
               tempbdv.wallvolt = 0;
               infile >> tempdouble;
               tempbdv.Jfreq = tempdouble;
            }
            else if(temp == "VOLTAGE") 
            {
               tempbdv.wallcur = 0;
               tempbdv.wallfloat = 0;
               tempbdv.wallcap = 0;
               tempbdv.wallvolt = 1;
               infile >> tempdouble;
               tempbdv.Vfreq = tempdouble;
            }
            else if(temp == "CAPACITOR") 
            {
               tempbdv.wallcap = 1;
               tempbdv.wallfloat = 0;
               tempbdv.wallcur = 0;
               tempbdv.wallvolt = 0;
               //infile >> tempdouble;
               //tempbdv.cap = tempdouble;
            }
            else
            {  
               tempbdv.wallcap = 0;
               tempbdv.wallfloat = 0;
               tempbdv.wallcur = 0;
               tempbdv.wallvolt = 0;
            }
            std::cout << std::endl << "Type:  " << tempbdv.clss << " " << tempbdv.walltype << std::endl;  
         } 
         else if(tempbdv.clss=="SOURCE") 
         {
            i_msh.perflag=0;
            infile >> temp;
            tempbdv.srctype = temp;
            std::cout << std::endl << "Type:  " << tempbdv.clss << " " << tempbdv.srctype  << std::endl;  
         } 
         else 
         {
            i_msh.perflag=0;
            tempbdv.wallflag = 0;
            std::cout << std::endl << "Type:  " << tempbdv.clss << std::endl;
         }

         for(l=0;l<nsp;l++)
         {
            infile >> temp >> temp;
            tempbdv.boundtype.push_back(temp);    //Read boundary type

            std::cout << sname[l]  <<" Density BC: " << tempbdv.boundtype[3*l] << "\n";

            if(tempbdv.boundtype[3*l] == "DIRICHLET" || tempbdv.boundtype[3*l] == "NEUMANN" )  //Read in strings of equations for DIRICHLET BC
            { 
               infile  >>  temp;
               tempbdv.bddens.push_back(temp);
               tempbdv.bddens_val.push_back(0.0);
               infile  >>  temp;
               tempbdv.bdddist.push_back(temp);
               tempbdv.nDN = tempbdv.nDN + 1;
            }
            else if(tempbdv.boundtype[3*l] == "PERIODIC")  
            { 
               tempbdv.bddens_val.push_back(0.0);
            } 

            infile >> temp >> temp;
            tempbdv.boundtype.push_back(temp); 

            std::cout << sname[l] << " Temperature BC: " << tempbdv.boundtype[1+3*l] << "\n";

            if(tempbdv.boundtype[1+3*l] == "DIRICHLET" || tempbdv.boundtype[1+3*l] == "NEUMANN")
            {  
               infile  >>  temp;
               tempbdv.bdtemp.push_back(temp);
               infile  >>  temp;
               tempbdv.bdthdist.push_back(temp);
               tempbdv.nDN = tempbdv.nDN + 1;
            }

            infile >> temp >> temp;
            tempbdv.boundtype.push_back(temp); 

            std::cout << sname[l] <<" Velocity BC: " << tempbdv.boundtype[2+3*l] << "\n";

            if(tempbdv.boundtype[2+3*l] == "DIRICHLET" || tempbdv.boundtype[2+3*l] == "NEUMANN")
            {  
               for(j=0;j<i_msh.vecdims;j++)
               {
                  infile  >>  temp;
                  tempbdv.bdvel.push_back(temp);
               }
               tempbdv.nDN = tempbdv.nDN + 1;
            }

         }

         if(sflag==0)
         {
            infile >> temp >> temp;
            tempbdv.boundtype.push_back(temp); 

            std::cout << "Potential BC: " << tempbdv.boundtype[3*nsp] << "\n";

            if(tempbdv.boundtype[3*nsp] == "DIRICHLET" || tempbdv.boundtype[3*nsp] == "NEUMANN")
            {  
               infile  >>  temp;
               tempbdv.bdphi = temp;
               tempbdv.bdphi_val = 0.0;
               tempbdv.nDN = tempbdv.nDN + 1;
            }
            else if(tempbdv.boundtype[3*nsp] == "PERIODIC")  
            { 
               infile  >>  temp;
               tempbdv.bdphi = temp;
               tempbdv.bdphi_val = 0.0;
               tempbdv.nDN = tempbdv.nDN + 1;
            } 
         }
         else if(sflag==2)
         {
            std::cout << "Efield BC: " << tempbdv.boundtype[3*nsp] << "\n";

            if(tempbdv.boundtype[3*nsp] == "DIRICHLET"  || tempbdv.boundtype[3*nsp] == "NEUMANN" )
            {  
               for(j=0;j<i_msh.vecdims;j++)
               {
                  infile  >>  temp;
                  tempbdv.bdE.push_back(temp);
               }
               tempbdv.nDN = tempbdv.nDN + 1;
            }
         }

         infile >> temp >> temp;
         tempbdv.boundtype.push_back(temp); 
         std::cout << "Bfield BC: " << tempbdv.boundtype[3*nsp+1] << "\n";

         if(tempbdv.boundtype[3*nsp+1] == "DIRICHLET" || tempbdv.boundtype[3*nsp+1] == "NEUMANN")
         {  
            for(j=0;j<i_msh.vecdims;j++)
            {
               infile  >>  temp;
               tempbdv.bdB.push_back(temp);
            }
            tempbdv.nDN = tempbdv.nDN + 1;
         }

         i_bdv.push_back(tempbdv);
         i_bnd.clearallbdv(tempbdv);

         /*std::cout << "\nVals:\n";
           std::cout << std::endl << i_bdv[i].boundname << ": \n";  //Name of boundary
           std::cout << std::endl << "Type:  " << i_bdv[i].clss << std::endl; 
           std::cout << sname[0]  <<" Density BC: " << i_bdv[i].boundtype[3*0] << "\n";
           std::cout << sname[1]  <<" Density BC: " << i_bdv[i].boundtype[3*1] << "\n";
           std::cout << sname[0] << " Temperature BC: " << i_bdv[i].boundtype[1+3*0] << "\n";
           std::cout << sname[1] << " Temperature BC: " << i_bdv[i].boundtype[1+3*1] << "\n";
           std::cout << sname[0] <<" Velocity BC: " << i_bdv[i].boundtype[2+3*0] << "\n";
           std::cout << sname[1] <<" Velocity BC: " << i_bdv[i].boundtype[2+3*1] << "\n";
           std::cout << "Potential BC: " << i_bdv[i].boundtype[3*nsp] << "\n";
           std::cout << "Efield BC: " << i_bdv[i].boundtype[3*nsp] << "\n";
           std::cout << "Bfield BC: " << i_bdv[i].boundtype[3*nsp+1] << "\n";*/

      }

   }

   /*
   if( sflag >= 10)          //....Euler(From CFD Project)....//
   {
      infile  >>  temp;
      infile  >>  temp  >> temp;

      std::cout << temp;

      if(temp == "RK4") i_svar.tdflag = 1;

      infile  >>  temp  >> temp; 
      std::cout << temp;
      if(temp == "LF")  i_svar.sdflag = 1;
      else if(temp == "CENTRAL")  i_svar.sdflag = 0;
      else if(temp == "MUSCL")  i_svar.sdflag = 2;
      else
      {
         std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
         exit(EXIT_FAILURE); 
      }

      std::cout << "\n\nSolver Type:\t";
      if(i_svar.tdflag==1) std::cout << "RK4" << "\t" ;
      else if(i_svar.sdflag==1) std::cout << "LF" ;
      else if(i_svar.sdflag==0) std::cout << "CENTRAL" ;
      else if(i_svar.sdflag==2) std::cout << "MUSCL" ;
      else
      {
         std::cout << "\nExiting...ERROR in SOLVER TYPE\n";
         exit(EXIT_FAILURE); 
      }

      infile  >>  temp;
      infile  >>  temp  >> i_svar.atmwght;  
      infile  >>  temp  >> i_svar.gam; 

      std::cout << "\n\nGas Properties:\t" << i_svar.atmwght << "\t" << i_svar.gam ;

      infile >> temp >> i_bnd.nbound;

      std::cout << "\n\nNumber of Boundaries: \t \n"  <<  i_bnd.nbound << std::endl ;

      i_bnd.nDirichlet = 0;

      for(i=0; i<i_bnd.nbound; i++)
      {
         infile >> temp;
         i_bnd.boundname.push_back(temp);

         std::cout << std::endl << i_bnd.boundname[i] << ": \n";


         for(k=0; k<i_msh.meshdims; k++)
         {
            infile >> temp;
            for(j=0; j<i_msh.meshdims; j++)
            {
               infile >> tempdouble;
               i_bnd.boundrange.push_back(tempdouble);
            }
         }

         infile >> temp >> temp;
         i_bnd.boundtype.push_back(temp); 

         std::cout << "\tDensity BC: " << i_bnd.boundtype[3*i] << "\n";

         if(i_bnd.boundtype[3*i] == "DIRICHLET")
         {  
            infile  >>  temp;
            i_bnd.bdrho.push_back(temp);
            i_bnd.nDirichlet = i_bnd.nDirichlet + 1;
         }

         infile >> temp >> temp;
         i_bnd.boundtype.push_back(temp); 

         std::cout << "\tVelocity BC: " << i_bnd.boundtype[3*i+1] << "\n";

         if(i_bnd.boundtype[3*i+1] == "DIRICHLET")
         {  
            for(j=0;j<i_msh.vecdims;j++)
            {
               infile  >>  temp;
               i_bnd.bdU.push_back(temp);
            }
            i_bnd.nDirichlet = i_bnd.nDirichlet + 1;
         }

         infile >> temp >> temp;
         i_bnd.boundtype.push_back(temp); 

         std::cout << "\tPressure BC: " << i_bnd.boundtype[3*i+2] << "\n";

         if(i_bnd.boundtype[3*i+2] == "DIRICHLET")
         {  
            infile  >>  temp;
            i_bnd.bdp.push_back(temp);
            i_bnd.nDirichlet = i_bnd.nDirichlet + 1;
         }
      }

   }
   */


   if(sflag<10)
   {
      infile  >>  temp;

      for(i=0;i<nsp;i++)  
      {
         infile  >>  temp  >>  temp;
         init_dens.push_back(temp);
         infile >>  temp;
         dens_dist.push_back(temp);
         infile >> tempdouble;
         dens_pert.push_back(tempdouble);
         infile  >> temp >> temp;
         init_Temp.push_back(temp);
         infile >> temp;
         therm_dist.push_back(temp);
         infile >> temp;
         for(j=0;j<i_msh.vecdims;j++)
         {
            infile >>  temp;
            init_U.push_back(temp);
         }
      }     

      infile  >>  temp;

      if(sflag==0)
      {
         infile >>  init_phi;
      }
      else if(sflag==2)
      {
         for(i=0; i<i_msh.vecdims; i++)
         {
            infile >>  temp;
            init_E.push_back(temp);
         }
      }

      infile  >> temp;

      if(temp == "CurrentLoop") //Read in current loop information to be used for magnetic field
      {
         i_svar.flag_CL = true; 
         infile >> loopradius >> loopcurrent >> loopcenter[0] >> loopcenter[1] >> loopcenter[2];
      }
      else
      {
         i_svar.flag_CL = false;
         for(i=0; i<i_msh.vecdims; i++)
         {
            infile >>  temp;
            init_B.push_back(temp);
         }
      }

      std::cout << "\nInitial Conditions: \n";

      for(i=0;i<nsp;i++)
      {
         std::cout << "\n\nSPECIES:   " << sname[i] << std::endl;
         std::cout << dens_dist[i] << " DENSITY:  " <<  init_dens[i] << "\n";
         std::cout << therm_dist[i] << " TEMPERATURE:  " <<  init_Temp[i] << "\n";
         std::cout << "VELOCITY:  ";
         for(j=0;j<i_msh.vecdims;j++) std::cout << "\t" <<  init_U[i*i_msh.vecdims+j];
         std::cout << std::endl;
      }

      std::cout << "FIELDS:\n" << std::endl;

      if(i_svar.flag_CL==true)
      {
         std::cout<< "\nCurrent loop: r=" << loopradius << ", I = " << loopcurrent << ",  (x,y,z) = (" << loopcenter[0] << "," << loopcenter[1] << "," << loopcenter[2] << ")\n";
      }
      else
      {
         std::cout << "\nMagnetic Field:\t" ;

         for(i=0;i<i_msh.vecdims;i++)
         {
            std::cout <<  init_B[i] << "\t";
         }
      }

      if(sflag==0)
      {
         std::cout << "\nPotential: \t" ;
         std::cout <<  init_phi << "\n" ;
      }
      else if(sflag==2)
      {
         std::cout << "\nEfield: \t" ;
         for(i=0;i<i_msh.vecdims;i++)
         {
            std::cout <<  init_E[i] << "\t" ;
         }
      }

      std::cout << std::endl;  



      //....Read in special regions....//


      infile  >>  temp  >> i_svar.nspcl;

      std::cout << "\n\nNumber of Special Regions: \t"  <<  i_svar.nspcl << std::endl ;

      for(i=0;i<i_svar.nspcl;i++)
      {
         tempspclv.nspcl = i_svar.nspcl;
         infile  >> temp >> tempspclv.spcltype ;

         if(tempspclv.spcltype=="SOURCE")
         {

            for(j=0;j<i_msh.meshdims;j++)
            {
               infile >> temp >> tempdouble;
               tempspclv.spclrange.push_back(tempdouble);
               infile >> tempdouble;
               tempspclv.spclrange.push_back(tempdouble);
            }

            for(j=0;j<i_svar.nsp;j++)
            {
               infile >> temp >> temp;
               tempspclv.spclFlux.push_back(temp);

               infile >> temp;
               tempspclv.spclddist.push_back(temp);

               infile >> temp >> temp;
               tempspclv.spclT.push_back(temp);

               infile >> temp;
               tempspclv.spclthdist.push_back(temp);

            }
         }
         else if(tempspclv.spcltype=="HEATING")
         {
            for(j=0;j<i_msh.meshdims;j++)
            {
               infile >> temp >> tempdouble;
               tempspclv.spclrange.push_back(tempdouble);
               infile >> tempdouble;
               tempspclv.spclrange.push_back(tempdouble);
            }

            infile >> temp >> tempdouble;
            tempspclv.spclJ = tempdouble;

            infile >> temp >> tempdouble;
            tempspclv.spclomega = tempdouble;
         }
         else if(tempspclv.spcltype=="EFIELD")
         {
            for(j=0;j<i_msh.meshdims;j++)
            {
               infile >> temp >> tempdouble;
               tempspclv.spclrange.push_back(tempdouble);
               infile >> tempdouble;
               tempspclv.spclrange.push_back(tempdouble);
            }

            infile >> temp >> tempdouble;
            tempspclv.spclE = tempdouble;

            infile >> temp >> tempdouble;
            tempspclv.spclomega = tempdouble;
         }

         i_spclv.push_back(tempspclv);
         clearallspclv(tempspclv);
      }

      for(j=0;j<i_svar.nspcl;j++)
      {
         std::cout << i_spclv[j].spcltype <<  " Region from "  << i_spclv[j].spclrange[0] << " to " << i_spclv[j].spclrange[1] << std::endl;
      }

   }
   /*
   else if(sflag>=10)
   {
      infile  >>  temp;
      infile  >>  temp  >>  init_rho;
      infile  >>  temp  >>  init_p;

      infile  >> temp;

      for(i=0;i<i_msh.vecdims;i++)     
      {
         infile >> temp;
         init_U.push_back(temp);
      }

      std::cout << "\nInitial Conditions: \n";
      std::cout << "Density: \t" <<  init_rho << "\n" << "Pressure: \t" << init_p << "\n"<<  "Velocity: \t";

      for(i=0;i<i_msh.vecdims;i++) std::cout << init_U[i] << "\t";
      std::cout << std::endl;
   }
   */

   infile >> temp;

   if(temp !="EOF") 
   {
      std::cout << "\nExiting...EOF ERROR in SOLVER INPUT\n";
      std::cout << "Expected EOF.  Found: " << temp << std::endl;
      exit(EXIT_FAILURE); 
   }

   infile.close();

}

//....Initialize Euler solver (CFD Class Project)....//
/*
   void initializeSolver::initializeEuler(flow & i_cont, mesh &i_mesh,boundary & i_bound, solverVars i_svar)
   {
   int i;

   i_cont.gam = i_svar.gam;

   initializeFluid(i_cont,i_mesh.cmesh,i_mesh.meshdims,i_mesh.vecdims);
   }
   */

//....Initialize PIC Solver (call routines for domain,particles,boundaries, etc)...//

void initializeSolver::initializePIC(std::vector<particles> &i_prt, fields & i_EM, std::vector<contnm> &i_cnt, mesh &i_mesh,std::vector<boundvars> & i_bdv, std::vector<spclvars> &i_spclv, solverVars &i_svar, writeOutput &i_wrt,MPIvars &i_mpiv)
{
   int i,j;
   int numprocs,procid; //MPI
   double temp;

   solver slvr;

   boundary tmp_bound;
   particles tmp_prt;
   contnm tmp_cnt;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   /* conditions are the same
   if(i_mesh.philoc==0) 
   {
      initializephi(i_EM.phi,i_mesh);   //CHG
      initializephitoE(i_EM.phi,i_EM.E,i_mesh);   //CHG
   }
   else if(i_mesh.philoc==1) 
   {
      initializephi(i_EM.phi,i_mesh);  //CHG
      initializephitoE(i_EM.phi,i_EM.E,i_mesh); //CHG
   }
   */
   initializephi(i_EM.phi,i_mesh);  //CHG
   initializephitoE(i_EM.phi,i_EM.E,i_mesh); //CHG

   if(i_mesh.Bloc==0) initializeBfield(i_EM.B,i_mesh.cmesh,i_mesh.meshdims,i_mesh.vecdims,i_svar.flag_CL); // cell centered
   else if(i_mesh.Bloc==1) initializeBfield(i_EM.B,i_mesh.pmesh,i_mesh.meshdims,i_mesh.vecdims,i_svar.flag_CL); // point (edge) centered

   if(i_mesh.meshdims==1) i_mesh.mesharea(i_EM);  //Set Cross-sectional area   CHG

   if(i_svar.rst>=1) initializeRestart(i_wrt,i_svar.totalTime,i_svar.nsp,i_svar.rst);  //Initialize starting from restart file

   for(i=0;i<nsp;i++)     // Iterate through different species
   {
      tmp_prt.pwght = i_svar.pwght; 
      tmp_cnt.pwght = i_svar.pwght;

#if SIMTYPE==3
      seedSingleParticle(tmp_prt,i_mesh,i_svar,i);
#else
      if(i_svar.rst==0) initializeParticles(tmp_prt,i_mesh,i_svar,i);  //CHG
      else initializeParticles(tmp_prt,i_mesh,i_svar,i,i_svar.rst);   
#endif

      if(i_mesh.philoc==1) initializeContinuum(tmp_cnt,tmp_prt,i_mesh.pmesh,i_mesh.meshdims,i_mesh.vecdims);//CHG
      else initializeContinuum(tmp_cnt,tmp_prt,i_mesh.cmesh,i_mesh.meshdims,i_mesh.vecdims);//CHG

      i_prt.push_back(tmp_prt);  
      i_cnt.push_back(tmp_cnt);  
      slvr.clearallParticles(tmp_prt);
      slvr.clearallFlow(tmp_cnt);
      i_cnt[i].nsp = nsp;
      i_prt[i].pid = i;

      for(j=0;j<i_mesh.cmesh.size();j++) 
      {
         temp = rand();
         temp = temp/RAND_MAX;
         i_prt[i].cseedcount.push_back(temp);
      }

      i_prt[i].fseedcount = 0.0;
   }

   if(i_mesh.philoc==1) initializeMPIvars(i_mpiv,i_mesh.pmesh,i_mesh.meshdims,i_mesh.vecdims,nsp);
   else initializeMPIvars(i_mpiv,i_mesh.cmesh,i_mesh.meshdims,i_mesh.vecdims,nsp);

   //writeOutput wrt(0,nsp);

   //for(j=0;j<nsp;j++) wrt.writeParticles(i_prt[j],i_EM,i_mesh,0.0); //  0


   if(i_svar.nct>0)   //..Initialize collision parameters 
   {  
      initializeNeutralBackground(tmp_cnt,i_mesh);
      i_cnt.push_back(tmp_cnt);
      slvr.clearallFlow(tmp_cnt);

      initializeCollisions(i_prt,i_mesh);
   }


   for(i=0;i<nsp;i++) i_mesh.initconnectpartandmesh(i_prt[i]); //Connect particles and mesh

   if(i_svar.rst>0 && numprocs>1) //MPI
   {
      for(j=0;j<i_svar.nsp;j++) slvr.redistributeparticles(i_prt[j],i_mesh.vecdims,i_mesh.meshdims);  //MPI
   }

   tmp_bound.initializeBoundary(i_mesh,i_bdv,i_prt); //Initialize boundary
   if(i_svar.rst>0)  tmp_bound.setParticleCounter(i_bdv,i_prt,i_svar.rst);   //Set particle count from restart file

   if(i_svar.nspcl>0) tmp_bound.initializeSpecialRegions(i_mesh,i_spclv,i_prt); //Initialize regions

   for(i=0;i<nsp;i++) i_mesh.connectpartandmesh(i_prt[i]); //Connect particles and mesh

   //for(i=0;i<nsp;i++) slvr.findvdf(i_prt[i],i_mesh,-1);
   //for(j=0;j<nsp;j++) wrt.writeParticles(i_prt[j],i_EM,i_mesh,0.0); //  1

   if(i_svar.q1dflag==1)
   {
      initializeQ1D(i_prt,i_EM,i_mesh,i_bdv);        //Initialize Quasi-1D 
   }

   //for(j=0;j<nsp;j++) wrt.writeParticles(i_prt[j],i_EM,i_mesh,0.0); //  2

   //for(j=0;j<nsp;j++) wrt.writePICField(i_mesh, i_cnt[j], i_EM,0.0); // 0
   MPI_Barrier(MPI_COMM_WORLD);
   slvr.weighContinuumPIC(i_prt, i_cnt, i_mesh,i_mpiv);
   //for(j=0;j<nsp;j++) wrt.writePICField(i_mesh, i_cnt[j], i_EM,0.0); // 0
   MPI_Barrier(MPI_COMM_WORLD);
   tmp_bound.applyContinuumBoundaryConditions(i_mesh,i_EM,i_cnt,i_bdv);
   //for(j=0;j<nsp;j++) wrt.writePICField(i_mesh, i_cnt[j], i_EM,0.0); // 0
   MPI_Barrier(MPI_COMM_WORLD);


   //slvr.thomas1Dpoisson(i_mesh,i_EM.phi,i_icont.N,i_econt.N,i_icont.charge,i_econt.charge);
#if SIMTYPE==0  || SIMTYPE==3
   slvr.poisson1D(i_mesh,i_EM.phi,i_cnt,i_bdv,i_prt,slvr.deltaT,0.0);   //..MB  Comment out for
#endif

   slvr.phitoE(i_EM.phi,i_EM.E,i_mesh);
   tmp_bound.applyEfieldBoundary(i_mesh,i_EM,i_bdv); //..MB Comment out for

   slvr.deltaT = i_svar.dt; 

   if(i_svar.rst==0) for(i=0;i<nsp;i++) slvr.updatePartVel(i_prt[i], i_EM, i_mesh, -i_mesh.mvscheme); //Initial half back step

   //std::cout << "\nAfter Backstep:\n";

   //for(j=0;j<nsp;j++) wrt.writeParticles(i_prt[j],i_EM,i_mesh,0.0); //  2
   //for(j=0;j<nsp;j++) wrt.writePICField(i_mesh, i_cnt[j], i_EM,0.0); // 0

   //for(i=0;i<nsp;i++) for(j=0;j<i_prt[i].pos.size();j++)  std::cout << i_prt[i].vel[j*i_mesh.vecdims] << "\t";
   //writeOutput wrt(0);  
   //wrt.writeParticles(i_prt[0],i_mesh.meshdims,i_mesh.vecdims);
   //wrt.writePICField(i_mesh, i_cnt[0], i_EM);

   //for(j=0;j<nsp;j++) wrt.writeParticles(i_prt[j],i_EM,i_mesh,0.0);

   std::cout << "\n\n........Finished Initializing Particles.........\n\n";

}

//....Initialize Particles (loading of particles)....//

void initializeSolver::initializeParticles(particles &i_part,mesh i_msh,solverVars i_svar,int spec)
{

   //.....Maxwellian set based on "Loading and Injection of Maxwellian Distributions in Particle Simulations" by Cartwright, Verboncoeur, and Birdsall....//

   int i,j,k,l;
   int npc;
   int numprocs,procid; //MPI
   long double temp1,temp2,temp3,tempen;
   double vupper,vlower,vtherm,i_dens,area,dnpc;
   bool cflag;
   std::string i_ddist,i_thdist;
   double i_pert;
   double tol=0.0e-10;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   i_part.spflag=spec;
   i_part.nsp = nsp;
   i_part.name = sname[spec];
   i_part.pwght = i_svar.pwght; 
   i_part.pnct = 0;

   std::vector<int> neigh;
   std::vector<double> i_pos;
   std::vector<double> i_vel;
   std::vector<std::string> i_initvel;
   std::vector<double> dx;
   std::vector<double> i_R,i_Rtot;

   for(i=0;i<i_msh.meshdims;i++)  i_pos.push_back(0.0);
   for(i=0;i<i_msh.meshdims;i++)  dx.push_back(0.0);
   for(i=0;i<i_msh.vecdims;i++)  i_vel.push_back(0.0);
   for(i=0;i<i_msh.vecdims;i++)  i_initvel.push_back("0.0");

   i_part.gp = 0;

   area = i_msh.areain;

   solver i_slvr;
   mathFunctions i_mth;

   i_part.mass = smass[spec];
   i_part.wmass = i_part.pwght*smass[spec];

   i_part.charge = scharge[spec];
   i_part.wcharge = i_part.pwght*scharge[spec];

   if(i_msh.mvscheme == 3)  i_part.mag = smag[spec]; // Set if Q1D
   else  i_part.mag = 1; // Set otherwise

   //i_part.cycIT = scycIT[spec];

   i_ddist = dens_dist[spec];
   i_thdist = therm_dist[spec];
   i_pert = dens_pert[spec];

   //....Seeding particles in space and velocity....//

   if(i_msh.pdist==0)  //Whole domain seeding
   {
      for(i=0;i<i_msh.vecdims;i++) i_initvel[i] = init_U[spec*i_msh.vecdims+i]; 

      std::cout << "\n\nInitializing " << sname[spec] << " particles over whole domain...\n";

      for(j=0;j<i_msh.meshdims;j++) i_pos[j] = 1.0;
      for(j=0;j<i_msh.meshdims;j++) dx[j] =  i_msh.meshlength[j];
      //area = 1.0;
      for(j=0;j<i_msh.meshdims;j++) area =  area*dx[j];

      i_dens = eqnparser(init_dens[spec],i_pos);
      for(j=0;j<i_msh.vecdims;j++) i_vel[j] = eqnparser(i_initvel[j],i_pos);
      vtherm = sqrt(eqnparser(init_Temp[spec],i_pos)*kb/(i_part.mass));

      //std::cout << std::endl << vtherm << std::endl;
      //std::cout << std::endl << i_vel[0] << "\t" << i_vel[1] << "\t" << i_vel[2] << std::endl;

      //npc = i_dens*area/i_part.pwght;   //..Using the inflow area
      //dnpc = i_dens*area/i_part.pwght;
      npc = (i_dens*area/i_part.pwght)/numprocs;   //..Using the inflow area  //MPI
      dnpc = (i_dens*area/i_part.pwght)/numprocs; //MPI

      if(npc==0) npc = 1;
      if((dnpc-npc)>0.5) npc = npc + 1;

      for(j=0;j<npc;j++) i_R.push_back((j+0.5)/npc);
      for(j=0;j<i_msh.vecdims*npc;j++) i_Rtot.push_back(0.0);

      //auto seed = time(NULL);
      auto seed = (time(NULL) + procid); //MPI

      for(j=0;j<i_msh.vecdims;j++)
      {
         //seed = (j+1)*time(NULL);
         seed = (j+1)*time(NULL) + procid; //MPI
         //std::random_shuffle(i_R.begin(),i_R.end());
         std::shuffle(i_R.begin(),i_R.end(),std::default_random_engine(seed));
         for(k=0;k<npc;k++) i_Rtot[i_msh.vecdims*k+j] = i_R[k];
      }

      vupper = 5.0;  //Upper velocity limit for distribution
      vlower = -5.0;  //Lower velocity limit for distribution

      temp1 = erf(1.5);
      temp1 = i_mth.erfinv(temp1);

      for(j=0;j<npc;j++)
      {
         if(i_ddist == "RANDOM") //Random spatial
         {
            for(k=0;k<i_msh.meshdims;k++)
            {
               temp1 =  rand();
               temp1 =  dx[k]*temp1/(RAND_MAX)+i_msh.meshstart[0];
#if SIMTYPE==1
               i_part.pos.push_back(0.0); // MB
#else
               i_part.pos.push_back(temp1);
#endif
            }
         }
         else if(i_ddist == "UNIFORM" || i_ddist == "PERTURBATION") //Uniform or perturbed spatial
         {
            for(k=0;k<i_msh.meshdims;k++)
            {
               temp1 =  (j+0.5)*(dx[k]/npc);
               temp1 =  temp1+i_pert*cos(2.0*pi*temp1/i_msh.meshlength[k])+tol+i_msh.meshstart[0];
#if SIMTYPE==1
               i_part.pos.push_back(0.0); // MB
#else
               i_part.pos.push_back(temp1);
#endif
            }
         }

         tempen=0.0;

         if(i_thdist == "MAXWELLIAN") //Maxwellian
         {
            for(k=0;k<i_msh.vecdims;k++)
            {
               temp1=rand();
               temp1=temp1/RAND_MAX;
               temp1=i_mth.erfinv(temp1*erf(vupper)+(1.0-temp1)*erf(vlower));
               temp1 = sqrt(2.0)*vtherm*temp1+i_vel[k];
               i_part.vel.push_back(temp1);
               tempen= tempen + (temp1)*(temp1);
            }
            i_part.en.push_back((0.5)*tempen*i_part.mass);   //Energy of individual particle, not group
         }
         else if(i_thdist == "QUIET") //Quiet Maxwellian
         {
            for(k=0;k<i_msh.vecdims;k++)
            {
               temp1=i_Rtot[j*i_msh.vecdims+k];
               temp1=i_mth.erfinv(temp1*erf(vupper)+(1.0-temp1)*erf(vlower));
               i_part.vel.push_back(sqrt(2.0)*vtherm*temp1+i_vel[k]);
               tempen= tempen + (sqrt(2.0)*vtherm*temp1+i_vel[k])*(sqrt(2.0)*vtherm*temp1+i_vel[k]);
            }
            i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle, not group
         }
         else if(i_thdist == "RANDOM") //RANDOM
         {
            for(k=0;k<i_msh.vecdims;k++)
            {
               temp1=rand()-RAND_MAX*0.5;
               temp1=2.0*temp1/RAND_MAX;
               if(k==0) i_part.vel.push_back(sqrt(2.0)*vtherm*temp1+i_vel[k]);
               else i_part.vel.push_back(sqrt(2.0)*sqrt(2.0)*vtherm*temp1+i_vel[k]);
               //i_part.vel.push_back(sqrt(2.0)*vtherm*temp1+i_vel[k]);
               tempen= tempen + (sqrt(2.0)*vtherm*temp1+i_vel[k])*(sqrt(2.0)*vtherm*temp1+i_vel[k]);
            }
            i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle, not group
         }
         else if(i_thdist == "UNIFORM") //Uniform
         {
            for(k=0;k<i_msh.vecdims;k++)
            {
               temp1=(j+0.5-0.5*npc)*2.0*sqrt(2.0)*vtherm/npc;
               i_part.vel.push_back(temp1+i_vel[k]);
               tempen= tempen + (temp1+i_vel[k])*(temp1+i_vel[k]);
            }
            i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle, not group
         }

      }
      i_R.clear();
      i_Rtot.clear();
   }
   else if(i_msh.pdist==1)   //Distribute particles cell by cell
   {
      for(i=0;i<i_msh.vecdims;i++) i_initvel[i] = init_U[spec*i_msh.vecdims+i]; 

      std::cout << "\n\nInitializing " << sname[spec] << " particles cell by cell...\n";

      for(i=0;i<(i_msh.cmesh.size()/i_msh.meshdims);i++)
      {
         cflag=true;

         for(j=0;j<i_msh.meshdims;j++)
         {
            if(i_msh.cmesh[i_msh.meshdims*i+j]<i_msh.meshbd[2*j] || i_msh.cmesh[i_msh.meshdims*i+j]>i_msh.meshbd[2*j+1]) cflag = false;
         }

         if(cflag == true)
         {
            for(j=0;j<i_msh.meshdims;j++) i_pos[j] = i_msh.cmesh[i_msh.meshdims*i+j];
            neigh = i_msh.pneighofc(i);
            for(j=0;j<i_msh.meshdims;j++) dx[j] =  fabs(i_msh.pmesh[i_msh.meshdims*neigh[2*j]+j]-i_msh.pmesh[i_msh.meshdims*neigh[2*j+1]+j]);
            area = i_msh.carea[i];
            for(j=0;j<i_msh.meshdims;j++) area =  area*dx[j];

            i_dens = eqnparser(init_dens[spec],i_pos);
            for(j=0;j<i_msh.vecdims;j++) i_vel[j] = eqnparser(i_initvel[j],i_pos);
            vtherm = sqrt(eqnparser(init_Temp[spec],i_pos)*kb/(i_part.mass));

            //npc = i_dens*area/i_part.pwght;
            //dnpc = i_dens*area/i_part.pwght;
            npc = (i_dens*area/i_part.pwght)/numprocs;  //MPI
            dnpc = (i_dens*area/i_part.pwght)/numprocs; //MPI


            if(npc==0) npc = 0;
            if((dnpc-npc)>0.5) npc = npc + 1;

            //std::cout << std::endl << i << "\t" << npc ;

            for(j=0;j<npc;j++) i_R.push_back((j+0.5)/npc);
            for(j=0;j<i_msh.vecdims*npc;j++) i_Rtot.push_back(0.0);

            //auto seed = time(NULL);
            auto seed = time(NULL)+procid;  //MPI
            for(j=0;j<i_msh.vecdims;j++)
            {
               //seed = (j+1)*time(NULL);
               seed = (j+1)*time(NULL)+procid; //MPI
               //std::random_shuffle(i_R.begin(),i_R.end());
               std::shuffle(i_R.begin(),i_R.end(),std::default_random_engine(seed));
               for(k=0;k<npc;k++) i_Rtot[i_msh.vecdims*k+j] = i_R[k];
            }

            vupper = 5.0;
            vlower = -5.0;

            temp1 = erf(1.5);
            temp1 = i_mth.erfinv(temp1);

            for(j=0;j<npc;j++)
            {
               if(i_ddist == "RANDOM") //Random spatial
               {
                  for(k=0;k<i_msh.meshdims;k++)
                  {
                     temp1 =  rand();
                     temp1 =  dx[k]*temp1/(RAND_MAX)+i_msh.pmesh[i_msh.meshdims*neigh[0]+k];
                     i_part.pos.push_back(temp1);
                  }
               }
               else if(i_ddist == "UNIFORM" || i_ddist == "PERTURBATION") //Uniform or perturbed spatial
               {
                  for(k=0;k<i_msh.meshdims;k++)
                  {
                     temp1 =  (j+0.5)*(dx[k]/npc)+i_msh.pmesh[i_msh.meshdims*neigh[0]+k];
                     temp1 =  temp1+i_pert*cos(2.0*pi*temp1/i_msh.meshlength[k]+tol);
                     i_part.pos.push_back(temp1);
                  }
               }

               tempen=0.0;

               if(i_thdist == "MAXWELLIAN") //Maxwellian velocity
               {
                  for(k=0;k<i_msh.vecdims;k++)
                  {
                     temp1=rand();
                     temp1=temp1/RAND_MAX;
                     temp1=i_mth.erfinv(temp1*erf(vupper)+(1-temp1)*erf(vlower));
                     i_part.vel.push_back(sqrt(2.0)*vtherm*temp1+i_vel[k]);
                     tempen= tempen + (sqrt(2.0)*vtherm*temp1+i_vel[k])*(sqrt(2.0)*vtherm*temp1+i_vel[k]);
                  }
                  i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle, not group
               }
               else if(i_thdist == "QUIET") //Quiet Maxwellian Velocity
               {
                  for(k=0;k<i_msh.vecdims;k++)
                  {
                     temp1=i_Rtot[j*i_msh.vecdims+k];
                     temp1=i_mth.erfinv(temp1*erf(vupper)+(1-temp1)*erf(vlower));
                     i_part.vel.push_back(sqrt(2)*vtherm*temp1+i_vel[k]);
                     tempen= tempen + (sqrt(2)*vtherm*temp1+i_vel[k])*(sqrt(2)*vtherm*temp1+i_vel[k]);
                  }
                  i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle not group
               }
               else if(i_thdist == "RANDOM") // RANDOM
               {
                  for(k=0;k<i_msh.vecdims;k++)
                  {
                     temp1=rand()-RAND_MAX*0.5;
                     temp1=2.0*temp1/RAND_MAX;
                     i_part.vel.push_back(sqrt(2.0)*vtherm*temp1+i_vel[k]);
                     tempen= tempen + (sqrt(2.0)*vtherm*temp1+i_vel[k])*(sqrt(2.0)*vtherm*temp1+i_vel[k]);
                  }
                  i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle, not group
               }
               else if(i_thdist == "UNIFORM") //Uniform
               {
                  for(k=0;k<i_msh.vecdims;k++)
                  {
                     temp1=(j+0.5-npc*0.5)*sqrt(2.0)*2.0*vtherm/npc;
                     i_part.vel.push_back(temp1+i_vel[k]);
                     tempen= tempen + (temp1+i_vel[k])*(temp1+i_vel[k]);
                  }
                  i_part.en.push_back((0.5)*tempen*i_part.mass);  //Energy of individual particle, not group
               }

            }
            i_R.clear();
            i_Rtot.clear();
         }
      }
   }

   std::cout << "Total " << sname[spec] << " particles:  " <<  i_part.pos.size()/i_msh.meshdims << std::endl;
}


//....Initialize particles from file....//

void initializeSolver::initializeParticles(particles &i_part,mesh i_msh,solverVars i_svar,int spec,int rst)
{

   int i,j,k,l;
   int numprocs,procid; //MPI
   long double tempdouble;
   bool cflag;
   std::string i_ddist,i_thdist;
   double i_pert;
   double tol=0.0e-10;
   int i_total_particles;
   int pindex;
   double entemp;

   std::vector<int> w_neigh;
   std::vector<double> w_weight;
   std::vector<double> enphi;
   std::vector<double> enE,intE;
   std::vector<double> w_pos;

   std::stringstream fname;
   std::string filename;

   MPI_Comm_size(MPI_COMM_WORLD,&numprocs);   //MPI
   MPI_Comm_rank(MPI_COMM_WORLD,&procid);  //MPI

   i_part.spflag=spec;
   i_part.nsp = nsp;
   i_part.name = sname[spec];
   i_part.pwght = i_svar.pwght; 
   i_part.pnct = 0;

   i_part.gp = 0;

   solver i_slvr;
   mathFunctions i_mth;

   i_part.mass = smass[spec];
   i_part.wmass = i_part.pwght*smass[spec];

   i_part.charge = scharge[spec];
   i_part.wcharge = i_part.pwght*scharge[spec];

   if(i_msh.mvscheme == 3)  i_part.mag = smag[spec]; // Set if Q1D
   else  i_part.mag = 1; // Set otherwise

   i_ddist = dens_dist[spec];
   i_thdist = therm_dist[spec];
   i_pert = dens_pert[spec];

   std::cout << "\n\tReading Particles...";
   fname << i_part.name.c_str() << "Particles" << rst << ".out";

   filename = fname.str();
   std::ifstream rdfile(filename.c_str(),std::ios_base::app | std::ios_base::out);  //MPI 

   //....NOTE:  CHANGE TO .out AND NEED TO DELETE TECPLOT HEADERS AND REPLACE WITH JUST NUMBER OF PARTICLES

   rdfile >> i_total_particles;

   if(procid==0)
   {
      if(i_msh.meshdims==1)
      {
         for(i = 0; i<i_total_particles; i++)
         {
            rdfile >> tempdouble;
            i_part.pos.push_back(tempdouble);
            rdfile >> tempdouble;
            for(j=0;j<i_msh.vecdims; j++) 
            {  
               rdfile >> tempdouble;
               i_part.vel.push_back(tempdouble);
            }
            rdfile >> tempdouble;
            i_part.en.push_back(tempdouble);
         }
      }
      else if(i_msh.meshdims==2)
      {
         for(i = 0; i<i_total_particles; i++)
         {
            rdfile >> tempdouble;
            i_part.pos.push_back(tempdouble);
            rdfile >> tempdouble;
            i_part.pos.push_back(tempdouble);
            for(j=0;j<i_msh.vecdims; j++) 
            {  
               rdfile >> tempdouble;
               i_part.vel.push_back(tempdouble);
            }
            rdfile >> tempdouble;
            i_part.en.push_back(tempdouble);
         }
      }
   }

   std::cout << "Total " << sname[spec] << " particles:  " <<  i_part.pos.size() << std::endl;
}

//....Initialize magnetic field....//

void initializeSolver::initializeBfield(std::vector<double> & i_EM_B, std::vector<double> i_mesh,int i_mesh_dims,int i_vec_dims,bool i_flag_CL)
{
   int i,j,k;
   int i_tot_cells;
   double temp1,temp2;
   std::vector<double> i_pos,i_B;

   for(i=0;i<i_vec_dims;i++) i_B.push_back(0.0);
   for(i=0;i<i_mesh_dims;i++) i_pos.push_back(0.0);

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   if(i_flag_CL==true) //Magnetic field for current loop
   {
      for(i=0; i<i_tot_cells; i++)
      {
         for(j=0;j<i_mesh_dims;j++) i_pos[j] = i_mesh[i_mesh_dims*i+j]; 
         i_B = loop_cartesian(i_pos,i_mesh_dims,i_vec_dims);
         for(j=0;j<i_vec_dims;j++) i_EM_B.push_back(i_B[j]); 
      }
   }
   else
   {
      for(i=0; i<i_tot_cells; i++)  //Magnetic field for equation
      {
         for(j=0; j<i_mesh_dims; j++)  i_pos[j] = i_mesh[i_mesh_dims*i+j];
         for(j=0; j<i_vec_dims; j++)   i_EM_B.push_back(eqnparser(init_B[j],i_pos));
      }
   }
}

//....Initialize electric field....//

void initializeSolver::initializeEfield(std::vector<double> &i_EM_E, std::vector<double> i_mesh,int i_mesh_dims, int i_vec_dims)
{
   int i,j,k;
   int i_tot_cells;
   double temp1,temp2;
   std::vector<double> i_pos;

   for(i=0;i<i_mesh_dims;i++) i_pos.push_back(0.0);

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   for(i=0; i<i_tot_cells; i++) //Initialize with equation
   {
      for(j=0;j<i_mesh_dims;j++) i_pos[j] = i_mesh[i_mesh_dims*i+j]; 
      for(j=0; j<i_vec_dims; j++)   i_EM_E.push_back(eqnparser(init_E[j],i_pos));
   }
}

//....Initialize electric potential....//

void initializeSolver::initializephi(std::vector<double> &i_EM_phi, const mesh &i_msh)
{
   int i,j,k;
   int i_tot_cells;
   double temp1,temp2;
   std::vector<double> i_pos;

   for(i=0;i<i_msh.meshdims;i++) i_pos.push_back(0.0);

   if(i_msh.philoc==0)
   {
      i_tot_cells = i_msh.cmesh.size()/i_msh.meshdims;

      for(i=0; i<i_tot_cells; i++) //Initialize with equation
      {
         for(j=0;j<i_msh.meshdims;j++) i_pos[j] = i_msh.cmesh[i_msh.meshdims*i+j]; 
         i_EM_phi.push_back(eqnparser(init_phi,i_pos));
      }
   }
   else if(i_msh.philoc==1)
   {
      i_tot_cells = i_msh.pmesh.size()/i_msh.meshdims;

      for(i=0; i<i_tot_cells; i++) //Initialize with equation
      {
         for(j=0;j<i_msh.meshdims;j++) i_pos[j] = i_msh.pmesh[i_msh.meshdims*i+j]; 
         i_EM_phi.push_back(eqnparser(init_phi,i_pos));
      }
   }
}

//....Initially relate electric potential and electric field....//

void initializeSolver::initializephitoE(std::vector<double> &i_EM_phi, std::vector<double> &i_EM_E, const mesh &i_mesh)
{
   int i,j,k;
   int num_cells;
   double temp1,temp2;

   mathFunctions mth;

   if(i_mesh.meshdims==1)
   {
      std::cout << "\t\nCalculating E-Field....";

      if(i_mesh.philoc==0)
      {
         num_cells = i_mesh.pmesh.size();

         for(i=0;i<num_cells;i++)
         {
            //i_EM_E.push_back(-mth.ddx1Dpwc(i_mesh,i_EM_phi,i));
            i_EM_E.push_back(-mth.ddx1Dpwc(i_mesh,i_mesh.deltax,i_EM_phi,i)); //EFF
            i_EM_E.push_back(0.0);
            i_EM_E.push_back(0.0);
         }
      }
      else if(i_mesh.philoc==1)
      {
         num_cells = i_mesh.pmesh.size();

         i_EM_E.push_back(0.0);
         i_EM_E.push_back(0.0);
         i_EM_E.push_back(0.0);
         for(i=1;i<(num_cells-1);i++)
         {
            i_EM_E.push_back(-mth.ddx1Dpwp(i_mesh,i_EM_phi,i));
            i_EM_E.push_back(0.0);
            i_EM_E.push_back(0.0);
         }
         i_EM_E.push_back(0.0);
         i_EM_E.push_back(0.0);
         i_EM_E.push_back(0.0);
      }
   } 

}

//....Initialize quasi 1D mode....//

void initializeSolver::initializeQ1D(std::vector<particles> &i_part, fields &i_EM, mesh &i_msh, std::vector<boundvars> &i_bdv)
{
   int i,j,k;
   int n;
   double moq;
   double temp1,temp2;
   std::vector<double> i_pos,int_B;

   mathFunctions i_mth;
   solver i_slvr;

   if(i_msh.Bloc==0)  //Initialize gradB Coefficient
   {
      n = i_msh.cintcells.size();

      for(i=0;i<i_bdv[0].cboundnum;i++) i_EM.gradBcoef.push_back(-(0.5)*i_mth.ddx1Dcwc(i_msh,i_EM.B,0,i_bdv[0].cboundnum));
      for(i=0; i<n; i++)  i_EM.gradBcoef.push_back((-0.5)*i_mth.ddx1Dcwc(i_msh,i_EM.B,0,i+i_bdv[0].cboundnum));
      for(i=0;i<i_bdv[1].cboundnum;i++) i_EM.gradBcoef.push_back(-(0.5)*i_mth.ddx1Dcwc(i_msh,i_EM.B,0,n-1+i_bdv[0].cboundnum));
   }
   else if(i_msh.Bloc==1)
   {
      n = i_msh.pintcells.size();

      for(i=0;i<i_bdv[0].pboundnum;i++) i_EM.gradBcoef.push_back(-(0.5)*i_mth.ddx1Dpwp(i_msh,i_EM.B,0,i_bdv[0].pboundnum));  //CHK
      for(i=0; i<n; i++)   i_EM.gradBcoef.push_back((-0.5)*i_mth.ddx1Dpwp(i_msh,i_EM.B,0,i+i_bdv[0].pboundnum));
      for(i=0;i<i_bdv[1].pboundnum;i++) i_EM.gradBcoef.push_back(-(0.5)*i_mth.ddx1Dpwp(i_msh,i_EM.B,0,n-1+i_bdv[0].pboundnum));

   }

   /*std::cout << std::endl << std::endl << "GradB:     "  << n << "\t" << i_EM.gradBcoef.size() << std::endl;
     for(i=0;i<i_EM.gradBcoef.size();i++) std::cout << i_EM.gradBcoef[i] << "\t";
     std::cout << std::endl << std::endl;
     */
   //for(i=0;i<i_msh.vecdims;i++) int_B.push_back(0.0);

   /*for(i=0;i<i_part[0].nsp;i++) //Set v_perp from v_r and v_theta
     {
     moq = fabs(i_part[i].mass/i_part[i].charge);

     if(i_part[i].mag==1)
     {
     for(j=0;j<i_part[i].pos.size();j++)  
     {
     i_part[i].vel[j*i_msh.vecdims+2] = sqrt(i_part[i].vel[j*i_msh.vecdims+2]*i_part[i].vel[j*i_msh.vecdims+2] + i_part[i].vel[j*i_msh.vecdims+1]*i_part[i].vel[j*i_msh.vecdims+1]);
     i_part[i].vel[j*i_msh.vecdims+1] = 0.0;
     }
     }
     }*/

}

//....Initialize fluid properties....//

void initializeSolver::initializeFluid(flow & i_cont, std::vector<double> i_mesh,int i_mesh_dims, int i_vec_dims)
{
   int i,j,k;
   int i_tot_cells;
   double Usq;
   std::vector<double> i_pos;

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   for(i=0;i<i_mesh_dims;i++) i_pos.push_back(0.0);

   for(i=0; i<i_tot_cells; i++)
   {
      for(j=0;j<i_mesh_dims;j++) i_pos[j] = i_mesh[i_mesh_dims*i+j]; 
      i_cont.rho.push_back(eqnparser(init_rho,i_pos));
      i_cont.p.push_back(eqnparser(init_p,i_pos));
   }
   for(i=0; i<i_tot_cells; i++)
   {
      for(j=0;j<i_mesh_dims;j++) i_pos[j] = i_mesh[i_mesh_dims*i+j]; 
      for(j=0; j<i_vec_dims; j++)  i_cont.U.push_back(eqnparser(init_U[j],i_pos));
      //   std::cout << "\n" << i_pos[0] << "\t" << i_cont.U[i] << "\n";
   }

   for(i=0;i<i_tot_cells; i++)
   {
      Usq = 0;
      for(j=0;j<i_vec_dims;j++) Usq = Usq + i_cont.U[i_vec_dims*i+j]*i_cont.U[i_vec_dims*i+j];   
      i_cont.En.push_back(i_cont.p[i]/(i_cont.gam-1.0)+i_cont.rho[i]*Usq/2.0);
   }


}


//....Initialize flow properties....//

void initializeSolver::initializeFlow(flow & i_cont, std::vector<double> i_mesh,int i_mesh_dims, int i_vec_dims)
{
   int i,j,k;
   int i_tot_cells;

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   for(i=0; i<i_tot_cells; i++)
   {
      i_cont.N.push_back(0.0);
      i_cont.T.push_back(0.0);
   }
   for(i=0; i<i_tot_cells; i++)
   {
      for(j=0; j<i_vec_dims; j++)  i_cont.U.push_back(0.0);
   }
}


//....Initialize flow properties for PIC solver....//

void initializeSolver::initializeFlowPIC(flow & i_cont, std::vector<double> i_mesh,int i_mesh_dims, int i_vec_dims)
{
   int i,j,k;
   int i_tot_cells;

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   for(i=0; i<i_tot_cells; i++)
   {
      i_cont.rho.push_back(0.0);
   }
}


//....Initialize continuum properties....//

void initializeSolver::initializeContinuum(contnm & i_cont, const particles & i_part, std::vector<double> i_mesh,int i_mesh_dims, int i_vec_dims)
{
   int i,j,k;
   int i_tot_cells;
   std::vector<double> i_pos;

   i_cont.spflag = i_part.spflag;
   i_cont.name = i_part.name;
   i_cont.mass = i_part.mass;
   i_cont.wmass = i_part.pwght*i_part.mass;
   i_cont.charge = i_part.charge;
   i_cont.wcharge = i_part.pwght*i_part.charge;
   i_cont.pwght = i_part.pwght; 

   std::cout << "\nInitializing " << i_cont.name << " continuum....\n";

   for(i=0;i<i_mesh_dims;i++) i_pos.push_back(0.0);

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   for(i=0; i<i_tot_cells; i++)
   {
      for(j=0;j<i_mesh_dims;j++) i_pos[j] = i_mesh[i_mesh_dims*i+j]; 
      i_cont.N.push_back(0.0);
      i_cont.En.push_back(0.0);
      i_cont.T.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_cont.Tvec.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_cont.Qvec.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_cont.U.push_back(0.0);
      i_cont.rho_back.push_back(eqnparser(rhoback,i_pos));
   }

}

void initializeSolver::initializeMPIvars(MPIvars & i_mpiv, std::vector<double> i_mesh,int i_mesh_dims, int i_vec_dims, int i_nsp)
{
   int i,j,k;
   int i_tot_cells;

   i_tot_cells = i_mesh.size()/i_mesh_dims;

   for(i=0; i<i_tot_cells*i_nsp; i++)
   {
      i_mpiv.Ntotal.push_back(0.0);
      i_mpiv.N_all.push_back(0.0);
      i_mpiv.Entotal.push_back(0.0);
      i_mpiv.En_all.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_mpiv.Tvectotal.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_mpiv.Qvectotal.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_mpiv.Tvec_all.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_mpiv.Qvec_all.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_mpiv.Utotal.push_back(0.0);
      for(j=0;j<i_vec_dims;j++) i_mpiv.U_all.push_back(0.0);
   }

}


void initializeSolver::initializeNeutralBackground(contnm & i_cont, const mesh &i_mesh)
{
   int i,j,k;
   int i_tot_cells;
   std::vector<double> i_pos;

   std::cout << "\nInitializing Neutral Background....\n";

   i_cont.spflag = i_cont.nsp;
   i_cont.name = "neutrals";
   i_cont.mass = 1.0;
   i_cont.wmass = 1.0;
   i_cont.charge = 1.0;
   i_cont.wcharge = 1.0;
   i_cont.pwght = 1.0; 

   std::cout << "\nInitializing " << i_cont.name << " continuum....\n";

   if(i_mesh.philoc==0 )  i_tot_cells = i_mesh.cmesh.size()/i_mesh.meshdims;
   else if(i_mesh.philoc==1) i_tot_cells = i_mesh.pmesh.size()/i_mesh.meshdims;

   for(i=0;i<i_mesh.meshdims;i++) i_pos.push_back(0.0);

   for(i=0; i<i_tot_cells; i++)
   {
      if(i_mesh.philoc==0) for(j=0;j<i_mesh.meshdims;j++) i_pos[j] = i_mesh.cmesh[i_mesh.meshdims*i+j]; 
      if(i_mesh.philoc==1) for(j=0;j<i_mesh.meshdims;j++) i_pos[j] = i_mesh.pmesh[i_mesh.meshdims*i+j]; 
      i_cont.N.push_back(eqnparser(neutdens,i_pos));
      i_cont.En.push_back(0.0);
      i_cont.T.push_back(eqnparser(neuttemp,i_pos));
      for(j=0;j<i_mesh.vecdims;j++) i_cont.U.push_back(0.0);
      i_cont.ioncount.push_back(0.0);
   }

   if(scrsstype[0]=="ARGON")  i_cont.R = 208.0; 
   else if(scrsstype[0]=="XENON")  i_cont.R = 63.324; 
   else if(scrsstype[0]=="HELIUM")  i_cont.R = 2077.0; 

   i_cont.NL = i_cont.N[0];
   i_cont.TL = i_cont.T[0];

}




//....Equation parser used to initialize code....//

double initializeSolver::eqnparser(std::string expression_string, std::vector<double> i_point) 
{
   int i_size = i_point.size();
   double i_x, i_y, i_z;
   exprtk::symbol_table<double> symbol_table;
   double time = 0.0;   


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


   symbol_table.add_variable("t",time);

   exprtk::expression<double> expression;
   expression.register_symbol_table(symbol_table);

   exprtk::parser<double> parser;
   parser.compile(expression_string,expression);

   return expression.value();
}


//....Current loop calculations for cartesian mesh....//

std::vector<double> initializeSolver::loop_cartesian(std::vector<double> centroid,int i_mesh_dims, int i_vec_dims)
{
   double i_x, i_y, i_z, i_B0, i_r, i_theta;
   double i_alpha, i_beta, i_gamma, i_Q, i_kay, i_K, i_E;
   std::vector<double>  i_B(3);
   double  i_Br;	

   mathFunctions mth;  

   if(i_mesh_dims == 2) centroid.push_back(0.0);

   i_x = centroid[0] - loopcenter[0];
   i_y = centroid[1] - loopcenter[1];
   i_z = centroid[2] - loopcenter[2];

   i_B0= loopcurrent*mu0/(2.0*loopradius);

   i_r = sqrt(i_y*i_y + i_z*i_z);
   i_theta = atan2(i_y,i_z);

   i_alpha = i_r/loopradius;
   i_beta = i_x/loopradius;
   i_gamma = i_x/i_r;

   if (i_alpha == 1.0 && i_x == 0.0)  i_alpha = i_alpha + 1E-3;

   i_Q = (1.0+i_alpha)*(1.0+i_alpha)+i_beta*i_beta;
   i_kay = sqrt(4.0*i_alpha/i_Q);

   i_K = mth.el_first_comp(i_kay);
   i_E = mth.el_second_comp(i_kay);

   i_B[0] = (i_B0/(pi*sqrt(i_Q)))*(i_E*(1.0-i_alpha*i_alpha-i_beta*i_beta)/(i_Q-4.0*i_alpha)+i_K) ; 
   i_Br = i_B0*i_gamma/(pi*sqrt(i_Q))*(i_E*(1.0+i_alpha*i_alpha+i_beta*i_beta)/(i_Q-4.0*i_alpha)-i_K);

   if (i_r == 0.0)  i_Br = 0.0;

   i_B[1] =  i_Br*sin(i_theta);
   i_B[2] =  i_Br*cos(i_theta);

   return i_B;
}

//....Seed single particle (for debugging)....//

void initializeSolver::seedSingleParticle(particles &i_part, mesh i_msh,solverVars i_svar, int spec)
{
   int i,j,k;
   double tol = 1.0e-12;
   std::vector<std::string> i_initvel;
   std::string i_ddist;
   double i_pert,entemp;


   for(i=0;i<i_msh.vecdims;i++) i_initvel.push_back("0.0");

   i_part.spflag=spec;
   i_part.nsp = nsp;
   i_part.name = sname[spec];
   i_part.pwght = i_svar.pwght; 

   i_part.mass = smass[spec];//i_svar.atmwght*amu;
   i_part.wmass = i_part.pwght*smass[spec];

   i_part.charge = scharge[spec];
   i_part.wcharge = i_part.pwght*scharge[spec];

   i_ddist = dens_dist[spec];
   i_pert = dens_pert[spec];

   if(i_msh.mvscheme == 3)  i_part.mag = smag[spec]; // Set if Q1D
   else  i_part.mag = 1; // Set otherwise

   for(i=0;i<i_msh.vecdims;i++) i_initvel[i] = init_U[spec*i_msh.vecdims+i]; 

   std::cout << "\n\nInitializing " << sname[spec] << "...\n";

   std::cout << "\nSeeding Single Particle....\n\n";

   entemp = 0.0;

   //for(i=0;i<1;i++)   i_part.pos.push_back(0.5*i_msh.meshlength[i]+i_pert);
   for(i=0;i<1;i++)   i_part.pos.push_back(0.0+i_pert);
   for(i=0;i<i_msh.vecdims;i++)    i_part.vel.push_back(eqnparser(i_initvel[i],i_part.pos));
   for(i=0;i<i_msh.vecdims;i++)    entemp = entemp + i_part.vel[i]*i_part.vel[i];
   i_part.en.push_back(entemp*0.5*i_part.mass);
}


//......Clears special region variables.......////

void initializeSolver::clearallspclv(spclvars &i_spclv)
{
   i_spclv.spclrange.clear();
   i_spclv.spcltype.clear();
   i_spclv.spcldens.clear();
   i_spclv.spclddist.clear();
   i_spclv.spclT.clear();
   i_spclv.spclthdist.clear(); 
   i_spclv.spclU.clear(); 
   i_spclv.spclFlux.clear(); 
   i_spclv.spclrange.clear(); 
}


//.......Initialize Collision parameters........///

void initializeSolver::initializeCollisions(std::vector<particles> &i_prt, mesh &i_msh)
{
   int i,j,k;
   double temp;

   for(i=0;i<nct;i++)
   {
      for(j=0;j<i_prt[0].nsp;j++)
      {
         if(scllsnname[i]==i_prt[j].name)
         { 
            i_prt[j].pnct += 1;
            i_prt[j].cllsntype.push_back(scllsntype[i]);
            i_prt[j].elcount= 0;
            i_prt[j].exccount = 0;
            i_prt[j].ioncount = 0;
            i_prt[j].ncolcount = 0.5;
            i_prt[j].crsssctn.push_back(scrsssctn[i]);

            if(i_prt[j].cllsntype[i_prt[j].pnct-1]=="NEUTRAL") //...Setup neutral collisions
            {
               i_prt[j].crsstype.push_back(scrsstype[i]);
               if(i_prt[j].crsstype[i_prt[j].pnct-1] != "CONSTANT")
               {
                  readcrsstable(i_prt[j].cllsnenergy,i_prt[j].en_crsstable,i_prt[j].el_crsstable,i_prt[j].inel_crsstable, i_prt[j].ion_crsstable, i_prt[j].name,i_prt[j].crsstype[i_prt[j].pnct-1],i_prt[j].cllsnmass,i_prt[j].cllsnB);
                  i_prt[j].crsssctn_max = *std::max_element(i_prt[j].el_crsstable.begin(), i_prt[j].el_crsstable.end()) + *std::max_element(i_prt[j].inel_crsstable.begin(), i_prt[j].inel_crsstable.end()) + *std::max_element(i_prt[j].ion_crsstable.begin(), i_prt[j].ion_crsstable.end());
               }

               /*i_prt[j].invmfp_max = 0.0; 
                 for(k=0;k<i_prt[j].cllsnenergy.size();k++)
                 {
                 temp = sqrt(2.0*i_prt[j].cllsnenergy[k]/i_prt[j].mass)*(i_prt[j].el_crsstable[k]+i_prt[j].inel_crsstable[k]+i_prt[j].ion_crsstable[k]);
                 if(temp>i_prt[j].invmfp_max) i_prt[j].invmfp_max = temp;
                 }*/

            }


            std::cout << std::endl <<"crssctn_max" <<  i_prt[j].crsssctn_max;// << "inv_mfp " << i_prt[j].invmfp_max  << std::endl;

            /*srand(time(NULL));

              for(k=0;k<i_msh.cmesh.size();k++) //..Initial counters for cell by cell probability (not used) 
              {
              temp = rand();
              temp = temp/RAND_MAX;
              i_prt[j].cllsncount.push_back(temp);
              }*/
         }
      }
   } 
   std::cout << std::endl << i_prt[0].pnct << "\t" << i_prt[1].pnct << std::endl;

}


//..Read in table for neutral cross sections..//

void initializeSolver::readcrsstable(std::vector<double> &i_cllsnenergy, std::vector<double> &i_en_crsstable, std::vector<double> &i_el_crsstable, std::vector<double> &i_inel_crsstable, std::vector<double> &i_ion_crsstable, std::string i_pname, std::string i_crsstype, double &i_cllsnmass, double &i_cllsnB)
{
   int i,j,k;
   std::string temp;
   double tempdouble;
   std::stringstream fname;
   std::string filename;
   std::string::size_type sz;

   fname << "CrossSectionData/" << i_pname << "_" << i_crsstype << "_crosssection_data.txt"; 
   filename = fname.str();

   std::cout << std::endl << filename << std::endl;

   std::ifstream rdfile(filename);

   rdfile >> i_cllsnB;

   rdfile >> i_cllsnmass;

   for(i=0;i<3;i++)  //..Read in energies consumed by collisions
   {
      rdfile >> tempdouble;
      i_cllsnenergy.push_back(tempdouble);
   }

   rdfile >> temp;

   while(temp!="EOF")
   {
      tempdouble = std::stod(temp,&sz);
      i_en_crsstable.push_back(tempdouble);

      rdfile >> tempdouble;
      i_el_crsstable.push_back(tempdouble);

      rdfile >> tempdouble;
      i_inel_crsstable.push_back(tempdouble);

      rdfile >> tempdouble;
      i_ion_crsstable.push_back(tempdouble);

      rdfile >> temp;
   }

   /*for(i=0;i<i_en_crsstable.size();i++) 
     {
     std::cout << "\n" <<  i_en_crsstable[i] << std::endl;
     std::cout << "\n" <<  i_el_crsstable[i] << std::endl;
     std::cout << "\n" <<  i_inel_crsstable[i] << std::endl;
     std::cout << "\n" <<  i_ion_crsstable[i] << std::endl;
     }*/


}

//..Initialize restarting from files

void initializeSolver::initializeRestart(writeOutput &i_wrt, double &i_time, int i_nsp, int i_rst)
{
   int i,j;
   int i_iter;
   std::stringstream fname;
   std::string filename;

   std::string name;

   std::cout << "\n\tReading Restart File...." ;

   name = "restart";

   fname <<  name.c_str() << i_rst << ".out";

   filename = fname.str();

   std::ifstream rdfile(filename.c_str());

   rdfile >> i_iter >> i_time;

   std::cout << std::endl << i_time << std::endl;

   for(i=0;i<i_nsp;i++)
   {
      rdfile >> i_wrt.part_fnum[i];
      rdfile >> i_wrt.cont_fnum[i];
      rdfile >> i_wrt.field_fnum[i];
      rdfile >> i_wrt.vdf_fnum[i];
      rdfile >> i_wrt.ps_fnum[i];
   }

   rdfile >> i_wrt.totcont_fnum >> i_wrt.rst_fnum;

   rdfile.close();
}

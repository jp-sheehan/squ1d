#include <iostream>
#include <stdio.h> 
#include <stdlib.h> 
#include <ctype.h>  
#include <unistd.h> 
#include "solvertype.h"
#include "initializesolver.h"
#include "solver.h"
#include "boundary.h"
#include "solverdefs.h"

/*...............................
   _____ ____       _______  
  / ___// __ \__  _<  / __ \
  \__ \/ / / / / / / / / / /
 ___/ / /_/ / /_/ / / /_/ /
/____/\___\_\__,_/_/_____/


           .--'''''''''--.
        .'      .---.      '.
       /    .-----------.    \
      /        .-----.        \
      |       .-.   .-.       |
      |      /   \ /   \      |
       \    | .-. | .-. |    /
        '-._| | | | | | |_.-'
            | '-' | '-' |
             \___/ \___/
          _.-'  /   \  `-._
        .' _.--|     |--._ '.
        ' _...-|     |-..._ '
               |     |
               '.___.'
                 | |
                _| |_
               /\( )/\
              /  ` '  \
             | |     | |
             '-'     '-'
             | |     | |
             | |     | |
             | |-----| |
          .`/  |     | |/`.
          |    |     |    |
          '._.'| .-. |'._.'
                \ | /
                | | |
                | | |
                | | |
               /| | |\
             .'_| | |_`.
             `. | | | .'
          .    /  |  \    .
         /o`.-'  / \  `-.`o\
        /o  o\ .'   `. /o  o\
        `.___.'       `.___.'

................................*/


int main(int argc, char *argv[])
{
   int sflag;
   int cflag,rflag,temp;

   std::cout << "Howdy, Welcome to Frans' PIC Solver \n \n";

   MPI_Init(&argc, &argv);

   
   //..Set runtime flags

   /*while((temp=getopt(argc,argv, "mcs")) !=-1)
   {
     switch (temp)
     {
       case 'm':    //..Magnetic Mirror/Bottle
         rflag = 1;
         break;
       case 'c':
         rflag = 2; //..Collision Tests
         break;
       case 's':
         rflag = 3; //..Single Particle
         break;
       case '?':
         rflag = 0;
         break;
       default:
         rflag = 0;
     }
   }*/

   cflag = COMPILER;
   rflag = SIMTYPE;

   if(cflag==0) std::cout << "\nCompiled with g++, running with LAPACK...\n\n";
   else if(cflag==1) std::cout << "\nCompiled with Intel, running with MKL...\n\n";

 
   if(rflag==0) std::cout << "\nRunning full simulation...\n\n";
   else if(rflag==1) std::cout << "\nRunning as magnetic mirror simulation...\n\n";
   else if(rflag==2) std::cout << "\nRunning as collision test simulation...\n\n";
   else if(rflag==3) std::cout << "\nRunning as single particle simulation...\n\n";
 
	
   initializeSolver pr;
   sType st;
   
   sflag = pr.readsolver();

   if(sflag==0) st.electrostatic1DPIC();
   else if(sflag==2) st.electromagnetic2DPIC();
   else if(sflag==10) st.euler1D();   
   else if(sflag==11) st.euler2D();   

   MPI_Finalize();

   return 0;
}

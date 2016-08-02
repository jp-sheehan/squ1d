#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <chrono>
#include "variables.h"
#include "constants.h"
#include "mesh.h"
#include "mathfunctions.h"
#include "time.h"
#include <ctime>
#include <iostream>
#include <iomanip>
#include "omp.h"

#ifndef STYPE_H
#define STYPE_H

class sType
{
	private:
 
	public:

        void electrostatic1DPIC() ;
        /* The following have not been implemented
        void euler1D() ;
        void euler2D() ;
        void electromagnetic2DPIC() ;
        */
};
#endif

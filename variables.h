#ifndef VARIABLE_H
#define VARIABLE_H

struct solverVars
{
   double lscale,tscale, cfl, dt, atmwght,gam;  //..Scales,timesep,atomic weight, gamma
   double totalTime;   //..Total Simualtion Time
   int  iter,p_iter;  //..Iterations
   bool flag_CL;  //..Current loop flag
   int  tdflag;   //..Temporal scheme (CFD)
   int  sdflag;  //..Spatial scheme flag (CFD)
   double  pwght;  //..Weight of particles
   int  nsp;  //..Number of species
   int  nbound;  //..Number of boundaries
   int  mvscheme;  //..Particle Mover scheme
   int  outavg,outcont,outpart,outvdf,outrst,outfinal,outfinalskip; //..Output info
   int  outvdf_flag;  //..Type of vdf output
   int  q1dflag;  //..Quasi1D flag
   int  nspcl;  //..Number of special regions
   int  nct;  //..Number of collision types
   int  rst;  //..Restart flag
   int  coul;  //..Coulomb Collision flag

};

struct particles
{
        std::string name;  //..Species name
	std::vector<double> pos;  //..Particle position
        std::vector<double> vel;  //..Particle velocity
        std::vector<double> en;   //..Particle energy
        std::vector<int> cell;  //..Nearest cell
        std::vector<int> pt;  //..Nearest point
        std::vector<double> cseedcount;  //..Counter for cell seed fraction
        std::vector<double> cllsncount; //..Counter for collision probability for each cell

        double mass;  //..Particle mass
        double wmass;  //..Particle weighted mass
        double charge;  //..Particle charge
        double wcharge;  //..Particle weighted charge
        double  pwght;  //..Particle weight
        double  nu_max;  //..Maximum Collision Frequency
        double  crsssctn_max;  //..Maximum Cross Section
        double  ncolcount;  //..Neutral collision count parameter
        double  cllsnmass;  //..Mass of neutral species for collision
        double  cllsnB;  //..Experimentally determine collision parameter B
        double fseedcount;  //..Counter for flux seed fraction
        bool  mag;  //..Particle magnetization Q1D
        int  spflag;  //..Species flag
        int  nsp;  //..Number of species
        int  gp;  //..Number of ghost particles
        int  pid;  //..Species id
        int  seedproc;  //..MPI  Processor for seeding
        int  elcount;   //..Elastic Collision Count
        int  exccount;  //..Excitation Collision Count
        int  ioncount;  //..Ionization Collision Count

        int  cycIT;  //..Subcycle Iteration

        int  pnct;  //..Number of collision type for particle
	std::vector<std::string>  cllsntype; //..Collision type
	std::vector<std::string>  crsstype; //..Cross-section type
	std::vector<double>  crsssctn; //..Cross-section size
	std::vector<double>  en_crsstable; //..Cross-section lookup table
	std::vector<double>  el_crsstable; //..Cross-section lookup table
	std::vector<double>  inel_crsstable; //..Cross-section lookup table
	std::vector<double>  ion_crsstable; //..Cross-section lookup table
	std::vector<double>  cllsnenergy; //..Collision energies

};

struct fields
{
        std::vector<double> phi;  //..Electric potential
	std::vector<double> E;  //..Electric field
	std::vector<double> B;  //..Magnetic field
	std::vector<double> gradBcoef;  //..Gradient of magnetic field coefficient
        int EM_flag;  //..
};


struct contnm
{
        std::string name;  //..Name of species
        double charge;  //..Charge of species
        double wcharge;  //..Weighted charge of species
        double mass;  //..Mass of species
        double wmass;  //..Weighted mass of species
        double  pwght;  //..Particle weight
        double  R;  //..Gas Constant (only used for Neutrals)
        double  NL;  //..Left Boundary Gas Density (for Neutrals only)
        double  TL;  //..Left Boundary Gas Temperature (for Neutrals only)
        int nsp;   //..Number of species
        int  spflag;  //..Species flag

	std::vector<double> N;   //..Number density
	std::vector<double> ioncount;   //..Number of Ionization Events
        std::vector<double> mom;  //..Macroscopic Momentum
	std::vector<double> U;   //..Macroscopic Velocity
	std::vector<double> T;   //..Temperature
	std::vector<double> Tvec;   //..Directional Temperature
	std::vector<double> Qvec;   //..Directional heat flow
	std::vector<double> p;  //..Pressure
	std::vector<double> En;  //..Energy
	std::vector<double> rho_back;  //..Background charge
};

struct MPIvars
{
	std::vector<double> Ntotal;   //..Number density (for MPI)
	std::vector<double> N_all;   //..Number density (for MPI)
	std::vector<double> Utotal;   //..Macroscopic Velocity (for MPI)
	std::vector<double> U_all;   //..Macroscopic Velocity (for MPI)
	std::vector<double> Tvectotal;   //..Directional Temperature (for MPI)
	std::vector<double> Qvectotal;   //..Heat flow Temperature (for MPI)
	std::vector<double> Tvec_all;   //..Directional Temperature (for MPI)
	std::vector<double> Qvec_all;   //..Heat flow Temperature (for MPI)
	std::vector<double> Entotal;  //..Energy(for MPI)
	std::vector<double> En_all;  //..Energy (for MPI)
};

struct flow
{
	std::vector<double> rho;  //..Charge density
	std::vector<double> N;  //..Number density
	std::vector<double> iN; //..Ion number density
	std::vector<double> eN;  //..Electron number density
        std::vector<double> mom;  //..Momentum
	std::vector<double> U;  //..Velocity
	std::vector<double> iU;  //..Ion velocity
	std::vector<double> eU;  //..Electron velocity
	std::vector<double> T;  //..Temperature
	std::vector<double> iT;  //..Ion Temperature
	std::vector<double> eT;  //..Electron Temperature
	std::vector<double> p;   //..Pressure
	std::vector<double> En;  //..Energy

        double gam;  //..Gamma
        int spflag;  //..species flag
};

struct boundvars
{
        int nbound;  //..Number of boundaries
        int nsp;  //..Number of species
        int nDN;  //..Number of Dirichlet and Neumann boundary conditions
        int wallflag;  //..Type of wall
        int wallfloat;  //..Float wall
        int wallcur;  //..Current source wall
        int wallcap;  //..Capacitor wall
        int wallvolt;  //..Voltage driven wall
        double Jfreq;  //..Current frequency wall
        double Vfreq;  //..Voltage frequency wall
        double cap;  //..capacitance
        std::string walltype;  //..Type of wall string
        std::string srctype;  //..Type of source

        std::vector<double> boundrange;  //..Range for boundary
        std::string boundname, clss;  //..Boundary name and class
        std::vector<std::string> boundtype;  //..Boundary type
        std::vector<std::string> bddens;  //..Density
        std::vector<std::string> bdddist; //..Dens. Distribution
        std::vector<std::string> bdtemp;  //..Temperature
        std::vector<std::string> bdthdist;  //..Velocity distribution
        std::vector<std::string> bdvel;  //..Velocity
        std::vector<std::string> bdE;  //..Electric field
        std::vector<std::string> bdB;  //..Magnetic field
        std::vector<std::string> bdrho;  //..Density
        std::vector<std::string> bdU;  //..Macroscopic Velocity
        std::vector<std::string> bdp; //..Pressure
        std::string bdphi;  //..Potential

        int cboundnum;  //..Cell center boundary number
        int pboundnum;  //..Corner point boundary number
        int fdir;  //..Flux Direction
        std::vector<int> cboundcells; //..Center boundary cells
        std::vector<int> pboundcells; //..Corner point boundary cells
        std::vector<int> partcount;  //..Particle count
        std::vector<double> bdneum; //..Neumann boundary conditions

        std::vector<double> bddens_val;  //..Density Value
        std::vector<double> bdU_val;  //..Macroscopic Velocity
        double bdphi_val;  //..Phi Value
        double sigma;  //..surface charge density
        double Ebd;  //..surface charge density

};

struct spclvars
{
	int nspcl;  //..Number of special regions
	double spclJ;  //..Total current
	double spclE;  //..Total Electric Field
	double spclomega;  //..Frequency
	std::string spcltype;   //..Type of Special Region
	double spclFlux_val;   //..Flux for special region
	double spclT_val;   //..Temperature for special region

	std::vector<std::string> spcldens;       //..Density of special region
	std::vector<std::string> spclddist;     //..Density distribution
	std::vector<std::string> spclT;        //..Temperature
	std::vector<std::string> spclthdist;  //..Velocity istribution
	std::vector<std::string> spclvel;    //..Velocity
	std::vector<std::string> spclU;      //..Velocity
	std::vector<std::string> spclFlux;  //..Number flux (particles/sec)

	std::vector<double> spclrange;    //..Range for special region
	std::vector<int> spclcells;   //..Cells for special region
	std::vector<int> spclpoints;   //..Points for special region

        double spclEperp;  //..Special region Efield

};


#endif

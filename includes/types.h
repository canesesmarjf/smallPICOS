#ifndef H_TYPES
#define H_TYPES

#include <iostream>
#include <vector>
#define ARMA_ALLOW_FAKE_GCC
#include <armadillo>
#include <string>
#include <cmath>

using namespace std;

// Physical constants
// =============================================================================
#define PRO_ZERO 1.0E-15	// Definition of zero in PROMETHEUS
#define F_E 1.602176E-19	// Electron charge in C (absolute value)
#define F_ME 9.109382E-31	// Electron mass in kg
#define F_MP 1.672621E-27	// Proton mass in kg
#define F_U 1.660538E-27	// Atomic mass unit in kg
#define F_KB 1.380650E-23	// Boltzmann constant in Joules/Kelvin
#define F_EPSILON 8.854E-12 // Vacuum permittivity in C^2/(N*m^2)
#define F_C 299792458.0 	// Light speed in m/s
#define F_MU (4*M_PI)*1E-7 	// Vacuum permeability in N/A^2

// Class to hold ion parameters:
// =============================================================================
class ions_params_TYP
{
public:
  int SPECIES;
  int N_CP; // Number of computational particles over entire simulation
  int pct_N_CP_Output; // percentage of N_CP recorded in output file
  int Z;
  double M;
};

// Class to hold the initial conditions PARAMETERS for each ion species:
// =============================================================================
class ions_IC_params_TYP
{
public:
    string fileName;
    double Tpar_offset;
    double Tpar_scale;
    double Tper_offset;
    double Tper_scale;
    double upar_offset;
    double upar_scale;
    double densityFraction;
    double mean_ai; // Mean particle weight

    string CP_fileName; // Filename that indicates where the computational particle density profile is to be found.
};

// Class to hold the boundary conditions PARAMETERS for each ion species:
// =============================================================================
struct ions_BC_params_TYP
{
    int type;
    double T;
    double E;
    double eta;
    double mean_x;
    double sigma_x;

    double G;
    string G_fileName; // file that stores time dependent source rate
    //int    G_NS;       // Size of file
};

// Class to hold the initial conditions PARAMETERS for electrons:
// =============================================================================
class electrons_IC_params_TYP
{
public:
    string fileName;
    double n_offset;
    double n_scale;
    double T_offset;
    double T_scale;
};

// Class to hold the initial conditions PARAMETERS for fields:
// =============================================================================
class fields_IC_params_TYP
{
public:
    string fileName;
    double Ex_offset;
    double Ex_scale;
    double Bx_offset;
    double Bx_scale;
};

// Class to store mesh PARAMETERS:
// =============================================================================
class mesh_params_TYP
{
public:
  double dx_norm;
  double x0;
  double r0_min;
  double r0_max;
  double Lx_min;
  double Lx_max;

  double ionSkinDepth;
  int Nx;
  double dx;
  double A0;
  double B0;

  // Methods:
  void getA0();
};

// Class to hold initial conditon PROFILES for each ion species:
// =============================================================================
class ions_IC_TYP
{
public:
    // Store profiles as given by H5 file:
    arma::vec x;
    arma::vec n;
    arma::vec Tpar;
    arma::vec Tper;
    arma::vec upar;
    arma::vec ncp_pdf; // PDF that defines initial computational particle density profile

    // Store profiles interpolated at the cell-centers of the mesh with ghost cells included:
    arma::vec x_mg;
    arma::vec n_mg;
    arma::vec Tpar_mg;
    arma::vec Tper_mg;
    arma::vec upar_mg;
    arma::vec ncp_pdf_mg;

    // IC profiles calculated in the code (with ghost cells included):
    arma::vec a_mg;
    arma::vec ncp_mg;
};

// Class to hold the initial conditions PROFILES for the electrons:
// =============================================================================
struct electrons_IC_TYP
{
    // Store profiles as given by H5 file:
    arma::vec x;
    arma::vec n;
    arma::vec T;

    // Store profiles interpolated at the cell-centers of the mesh:
    arma::vec x_mg;
    arma::vec n_mg;
    arma::vec T_mg;
};

// Class to hold the initial conditions PROFILES for the fields:
// =============================================================================
struct fields_IC_TYP
{
    // Store profiles as given by H5 file:
    arma::vec x;
    arma::vec Bx;
    arma::vec Ex;
    arma::vec dBx;
    arma::vec ddBx;

    // Store profiles interpolated at the cell-centers of the mesh:
    arma::vec x_mg;
    arma::vec Bx_mg;
    arma::vec Ex_mg;
    arma::vec dBx_mg;
    arma::vec ddBx_mg;
};

// Class to hold all initial condition PROFILES:
// =============================================================================
class IC_TYP
{
public:
    vector<ions_IC_TYP> ions;
    electrons_IC_TYP electrons;
    fields_IC_TYP fields;
};

// Class to represent each simulated ion species:
// =============================================================================
class ions_TYP
{

  // Notes:
  // Consider double get_K();
  // Consider double get_NR(); Keep in mind that this calculation will be local to the MPI
  // To get the real N_R we would need to reduce over all MPI process and then broadcast

public:
    int SPECIES;   // 0: tracer 1: GC particles
    uint N_CP;     // Number of computational particles in ENTIRE simulation
    uint N_CP_MPI; // Number of computational particles per MPI process
    double K;      // Distribution function normalization constant K = N_R/N_SP

    double N_R;    // Number of real particles in ENTIRE simulation (Dynamic)
    double N_SP;   // Number of super-particles in ENTIRE simulation (Dynamic)

    double Z;      // Atomic number
    double M;      // Mass
    double Q;      // Charge

    // Variables to control computational particle output:
    double pct_N_CP_Output; // Fraction of computational particles from ENTIRE simulation to record in output file
    int N_CP_Output;        // Number of computational particles from ENTIRE simulation to record in output file
    int N_CP_Output_MPI;    // Number of computational particles per MPI to record in output file

    // Particle-defined quantities:
    // ============================
    // Each MPI process has its own ensemble of particles and associated attributes as recorded below:
    // Attributes:
    arma::vec x_p;
    arma::mat v_p;
    arma::vec a_p;

    // Nearest grid point
    arma::ivec mn;

    // Fields:
    arma::vec Ex_p;
    arma::vec Bx_p;
    arma::vec dBx_p;
    arma::vec ddBx_p;

    // Moments: (Needed for collision operator)
    arma::vec n_p;
    arma::vec nv_p;
    arma::vec Tpar_p;
    arma::vec Tper_p;
    arma::vec Te_p;

    // Assignment function:
    arma::vec wxl;
    arma::vec wxc;
    arma::vec wxr;

    // Mesh-defined quantities:
    // ============================
    // These quantities record the values reduced over all MPI processes.
    // Each MPI process will have a partial mesh-defined quantity which then needs to be reduced over all MPI processes to get the final value for the mesh-defined quantity.
    arma::vec n_m;
    arma::vec nv_m;
    arma::vec P11_m;
    arma::vec P22_m;
    arma::vec Tpar_m;
    arma::vec Tper_m;

    arma::vec ncp_m;
    arma::vec ncp_m_MPI; // Computational particle density in each MPI process
};

// Class to represent electrons species:
// =============================================================================
class electrons_TYP
{
public:

    // Mesh-defined temperature:
    arma::vec Te_m;

    // Constructor:
    electrons_TYP(){};

    // Destructor:
    ~electrons_TYP(){};
};

// Class to represent electromagnetic fields in the simulation:
// =============================================================================
class fields_TYP
{
public:

    arma::vec Ex_m;
    arma::vec Bx_m;
    arma::vec dBx_m;
    arma::vec ddBx_m;

    arma::vec Am;

    fields_TYP(){};

    //fields_TYP(unsigned int N) : Ex_m(N), Bx_m(N), dBx_m(N), ddBx_m(N) {};

    ~fields_TYP(){};

    //void zeros(unsigned int N);
    //void fill(double A);
    void getAm(double A0, double B0);
};

//  Define structure to store characteristic values for the normalization:
// =============================================================================
struct CV_TYP
{
	double ne;
	double Te;
	double B;
	double Tpar;
	double Tper;

	CV_TYP()
	{
		ne   = 0;
		Te   = 0;
		B    = 0;
		Tpar = 0;
		Tper = 0;
	}
};

//  Define structure to store switches that control physics modules:
// =============================================================================
struct SW_TYP
{
	int EfieldSolve;
	int BfieldSolve;
	int Collisions;
	int RFheating;
	int linearSolve;
	int advancePos;

	SW_TYP()
	{
		EfieldSolve   = 0;
		BfieldSolve   = 0;
		Collisions    = 0;
		RFheating     = 0;
		linearSolve   = 0;
		advancePos    = 0;
	}

};

// Class to store position of mesh cell-centers:
class mesh_TYP
{
public:
  int Nx;
  double dx;
  arma::vec xm;
  arma::vec xmg;
  arma::vec Am;
  arma::vec Bxm;

  // Overloaded constructor:
  // mesh_TYP(params_TYP * params);

};

// Struct to store simulation PARAMETERS:
// =============================================================================
struct params_TYP
{
    // List of variables in the outputs
    vector<string> outputs_variables;

    // Select method for particle RK4 integrator
    int advanceParticleMethod;

    //Control parameters for the simulation:
    // Path to save the outputs
    string PATH;

    int argc;
    char **argv;

    // Flag for using a quiet start
    // bool quietStart;

    double smoothingParameter;
    double simulationTime;
    double currentTime = 0;
    int timeIterations;

    // Consider eliminating one of the following:
    // Or call the 2nd one dt_norm
    double DT;//Time step
    double DTc;//Cyclotron period fraction.

    // Needs redefinition as we dont use the cyclotron period as reference anymore:
    double outputCadence;//Save variables each "outputCadence" times the background ion cycloperiod.

    // Mesh parameters:
    mesh_params_TYP mesh_params;

    // Ions properties:
    int numberOfParticleSpecies; // This species are evolved self-consistently with the fields
    int numberOfTracerSpecies; // This species are not self-consistently evolved with the fields

    // Ion parameters:
    vector<ions_params_TYP> ions_params;

    // Initial conditions parameters:
    electrons_IC_params_TYP electrons_IC;
    fields_IC_params_TYP fields_IC;
    vector<ions_IC_params_TYP> ions_IC;

    // Boundary condition parameters:
    vector<ions_BC_params_TYP> ions_BC;

    // Simulation Characterstic values:
    CV_TYP CV;

    // Simulation switches:
    SW_TYP SW;

    // RF operator conditions:
    //RF_TYP RF;

    int filtersPerIterationFields;
    int filtersPerIterationIons;

    // How are these used?
    //double ionLarmorRadius;
    //double ionSkinDepth;
    //double ionGyroPeriod;

    //double DrL;
    //double dp;

    // How are the following used?
    //int checkStability;
    //int rateOfChecking;//Each 'rateOfChecking' iterations we use the CLF criteria for the particles to stabilize the simulation

    // MPI parameters
    // mpi_params_TYP mpi;

    // Error codes
    map<int,string> errorCodes;

    //  Methods:
    void getCharacteristicIonSkinDepth();
    void get_Nx_dx(double ionSkinDepth);

    // Constructor
    params_TYP(){};

};

#endif

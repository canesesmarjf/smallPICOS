#ifndef H_INITIALIZE
#define H_INITIALIZE

// Intrinsic header files:
// =======================
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <typeinfo>
#include <sstream>

// Armadillo header:
// =================
#define ARMA_ALLOW_FAKE_GCC
#define ARMA_USE_HDF5
#include <armadillo>

// User-defined headers:
// =====================
#include "types.h"
//#include "quietStart.h"
//#include "PIC.h"
//#include "mpi_main.h"

using namespace std;
using namespace arma;

class init_TYP
{
    vector<string> split(const string& str, const string& delim);

    map<string,string> readTextFile(string * inputFile);

    public:

    init_TYP(params_TYP * params, int argc, char* argv[]);

    void read_inputFile(params_TYP * params);

    void read_ionsPropertiesFile(params_TYP * params);

    void read_IC_profiles(params_TYP * params, mesh_TYP * mesh, IC_TYP * IC);

    void interpolate_IC_profiles(params_TYP * params, mesh_TYP * mesh, IC_TYP * IC);

    void create_mesh(params_TYP * params, mesh_TYP * mesh);

    void initialize_fields(params_TYP * params, IC_TYP * IC, fields_TYP * fields);

    void initialize_electrons(params_TYP * params, IC_TYP * IC, electrons_TYP * electrons);

    void calculate_IC_particleWeight(params_TYP * params, IC_TYP * IC, vector<ions_TYP> * IONS);

    void initialize_ions(params_TYP * params, IC_TYP * IC, mesh_TYP * mesh, vector<ions_TYP> * IONS);

    void allocate_meshDefinedIonArrays(params_TYP * params, ions_TYP * IONS);

    void allocate_particleDefinedIonArrays(params_TYP * params, ions_TYP * IONS);
};

#endif

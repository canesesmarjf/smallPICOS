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

    void readInputFile(params_TYP * params);

    void readIonPropertiesFile(params_TYP * params);

    void readInitialConditionProfiles(params_TYP * params, IC_TYP * IC);

    void calculateGlobalQuantities(params_TYP * params, IC_TYP * IC, vector<ions_TYP> * IONS);

};

#endif

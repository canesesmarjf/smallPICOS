#include<iostream>
#include<vector>

#define ARMA_ALLOW_FAKE_GCC

#include<armadillo>
#include<cmath>
#include <ctime>
#include <utility>

// User defined header files:
#include "types.h"
#include "initialize.h"

// Namespaces used in code:
using namespace std;
using namespace arma;

// Constants:
# define PI arma::datum::pi

int main(int argc, char* argv[])
{
    // Input parameters for simulation:
    params_TYP params;

    // Ion species vector:
    vector<ions_TYP> IONS;

    // Electron species object:
    electrons_TYP electrons;

    // Electromagnetic fields:
    fields_TYP fields;

    // Initialization object:
    init_TYP init(&params, argc, argv);

    // Object to hold IC condition profiles:
    IC_TYP IC;

    // Read input files:
    // =========================================================================
    // 1- Read "input_file.input" into "params"
    // 2- Read "ions_properties.ion" into "params":
    // 3- Read IC profiles from external files into "IC":
    init.readInputFile(&params);
    init.readIonPropertiesFile(&params);
    init.readInitialConditionProfiles(&params, &IC);

}

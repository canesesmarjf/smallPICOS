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
#include "output_H5.h"
#include "PIC.h"

// Namespaces used in code:
using namespace std;
using namespace arma;

// Constants:
# define PI arma::datum::pi

// Notes:
// - DTc could be changed to dt_norm
// - DT could be changed to dt
// - The use of FS could be eliminated and just have either electrons_TYP and ions_TYP have a method to make the spatial and temporal checks.
// Reconsider the use of CV.Te, CV.Tpar and CV.Tper in the simulation. Are they really needed?

int main(int argc, char* argv[])
{
    // Initialize MPI process:
    // =========================================================================
    //MPI_Init(&argc, &argv);

    // MPI object to hold topology information:
    // MPI_MAIN_TYP mpi_main;

    // Input parameters for simulation:
    params_TYP params;

    // Ion species vector:
    vector<ions_TYP> IONS;

    // Electron species object:
    electrons_TYP electrons;

    // Electromagnetic fields:
    fields_TYP fields;

    // Characteristic scales:
    // CS_TYP CS;

    // Mesh object:
    mesh_TYP mesh;

    // Object to hold IC condition profiles:
    IC_TYP IC;

    // UNITS object:
    // units_TYP units;

    // Initialization object:
    init_TYP init(&params, argc, argv);

    // =========================================================================
    // - Read "input_file.input" into "params"
    init.read_inputFile(&params);

    // - Read "ions_properties.ion" into "params":
    init.read_ionsPropertiesFile(&params);

    // - Create mesh using input parameters:
    init.create_mesh(&params,&mesh);

    // Create MPI topology:
    // mpi_main.createMPITopology(&params);

    // - Read IC profiles from external files into "IC":
    init.read_IC_profiles(&params,&mesh,&IC);

    // - Interpolate IC profiles to mesh grid:
    init.interpolate_IC_profiles(&params,&mesh,&IC);

    // Calculate particle weight initial condition profile:
    init.calculate_IC_particleWeight(&params,&IC,&IONS);

    // - Initialize fields using IC field profiles:
    init.initialize_fields(&params,&IC,&fields);

    // - Initialize electrons using IC profiles:
    init.initialize_electrons(&params,&IC,&electrons);

    // - Initialize ions using IC profiles:
    init.initialize_ions(&params,&IC,&mesh,&IONS);

    // - Define characteristic scales:
    // units.defineCharacteristicScalesAndBcast(&params, &IONS, &CS);

    // - HDF object constructor and create "main.h5"
    // HDF_TYP HDF(&params, &FS, &IONS);
    HDF_TYP HDF;

    // - Define time step based on ion CFL condition:
    // units.defineTimeStep(&params, &IONS);

    // - Normalize "params", "IONS", "electrons", "fields" using "CS" (NEED to add mesh)
    // units.normalizeVariables(&params, &IONS, &electrons, &fields, &CS);

    // - Create variables for tracking similation time:
    // double t1 = 0.0;
    // double t2 = 0.0;
    params.currentTime = 0.0;
    // int outputIterator = 0;
    // int numberOfIterationsForEstimator = 1000;

    // - Create EM solver solver:
    // fields_solver_TYP fields_solver(&params, &CS);

    // #
    // - Create PIC object:
    PIC_TYP PIC(&params, &mesh, &fields, &IONS, &electrons);

    // - Create RF operator object:
    // RF_Operator_TYP RF_operator(&params,&CS,&fields,&IONS);

    // #
    // - Save 1st output:
    // HDF.saveOutputs(&params, &IONS, &electrons, &fields, &CS, 0, 0);
    string fileName = "file_1.h5";
    HDF.saveData(fileName,&params,&fields,&IONS);

    // - Start timing simulations:
    // t1 = MPI_Wtime();

    // - Start time iterations:
    // for(int tt=0; tt<params.timeIterations; tt++)
    {

      // - Advance particles and re-inject:
      // if (params.SW.advancePos == 1)
      {
          // - Advance particle position and velocity to level X^(N+1):
          // PIC.advanceParticles(&params,&mesh,&fields, &IONS);

          // - Re-inject particles that leave computational domain:
          // particleBC.applyParticleReinjection(&params,&CS,&fields,&IONS);

          // #
          // - Assign cell:
          // PIC.assignCell_AllSpecies(&params,&mesh,&IONS);

          // #
          // - Interpolate all fields:
          // PIC.interpolateFields_AllSpecies(&params,&IONS,&fields);

          // #
          // - Interpolate electron temperature:
          // PIC.interpolateElectrons_AllSpecies(&params,&IONS,&electrons);
      }

      // #
      // - Calculate ion moments:
      // PIC.extrapolateMoments_AllSpecies(&params,&CS,&fields,&IONS);

      // - Apply collision operator:

      // - Apply RF operator:

      // - Calculate new electric field:

      // - Advance time:
      // params.currentTime += params.dt*CS.time;

      // #
      // - Save data:

      // - Estimate simulation time:
      /*
      if(tt == numberOfIterationsForEstimator)
      {
          t2 = MPI_Wtime();

          double estimatedSimulationTime = ( (double)params.timeIterations*(t2 - t1)/(double)numberOfIterationsForEstimator )/60.0;

          //if(params.mpi.MPI_DOMAIN_NUMBER == 0)
          {
              cout << "ESTIMATED TIME OF COMPLETION: " << estimatedSimulationTime <<" MINUTES" << endl;
          }
      }*/

    } // End time iterations

  //- Finalize MPI communication:

  return(0);
}

// NOTES:
// - Save data to H5 files for post-process:
// The main variables to save are:
// x_p, v_p, a_p, n_m
// hdf5.save(&params,&mesh,&electrons,&IONS);

// - Resample distribution:

// - Save data to H5 files for post-process:
// hdf5.save(&params,&mesh,&electrons,&IONS);

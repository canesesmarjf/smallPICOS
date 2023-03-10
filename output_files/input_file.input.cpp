// Input data:
// =============================================================================
mpisForFields               2
//quietStart                  1
numberOfParticleSpecies     2
numberOfTracerSpecies       0
advanceParticleMethod       1

// Characteristic values:
// =============================================================================
CV_ne                       5.0E19
CV_Te                       15
CV_B                        1.5
CV_Tpar                     15
CV_Tper                     15

// Simulation time:
// =============================================================================
DTc                         0.5
simulationTime              25000

// Switches:
// =============================================================================
SW_EfieldSolve              1
SW_BfieldSolve              0
SW_Collisions               1
SW_RFheating                1
SW_advancePos               1
//SW_linearSolve             0

// Geometry:
// =============================================================================
//GEO_dp                    0.4243
GEO_dx_norm                 0.4243
GEO_x0                      0.0
GE0_r0_min                  0.0
GEO_r0_max                  0.05
GEO_Lx_min                  -3
GEO_Lx_max                  +3

// Fields initial conditions:
// =============================================================================
IC_fields_fileName          fields_IC.h5
IC_Ex_offset                0.0
IC_Ex_scale                 1.0
IC_Bx_offset                0.0
IC_Bx_scale                 1.0
//IC_uniformBfield            0
//IC_BX                       0.2
//IC_BY                       0.0
//IC_BZ                       0.0
//IC_BX_NX                    200
//IC_BX_fileName              MPEX_B_norm_PICOS_scenario_14.txt

// Electrons initial conditions:
// =============================================================================
IC_electrons_fileName       electrons_IC.h5
IC_ne_offset                0.0
IC_ne_scale                 1.0
IC_Te_offset                0.0
IC_Te_scale                 1.0

// RF operator:
// =============================================================================
RF_Prf                      50E3
RF_n_harmonic               1
RF_freq                     8.385E6
RF_x1                       4.0
RF_x2                       6.5
RF_t_ON                     12000
RF_t_OFF                    20000
RF_kpar                     20
RF_kper                     100
RF_handedness               -1
RF_Prf_fileName             Prf_profile.txt
RF_Prf_NS                   200

// Output variables:
// =============================================================================
outputCadence               500
outputs_variables           {X_p,V_p,a_p,BX_p,BX_m,n_m,Tpar_m,Tper_m,u_m,EX_m}

// Data smoothing:
// =============================================================================
smoothingParameter          1.0E-4
filtersPerIterationFields   2
filtersPerIterationIons     2

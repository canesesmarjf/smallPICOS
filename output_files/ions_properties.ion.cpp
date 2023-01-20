
// Specify file with ion profiles:
// =============================================================================
IC_ionProfiles_fileName       Ions_IC.h5

// =============================================================================
// *****************************     SPECIES 1     *****************************
// =============================================================================

// Species 1, general:
// =============================================================================
SPECIES1                      1
NPC1                          2500
pctSupPartOutput1             100
Z1                            1
M1                            2.0

// Species 1, initial conditions:
// =============================================================================
IC_Tpar_offset_1              0.0
IC_Tpar_scale_1               1.0

IC_Tper_offset_1              0.0
IC_Tper_scale_1               1.0

IC_upar_offset_1              0.0
IC_upar_scale_1               1.0

IC_densityFraction_1          0.8

// Species 1, boundary conditions:
// =============================================================================
// 1: Warm plasma source, 2: NBI, 3: Periodic, 4: simple-reinjection
BC_type_1                     1

BC_T_1                        15
BC_E_1                        0
BC_eta_1                      45
BC_mean_x_1                   0
BC_sigma_x_1                  0.4

BC_G_1                        1E22
BC_G_fileName_1               G_profile.txt
BC_G_NS_1                     200

// =============================================================================
// *****************************     SPECIES 2     *****************************
// =============================================================================

// Species 2, general:
// =============================================================================
SPECIES2                      1
NPC2                          2500
pctSupPartOutput1             100
Z2                            1
M2                            3.0

// Species 2, initial conditions:
// =============================================================================
IC_Tpar_offset_2               0.0
IC_Tpar_scale_2                1.0

IC_Tper_offset_2               0.0
IC_Tper_scale_2                1.0

IC_upar_offset_2               0.0
IC_upar_scale_2                1.0

IC_densityFraction_2           0.2

// Species 2, boundary conditions:
// =============================================================================
// 1: Warm plasma source, 2: NBI, 3: Periodic, 4: simple-reinjection
BC_type_2                      1

BC_T_2                         15
BC_E_2                         0
BC_eta_2                       45
BC_mean_x_2                    0
BC_sigma_x_2                   0.4

BC_G_2                         1E22
BC_G_fileName_2                G_profile.txt
BC_G_NS_2                      200

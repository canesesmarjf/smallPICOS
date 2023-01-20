
#include <initialize.h>
#include <exception>

// #########################################################################################################

// Function to split strings:
vector<string> init_TYP::split(const string& str, const string& delim)
{
    vector<string> tokens;
    size_t prev = 0, pos = 0;
    do
    {
        pos = str.find(delim, prev);

        if (pos == string::npos)
        {
            pos = str.length();
        }

        string token = str.substr(prev, pos-prev);

        if (!token.empty())
        {
            tokens.push_back(token);
        }

        prev = pos + delim.length();
    }
    while (pos < str.length() && prev < str.length());

    return tokens;
}

// #########################################################################################################

// Function to read and load data from inputfile.input:
map<string,string> init_TYP::readTextFile(string * inputFile)
{
    // Create stream object:
    // =====================
    fstream reader;

    // Create map object:
    // ==================
    std::map<string,string> readMap;

    // Open input file using reader object:
    // ====================================
    reader.open(inputFile->data(),ifstream::in);

    // Handle error:
    // =============
    if (!reader)
    {
        //MPI_Barrier(MPI_COMM_WORLD);

        cerr << "PICOS++ ERROR: The input file couldn't be opened." << endl;
        std::terminate();
        //MPI_Abort(MPI_COMM_WORLD, -101);
    }

    // Parse through file:
    // ===================
    string lineContent;
    vector<string> keyValuePair;
    while ( reader.good() )
    {
        // Read entire line:
        getline(reader,lineContent);

        // Search for comment symbol:
        size_t commentCharPos = lineContent.find("//",0);

        // Check for comment symbol:
        if (commentCharPos == 0 || lineContent.empty())
        {
            // Skip line
        }
        else
        {
            // Get value pair:
            keyValuePair = split(lineContent," ");

            // Update map:
            readMap[ keyValuePair[0] ] = keyValuePair[1];
        }
    }

    // Close stream object:
    // ===================
    reader.close();

    // Return map:
    // ==========
    return readMap;
}

// #########################################################################################################
// Constructor:
init_TYP::init_TYP(params_TYP * params, int argc, char* argv[])
{
    // Get RANK and SIZE of nodes within COMM_WORLD:
    // =============================================
    //MPI_Comm_size(MPI_COMM_WORLD, &params->mpi.NUMBER_MPI_DOMAINS);
    //MPI_Comm_rank(MPI_COMM_WORLD, &params->mpi.MPI_DOMAIN_NUMBER);

    // Error codes:
    // ============
    //params->errorCodes[-100] = "Odd number of MPI processes";
    //params->errorCodes[-102] = "MPI's Cartesian topology could not be created";
    params->errorCodes[-103] = "Grid size violates assumptions of hybrid model for the plasma -- DX smaller than the electron skin depth can not be resolved";
    params->errorCodes[-106] = "Inconsistency in inital ion's velocity distribution function";

    // Program information:
    // ===========================
    //if (params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
        cout << "* PICOS++, a 1D-2V GC hybrid PIC code for Open plasma Systems           *" << endl;
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *" << endl;
        cout << endl;
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    // Arguments and paths to main function:
    // =====================================
    params->PATH = argv[2];
    params->argc = argc;
    params->argv = argv;

    // Check number of MPI domains:
    // ============================
    /*if( fmod( (double)params->mpi.NUMBER_MPI_DOMAINS, 2.0 ) > 0.0 )
    {
        MPI_Barrier(MPI_COMM_WORLD);

        if(params->mpi.MPI_DOMAIN_NUMBER == 0)
        {
            cerr << "PICOS++ ERROR: The number of MPI processes must be an even number." << endl;
        }

        MPI_Abort(MPI_COMM_WORLD,-100);
    }*/

    // Stream date when simulation is started:
    // =======================================
    //if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        time_t current_time = std::time(NULL);
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
        cout << "STARTING " << params->argv[1] << " SIMULATION ON: " << std::ctime(&current_time) << endl;
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * " << endl;
    }

}

// #########################################################################################################

// Temporary function created to check the contents of map:
void CheckMap(std::map<string,string> parametersStringMap, string key)
{
  std::map<string,string>::iterator it;
  it = parametersStringMap.find(key);
  if (it != parametersStringMap.end())
  {
      std::cout << "value for key '"<< it->first << "' is '" << it->second << "'" << std::endl;
  }
  else
  {
      std::cout << "key not found" << std::endl;
  }
}

// #########################################################################################################

// Populate params with data from input file:
void init_TYP::readInputFile(params_TYP * params)
{
    //MPI_Barrier(MPI_COMM_WORLD);

    // Print to terminal:
    // ==================
    //if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
        cout << "READING INPUT FILE ..." << endl;
    }

    // Get name of path to input file:
    // ===============================
    string name;
    if(params->argc > 3)
    {
        string argv(params->argv[3]);
        name = "input_files/input_file_" + argv + ".input";
        params->PATH += "/" + argv;
    }
    else
    {
        name = "input_files/input_file.input";
        params->PATH += "/";
    }

    // Read input file and assemble map:
    // ================================
  	std::map<string,string> parametersStringMap;
  	parametersStringMap = readTextFile(&name);

    // Create HDF5 folders if they don't exist:
    // ========================================
  	// if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
  		string mkdir_outputs_dir = "mkdir " + params->PATH;
  		const char * sys = mkdir_outputs_dir.c_str();
  		int rsys = system(sys);

  		string mkdir_outputs_dir_HDF5 = mkdir_outputs_dir + "/HDF5";
  		sys = mkdir_outputs_dir_HDF5.c_str();
  		rsys = system(sys);
  	}

    // Populate "params" with data from input file:
    // ============================================
    //params->mpi.MPIS_FIELDS  = stoi( parametersStringMap["mpisForFields"] );

    // if(stoi( parametersStringMap["quietStart"] ) == 1)
    // {
    //     params->quietStart = true;
    // }
    // else
    // {
    //     params->quietStart = false;
    // }

    params->numberOfParticleSpecies = stoi( parametersStringMap["numberOfParticleSpecies"] );
    params->numberOfTracerSpecies   = stoi( parametersStringMap["numberOfTracerSpecies"] );
    params->advanceParticleMethod   = stoi( parametersStringMap["advanceParticleMethod"] );

    // Characteristic values:
    // -------------------------------------------------------------------------
    params->CV.ne   = stod( parametersStringMap["CV_ne"] );
    params->CV.Te   = stod( parametersStringMap["CV_Te"] )*F_E/F_KB; // In Kelvin
    params->CV.B    = stod( parametersStringMap["CV_B"] );
    params->CV.Tpar = stod( parametersStringMap["CV_Tpar"] )*F_E/F_KB; // In Kelvin
    params->CV.Tper = stod( parametersStringMap["CV_Tper"] )*F_E/F_KB; // In Kelvin

    // Simulation time:
    // -------------------------------------------------------------------------
    params->DTc            = stod( parametersStringMap["DTc"] );
    params->simulationTime = std::stod( parametersStringMap["simulationTime"] );

    // Switches:
    // -------------------------------------------------------------------------
    params->SW.EfieldSolve   = stoi( parametersStringMap["SW_EfieldSolve"] );
    params->SW.BfieldSolve   = stoi( parametersStringMap["SW_BfieldSolve"] );
    params->SW.Collisions    = stoi( parametersStringMap["SW_Collisions"] );
    params->SW.RFheating     = stoi( parametersStringMap["SW_RFheating"] );
    params->SW.advancePos    = stoi( parametersStringMap["SW_advancePos"] );
    //params->SW.linearSolve   = stoi( parametersStringMap["SW_linearSolve"] );

    // Magnetic field initial conditions:
    // -------------------------------------------------------------------------
    params->fields_IC.fileName  = parametersStringMap["IC_fields_fileName"];
    params->fields_IC.Ex_offset = stod( parametersStringMap["IC_Ex_offset"] );
    params->fields_IC.Ex_scale  = stod( parametersStringMap["IC_Ex_scale"] );
    params->fields_IC.Bx_offset = stod( parametersStringMap["IC_Bx_offset"] );
    params->fields_IC.Bx_scale  = stod( parametersStringMap["IC_Bx_scale"] );

    // Electron initial conditions:
    // -------------------------------------------------------------------------
    params->electrons_IC.fileName  = parametersStringMap["IC_electrons_fileName"];
    params->electrons_IC.n_offset = stod(parametersStringMap["IC_n_offset"]);
    params->electrons_IC.n_scale  = stod(parametersStringMap["IC_n_scale"]);
    params->electrons_IC.T_offset = stod(parametersStringMap["IC_T_offset"]);
    params->electrons_IC.T_scale  = stod(parametersStringMap["IC_T_scale"]);

    // Geometry:
    // -------------------------------------------------------------------------
    params->geometry.dx_norm = stod(parametersStringMap["GEO_dx_norm"]);
    params->geometry.x0      = stod(parametersStringMap["GEO_x0"]);
    params->geometry.r0_min  = stod(parametersStringMap["GE0_r0_min"]);
    params->geometry.r0_max  = stod(parametersStringMap["GEO_r0_max"]);
    params->geometry.Lx_min  = stod(parametersStringMap["GEO_Lx_min"]);
    params->geometry.Lx_max  = stod(parametersStringMap["GEO_Lx_max"]);

    // RF parameters
    // -------------------------------------------------------------------------
    // params->RF.Prf        = stod( parametersStringMap["RF_Prf"] );
    // params->RF.n_harmonic = stoi( parametersStringMap["RF_n_harmonic"] );
    // params->RF.freq       = stod( parametersStringMap["RF_freq"]);
    // params->RF.x1         = stod( parametersStringMap["RF_x1"]  );
    // params->RF.x2         = stod( parametersStringMap["RF_x2"]  );
    // params->RF.t_ON       = stod( parametersStringMap["RF_t_ON"]  );
    // params->RF.t_OFF      = stod( parametersStringMap["RF_t_OFF"]  );
    // params->RF.kpar       = stod( parametersStringMap["RF_kpar"]);
    // params->RF.kper       = stod( parametersStringMap["RF_kper"]);
    // params->RF.handedness = stoi( parametersStringMap["RF_handedness"]);
    // params->RF.Prf_NS     = stoi( parametersStringMap["RF_Prf_NS"] );
    // params->RF.Prf_fileName = parametersStringMap["RF_Prf_fileName"];

    // Output variables:
    // -------------------------------------------------------------------------
    params->outputCadence           = stod( parametersStringMap["outputCadence"] );
    string nonparsed_variables_list = parametersStringMap["outputs_variables"].substr(1, parametersStringMap["outputs_variables"].length() - 2);
    params->outputs_variables       = split(nonparsed_variables_list,",");

    // Data smoothing:
    // -------------------------------------------------------------------------
    params->smoothingParameter        = stod( parametersStringMap["smoothingParameter"] );
    params->filtersPerIterationFields = stoi( parametersStringMap["filtersPerIterationFields"] );
    params->filtersPerIterationIons   = stoi( parametersStringMap["filtersPerIterationIons"] );

    // Here is where we need to read the ions_properties file and populate the ion initial condition parameters class:

    // MPI_Barrier(MPI_COMM_WORLD);

    // if(params->mpi.MPI_DOMAIN_NUMBER == 0)
    {
        cout << "READING INPUT FILE COMPLETED" << endl;
        cout << "* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *\n";
    }
}

// #########################################################################################################

// Read and populate ion parameters input file:
void init_TYP::readIonPropertiesFile(params_TYP * params)
{
  // MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  // if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
    cout << "* * * * * * * * * * * * LOADING ION PARAMETERS * * * * * * * * * * * * * * * * * *\n";
    cout << "+ Number of ion species: " << params->numberOfParticleSpecies << endl;
    cout << "+ Number of tracer species: " << params->numberOfTracerSpecies << endl;
  }

  // Assemble path to "ion_properties.ion":
  // ======================================
  string name;
  if(params->argc > 3)
  {
    string argv(params->argv[3]);
    name = "input_files/ions_properties_" + argv + ".ion";
  }
  else
  {
    name = "input_files/ions_properties.ion";
  }

  // Read data from "ion_properties.ion" into "parametersMap":
  // =========================================================
  std::map<string,string> parametersMap;
  parametersMap = readTextFile(&name);

  // Determine the total number of ION species:
  // =========================================
  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

  // Allocate memory to vector holding ion IC and BC condition PARAMETERS:
  // ==================================================================
  params->ions_IC.resize(totalNumSpecies);
  params->ions_BC.resize(totalNumSpecies);
  //IONS->resize(totalNumSpecies);

  // Populate ion IC and BC conditions PARAMETERS in params:
  // ======================================================
  for(int ss = 0; ss<totalNumSpecies; ss++)
  {
    params->ions_IC[ss].fileName = parametersMap["IC_ionProfiles_fileName"];

    string name;
    stringstream kk;
    kk << ss + 1;

    // General:
    name = "SPECIES" + kk.str();
    params->ions_IC[ss].SPECIES = stoi(parametersMap[name]);
    //IONS->at(ss).SPECIES         = stoi(parametersMap[name]);

    name = "N_CP" + kk.str();
    params->ions_IC[ss].N_CP = stoi(parametersMap[name]);
    //IONS->at(ss).N_CP         = stoul(parametersMap[name]);

    name = "pct_N_CP_Output" + kk.str();
    params->ions_IC[ss].pct_N_CP_Output = stoi(parametersMap[name]);

    name = "Z" + kk.str();
    params->ions_IC[ss].Z = stoi(parametersMap[name]);
    //IONS->at(ss).Z         = stoi(parametersMap[name]);

    name = "M" + kk.str();
    params->ions_IC[ss].M = stoi(parametersMap[name]);
    //IONS->at(ss).M         = stoi(parametersMap[name]);

    // IC:
    name = "IC_Tpar_offset" + kk.str();
    params->ions_IC[ss].Tpar_offset = stod(parametersMap[name]);

    name = "IC_Tpar_scale" + kk.str();
    params->ions_IC[ss].Tpar_scale = stod(parametersMap[name]);

    name = "IC_Tper_offset" + kk.str();
    params->ions_IC[ss].Tper_offset = stod(parametersMap[name]);

    name = "IC_Tper_scale" + kk.str();
    params->ions_IC[ss].Tper_scale = stod(parametersMap[name]);

    name = "IC_upar_offset" + kk.str();
    params->ions_IC[ss].upar_offset = stod(parametersMap[name]);

    name = "IC_upar_scale" + kk.str();
    params->ions_IC[ss].upar_scale = stod(parametersMap[name]);

    name = "IC_densityFraction" + kk.str();
    params->ions_IC[ss].densityFraction = stod(parametersMap[name]);

    // BC:
    name = "BC_type" + kk.str();
    params->ions_BC[ss].type = stod(parametersMap[name]);

    name = "BC_T" + kk.str();
    params->ions_BC[ss].T = stod(parametersMap[name]);

    name = "BC_E" + kk.str();
    params->ions_BC[ss].E = stod(parametersMap[name]);

    name = "BC_eta" + kk.str();
    params->ions_BC[ss].eta = stod(parametersMap[name]);

    name = "BC_mean_x" + kk.str();
    params->ions_BC[ss].mean_x = stod(parametersMap[name]);

    name = "BC_sigma_x" + kk.str();
    params->ions_BC[ss].sigma_x = stod(parametersMap[name]);

    name = "BC_G" + kk.str();
    params->ions_BC[ss].G = stod(parametersMap[name]);

    name = "BC_G_fileName" + kk.str();
    params->ions_BC[ss].G_fileName = parametersMap[name];
  }

  // Print to terminal:
  // ==================
  // if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * * ION PARAMETERS LOADED * * * * * * * * * * * * * * * * * *\n";
  }

  // MPI_Barrier(MPI_COMM_WORLD);

}

// #########################################################################################################

// Read "ion properties" file and populate IC object with profiles:
void init_TYP::readInitialConditionProfiles(params_TYP * params, IC_TYP * IC)
{
  // MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  // if(params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << "* * * * * * * * * * * *  LOADING INITIAL CONDITION PROFILES  * * * * * * * * * * * * * * * * * *" << endl;
  }

  // Directory where profile data is located:
  string directory = "input_files/";

  // Electron profiles:
  // =====================================================
  // Assemble full path to HDF5 file with profiles:
  string fileName  = params->electrons_IC.fileName;
  string fullPath = directory + fileName;

  // temp containers:
  arma::vec y;
  double offset;
  double scale;

  // Load and assign data:
  y.load(arma::hdf5_name(fullPath,"x"));
  IC->electrons.x = y;

  y.load(arma::hdf5_name(fullPath,"n"));
  offset = params->electrons_IC.n_offset;
  scale  = params->electrons_IC.n_scale;
  IC->electrons.n = offset + scale*y;

  y.load(arma::hdf5_name(fullPath,"T"));
  offset = params->electrons_IC.T_offset;
  scale  = params->electrons_IC.T_scale;
  IC->electrons.T = offset + scale*y;

  // Fields profiles:
  // ==================================================
  fileName  = params->fields_IC.fileName;
  fullPath = directory + fileName;

  // Load and assign data:
  y.load(arma::hdf5_name(fullPath,"x"));
  IC->fields.x = y;

  y.load(arma::hdf5_name(fullPath,"Bx"));
  offset = params->fields_IC.Bx_offset;
  scale  = params->fields_IC.Bx_scale;
  IC->fields.Bx = offset + scale*y;

  y.load(arma::hdf5_name(fullPath,"Ex"));
  offset = params->fields_IC.Ex_offset;
  scale  = params->fields_IC.Ex_scale;
  IC->fields.Ex = offset + scale*y;

  // Compute 1st and 2nd derivatives of Bx:
  double dx = IC->fields.x(2) - IC->fields.x(1);
  int nx    = IC->fields.Bx.n_elem;
  arma::vec B(nx,1);
  arma::vec dB(nx,1);
  arma::vec ddB(nx,1);

  B = IC->fields.Bx;
  dB.subvec(1,nx-2) = (B.subvec(2,nx-1) - B.subvec(0,nx-3))/(2*dx);
  dB(0)    = dB(1);
  dB(nx-1) = dB(nx-2);
  IC->fields.dBx = dB;

  dB = IC->fields.dBx;
  ddB.subvec(1,nx-2) = (dB.subvec(2,nx-1) - dB.subvec(0,nx-3))/(2*dx);
  ddB(0)    = ddB(1);
  ddB(nx-1) = ddB(nx-2);
  IC->fields.ddBx = ddB;

  if (0)
  {
    IC->fields.x.save("x.txt",arma::raw_ascii);
    IC->fields.Bx.save("Bx.txt",arma::raw_ascii);
    IC->fields.dBx.save("dBx.txt",arma::raw_ascii);
    IC->fields.ddBx.save("ddBx.txt",arma::raw_ascii);
  }

  // Ions profiles:
  // ==================================================
  // Determine the total number of ION species:
  int totalNumSpecies(params->numberOfParticleSpecies + params->numberOfTracerSpecies);

  // Allocate memory:
  IC->ions.resize(totalNumSpecies);

  // Assemble full path to HDF5 file with profiles:
  fileName  = params->ions_IC.at(0).fileName;
  fullPath = directory + fileName;

  // Load data from H5 file into containers:
  for (int ss = 0; ss < totalNumSpecies; ss++)
  {
    string dataset;
    string group;
    stringstream kk;
    kk << ss + 1;
    group = "ions_" + kk.str();

    dataset = group + "/x";
    y.load(arma::hdf5_name(fullPath,dataset));
    IC->ions.at(ss).x = y;

    scale  = params->ions_IC.at(ss).densityFraction;
    IC->ions.at(ss).n = scale*IC->electrons.n;

    dataset = group + "/Tpar";
    y.load(arma::hdf5_name(fullPath,dataset));
    offset = params->ions_IC.at(ss).Tpar_offset;
    scale  = params->ions_IC.at(ss).Tpar_scale;
    IC->ions.at(ss).Tpar = offset + scale*y;

    dataset = group + "/Tper";
    y.load(arma::hdf5_name(fullPath,dataset));
    offset = params->ions_IC.at(ss).Tper_offset;
    scale  = params->ions_IC.at(ss).Tper_scale;
    IC->ions.at(ss).Tper = offset + scale*y;

    dataset = group + "/upar";
    y.load(arma::hdf5_name(fullPath,dataset));
    offset = params->ions_IC.at(ss).upar_offset;
    scale  = params->ions_IC.at(ss).upar_scale;
    IC->ions.at(ss).upar = offset + scale*y;
  }

  IC->ions.at(0).n.save("n0.txt",arma::raw_ascii);
  IC->ions.at(1).n.save("n1.txt",arma::raw_ascii);

}

// #########################################################################################################

void init_TYP::calculateGlobalQuantities(params_TYP * params, IC_TYP * IC, vector<ions_TYP> * IONS)
{
  //MPI_Barrier(MPI_COMM_WORLD);

  // Print to terminal:
  // ==================
  //if (params->mpi.MPI_DOMAIN_NUMBER == 0)
  {
      cout << endl << "* * * * * * * * * * * CALCULATING GLOBAL QUANTITIES * * * * * * * * * * * * * * * * *\n";
  }
}

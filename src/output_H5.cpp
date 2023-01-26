#include "output_H5.h"
using namespace arma;

// ===========================================================================================
void HDF_TYP::H5_writeToFile(string fileName,string groupName,string datasetName,arma::vec * v,int access)
{
  // Declare HDF5 objects:
    // ===========================================================================
    H5::H5File * file;
    H5::Group * group;
    H5::DataSpace * space;
    H5::DataSet * dataset;

    // Create file:
    // ===========================================================================
    if (access == 0) // New file
      file = new H5File(fileName,H5F_ACC_TRUNC);
    else if (access == 1) // Existing file
      file = new H5File(fileName,H5F_ACC_RDWR);

    // Create/open group:
    // ===========================================================================
    if (access == 0) // New file
      group = new H5::Group(file->createGroup(groupName));
    else if (access == 1) // Existing file
    {
      if(!file->exists(groupName))
        group = new H5::Group(file->createGroup(groupName));
      else
        group = new H5::Group(file->openGroup(groupName));
    }

    // Create dataspace:
    // ===========================================================================
    int rank = 1;
    hsize_t dims[1] = {v->n_elem};
    space = new H5::DataSpace(rank,dims);

    // Create dataset:
    // ===========================================================================
    dataset = new H5::DataSet(group->createDataSet(datasetName,PredType::NATIVE_DOUBLE,*space));


    // Write data:
    // ===========================================================================
    dataset->write(v->memptr(),PredType::NATIVE_DOUBLE);

    // Release memory:
    // ===========================================================================
    delete file;
    delete group;
    delete space;
    delete dataset;
}

// ===========================================================================================
void HDF_TYP::H5_writeToFile(string fileName,string groupName,string datasetName,arma::mat * m,int access)
{
  // Declare HDF5 objects:
    // ===========================================================================
    H5::H5File * file;
    H5::Group * group;
    H5::DataSpace * space;
    H5::DataSet * dataset;

    // Create file:
    // ===========================================================================
    if (access == 0) // New file
      file = new H5File(fileName,H5F_ACC_TRUNC);
    else if (access == 1) // Existing file
      file = new H5File(fileName,H5F_ACC_RDWR);

    // Create/open group:
    // ===========================================================================
    if (access == 0) // New file
      group = new H5::Group(file->createGroup(groupName));
    else if (access == 1) // Existing file
    {
      if(!file->exists(groupName))
        group = new H5::Group(file->createGroup(groupName));
      else
        group = new H5::Group(file->openGroup(groupName));
    }

    // Create dataspace:
    // ===========================================================================
    int rank = 2;
    hsize_t dims[2] = {m->n_cols,m->n_rows};
    space = new H5::DataSpace(rank,dims);

    // Create dataset:
    // ===========================================================================
    dataset = new H5::DataSet(group->createDataSet(datasetName,PredType::NATIVE_DOUBLE,*space));

    // Write data:
    // ===========================================================================
    dataset->write(m->memptr(),PredType::NATIVE_DOUBLE);

    // Release memory:
    // ===========================================================================
    delete file;
    delete group;
    delete space;
    delete dataset;
}

#include <cstdio>
#include <cstdlib>
#include <hdf5.h>

/* Check return values from HDF5 routines */
#define CHKERR(status,name) if ( status < 0 ) \
     fprintf( stderr, "Error: nonzero status (%d) in %s\n", status, name )

//----------------------------------------------------------------------------
// Opens an HDF5 file, determines the matrix dimensions, allocates
// memory for the matrix, and reads the matrix from the file.

void readMatrix(
    const char* filename, // in  - name of HDF5 file
    const char* path,     // in  - path within HDF5 file to matrix data
    double** a,           // out - pointer to pointer to matrix data
    int* rows,            // out - pointer to number of rows in matrix
    int* cols             // out - pointer to number of columns in matrix
    )
{
    hid_t file_id;         // HDF5 id for file
    hid_t dataspace_id;    // HDF5 id for dataspace in file
    hid_t dataset_id;      // HDF5 id for dataset in file
    hid_t memspace_id;     // HDF5 id for dataset in memory
    hsize_t* dims;         // matrix dimensions
    herr_t status;         // HDF5 return code
    int ndim;              // number of dimensions in HDF5 dataset

    // Open existing HDF5 file
    file_id = H5Fopen( filename, H5F_ACC_RDONLY, H5P_DEFAULT );
    if ( file_id < 0 ) exit( EXIT_FAILURE );

    // Open dataset in file
    dataset_id = H5Dopen( file_id, path, H5P_DEFAULT );
    if ( dataset_id < 0 ) exit( EXIT_FAILURE );

    // Determine dataset parameters
    dataspace_id = H5Dget_space( dataset_id );
    ndim = H5Sget_simple_extent_ndims( dataspace_id );
    dims = new hsize_t [ndim];

    // Get dimensions for dataset
    ndim = H5Sget_simple_extent_dims( dataspace_id, dims, NULL );
    if ( ndim != 2 )
    {
        fprintf( stderr, "Expected dataspace to be 2-dimensional " );
        fprintf( stderr, "but it appears to be %d-dimensional\n", ndim );
        exit( EXIT_FAILURE );
    }

    // Store matrix dimensions in parameter locations
    *rows = dims[0];
    *cols = dims[1];

    // Create memory dataspace
    memspace_id = H5Screate_simple( ndim, dims, NULL );
    if ( memspace_id < 0 ) exit( EXIT_FAILURE );

    // Allocate memory for matrix and read data from file
    *a = new double [dims[0] * dims[1]];
    status = H5Dread( dataset_id, H5T_NATIVE_DOUBLE, memspace_id,
                      dataspace_id, H5P_DEFAULT, *a );
    CHKERR( status, "H5Dread()" );

    // Close all remaining HDF5 objects
    CHKERR( H5Sclose( memspace_id ), "H5Sclose()" );
    CHKERR( H5Dclose( dataset_id ), "H5Dclose()" );
    CHKERR( H5Sclose( dataspace_id ), "H5Sclose()" );
    CHKERR( H5Fclose( file_id ), "H5Fclose()" );

    // Clean up
    delete [] dims;
}

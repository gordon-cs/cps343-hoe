#ifndef __READMATRIXSEQ_H__
#define __READMATRIXSEQ_H__

void readMatrix(
    const char* filename, // in  - name of HDF5 file
    const char* path,     // in  - path within HDF5 file to matrix data
    double** A,           // out - pointer to pointer to matrix data
    int* rows,            // out - pointer to number of rows in matrix
    int* cols             // out - pointer to number of columns in matrix
    );

void destroyMatrix(
    double* A             // in  - pointer to matrix data
    );

#endif // __READMATRIXSEQ_H__

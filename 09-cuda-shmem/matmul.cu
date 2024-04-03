// $Smake: nvcc -DBS=32 -O3 -o %F %f wtime.c
//
// Demonstrates use of device shared memory in matrix-matrix multiplication.
//
// Jonathan Senning <jonathan.senning@gordon.edu>
// Department of Mathematics and Computer Science
// Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
// Spring 2016, 2018.

#include <cstdio>
#include <cuda.h>
#include "wtime.h"

#define IDX(i,j,n) ((i)*(n)+j) // row major

#if !defined(BS)
const int BlockDim = 16;
#else
const int BlockDim = BS;  // needs to be 32 or less
#endif

const int MaxSizeToDisplay = 25;

typedef float FLOAT;
//typedef double FLOAT;

//----------------------------------------------------------------------------

// Matrix-matrix kernel (global memory)
__global__ void matmulGlobal(FLOAT* c, FLOAT* a, FLOAT* b, int n)
{
    // element of matrix c to compute
    const int row = blockIdx.y * blockDim.y + threadIdx.y;
    const int col = blockIdx.x * blockDim.x + threadIdx.x;
    if (col >= n || row >= n) return; // nothing to do

    FLOAT sum = (FLOAT) 0.0;
    for (int k = 0; k < n; k++)
    {
        sum += a[IDX(row,k,n)] * b[IDX(k,col,n)];
    }
    c[IDX(row,col,n)] = sum;
}

//----------------------------------------------------------------------------

// Matrix-matrix kernel (shared memory)
__global__ void matmulShared(FLOAT* c, FLOAT* a, FLOAT* b, int n)
{
    __shared__ FLOAT s_a[BlockDim][BlockDim];
    __shared__ FLOAT s_b[BlockDim][BlockDim];

    // element of matrix c to compute
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    int col = blockIdx.x * blockDim.x + threadIdx.x;

    // loop over row of blocks in matrix a and column of blocks
    // in matrix b; storing blocks in shared mem, 
    FLOAT sum = (FLOAT) 0.0;
    const int numBlocks = (n + BlockDim - 1) / BlockDim;
    for (int m = 0; m < numBlocks; m++)
    {
        // copy block from matrix to shared memory
        s_a[threadIdx.y][threadIdx.x] = 0.0;
        s_b[threadIdx.y][threadIdx.x] = 0.0;
        int r = m * blockDim.y + threadIdx.y;
        int c = m * blockDim.x + threadIdx.x;
        if (row < n && c < n) s_a[threadIdx.y][threadIdx.x] = a[IDX(row,c,n)];
        if (r < n && col < n) s_b[threadIdx.y][threadIdx.x] = b[IDX(r,col,n)];
        __syncthreads();

        // compute this part of row-column product
        for (int k = 0; k < BlockDim; k++)
        {
            sum += s_a[threadIdx.y][k] * s_b[k][threadIdx.x];
        }
        __syncthreads();
    }

    // all done; store computed element in matrix c
    if (row < n && col < n) c[IDX(row,col,n)] = sum;
}

//-----------------------------------------------------------------------------

// Check CUDA function return error code
void cudaChkErr(cudaError_t code)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr, "CUDA ERROR: %s\n", cudaGetErrorString(code));
        exit(EXIT_FAILURE);
    }
}

//----------------------------------------------------------------------------

// Fill matrix with reasonable values
void initializeMatrix(FLOAT* a, int m, int n, FLOAT sf = 1.0)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[IDX(i,j,n)] = sf * (-1.0 * i + j);
            //a[IDX(i,j,n)] = sf * ((i+1) * 10 + j+1);
        }
    }
}

//----------------------------------------------------------------------------

// Display matrix contents
void dumpMatrix(FLOAT* a, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        printf("[");
        for (int j = 0; j < n; j++)
        {
            printf(" %8.2f", a[IDX(i,j,n)]);
        }
        printf("]\n");
    }
}

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    double t0, t1;       // timing variables

    // Read matrix dimension from command line
    int n = 4;
    if (argc > 1) n = atoi(argv[1]);
    if (n <= 0) n = 4; // safety check
    printf("Matrix-matrix product with %dx%d matrices.\n", n, n);
    printf("BlockDim = %d; Threads per block = %d\n", BlockDim,
           BlockDim * BlockDim);

    // Declare and allocate memory for matrices
    FLOAT* a = new FLOAT [n * n];
    FLOAT* b = new FLOAT [n * n];
    FLOAT* c = new FLOAT [n * n];  // C = A * B

    // Initialize and display matrices (if small enough)
    initializeMatrix(a, n, n, 1.0);
    initializeMatrix(b, n, n, 2.0);
    if (n <= MaxSizeToDisplay)
    {
        printf("A =\n");
        dumpMatrix(a, n, n);
        printf("\nB =\n");
        dumpMatrix(b, n, n);
    }

    // Declare and allocate memory for matrices on device
    size_t matrixSize = n * n * sizeof(FLOAT);
    FLOAT* d_a;  // device memory for first factor
    FLOAT* d_b;  // device memory for second factor
    FLOAT* d_c;  // device memory for product
    cudaChkErr(cudaMalloc(&d_a, matrixSize));
    cudaChkErr(cudaMalloc(&d_b, matrixSize));
    cudaChkErr(cudaMalloc(&d_c, matrixSize));

    // Initialize matrices on device
    t0 = wtime();
    cudaChkErr(cudaMemcpy(d_a, a, matrixSize, cudaMemcpyHostToDevice));
    cudaChkErr(cudaMemcpy(d_b, b, matrixSize, cudaMemcpyHostToDevice));
    t1 = wtime();
    double data_transfer_time = t1 - t0;

    // Prepare for kernel launches: use 2D grid
    dim3 blockDim(BlockDim, BlockDim);
    dim3 gridDim((n + blockDim.x - 1) / blockDim.x,
                 (n + blockDim.y - 1) / blockDim.y); 

    // Compute product using global-memory-only kernel
    t0 = wtime();
    matmulGlobal<<<gridDim, blockDim>>>(d_c, d_a, d_b, n);
    cudaChkErr(cudaDeviceSynchronize()); // wait for kernel to finish
    cudaChkErr(cudaGetLastError());      // check for any errors in kernel
    t1 = wtime();
    double global_kernel_time = t1 - t0;

    // Copy result from device to host
    t0 = wtime();
    cudaChkErr(cudaMemcpy(c, d_c, matrixSize, cudaMemcpyDeviceToHost));
    t1 = wtime();
    data_transfer_time += (t1 - t0);
    if (n <= MaxSizeToDisplay)
    {
        printf("\n(Global Memory Only) A*B =\n");
        dumpMatrix(c, n, n);
    }

    // Compute product using shared-memory kernel
    t0 = wtime();
    matmulShared<<<gridDim, blockDim>>>(d_c, d_a, d_b, n);
    cudaChkErr(cudaDeviceSynchronize()); // wait for kernel to finish
    cudaChkErr(cudaGetLastError());      // check for any errors in kernel
    t1 = wtime();
    double shared_kernel_time = t1 - t0;

    // Copy result from device to host
    cudaChkErr(cudaMemcpy(c, d_c, matrixSize, cudaMemcpyDeviceToHost));
    if (n <= MaxSizeToDisplay)
    {
        printf("\n(with shared memory) A*B =\n");
        dumpMatrix(c, n, n);
    }

    // Report times and speedup
    printf("Data transfer time = %f sec\n", data_transfer_time);
    printf("Global kernel time = %f sec\n", global_kernel_time);
    printf("Shared kernel time = %f sec\n", shared_kernel_time);
    printf("Speedup = %6.2f\n", global_kernel_time / shared_kernel_time);

    // All done; "let my people go!"
    cudaChkErr(cudaFree(d_a));
    cudaChkErr(cudaFree(d_b));
    cudaChkErr(cudaFree(d_c));
    delete [] a;
    delete [] b;
    delete [] c;

    return 0;
}

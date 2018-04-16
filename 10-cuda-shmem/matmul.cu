// $Smake: nvcc -Xptxas -v -arch=sm_30 -O2 -o %F %f wtime.c

#include <cstdio>
#include <cuda.h>
#include "wtime.h"

#define IDX(i,j,n) ((i)*(n)+j)

#if !defined(BS)
const int BlockSize = 16;
#else
const int BlockSize = BS;  // normally 32 or less
#endif

const int MaxSizeToDisplay = 25;

typedef float FLOAT;
//typedef double FLOAT;

//-----------------------------------------------------------------------------

void cudaChkErr(
    cudaError_t code,  // value returned by CUDA runtime function
    int tag = -1     // optional tag; used to help identify call with error
    )
//
// Checks code returned by CUDA runtime function.  If not success, an error
// message is printed.  An optional second parameter is also printed if
// non-negative -- this can be used to help identify which function call
// was responsible for the error.
//
{
    if ( code != cudaSuccess )
    {
        fprintf( stderr, "CUDA ERROR: %s\n", cudaGetErrorString( code ) );
        if ( tag >= 0 ) fprintf( stderr, "tag = %d\n", tag );
        exit( EXIT_FAILURE );
    }
}

//----------------------------------------------------------------------------

// matrix-matrix kernel (global memory)
__global__ void matmulGlobal( FLOAT* c, FLOAT* a, FLOAT* b, int n )
{
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;
    if ( col < n && row < n )
    {
	FLOAT sum = 0.0;
	for ( int k = 0; k < n; k++ )
	{
	    sum += a[IDX(row,k,n)] * b[IDX(k,col,n)];
	}
	c[IDX(row,col,n)] = sum;
    }
}

//----------------------------------------------------------------------------

// matrix-matrix kernel (shared memory)
__global__ void matmulShared( FLOAT* c, FLOAT* a, FLOAT* b, int n )
{
    // element of matrix c to compute
    int col = blockIdx.x * blockDim.x + threadIdx.x;
    int row = blockIdx.y * blockDim.y + threadIdx.y;

    // loop over blocks from block row of matrix a and
    // block column of matrix b.
    FLOAT sum = 0.0;
    int numBlocks = ( n + BlockSize - 1 ) / BlockSize;
    for ( int m = 0; m < numBlocks; m++ )
    {
	// copy block from matrix to shared memory
	__shared__ FLOAT a_s[BlockSize][BlockSize];
	__shared__ FLOAT b_s[BlockSize][BlockSize];
	int c = m * BlockSize + threadIdx.x;
	int r = m * BlockSize + threadIdx.y;
	a_s[threadIdx.y][threadIdx.x] = a[IDX(row,c,n)];
	b_s[threadIdx.y][threadIdx.x] = b[IDX(r,col,n)];
	__syncthreads();

	// length of this part of row-column product is BlockSize
	// except for last block when it may be smaller
	int sliceLen = ( m == numBlocks - 1 ? n - m * BlockSize : BlockSize );

	// compute this part of row-column product
	for ( int k = 0; k < sliceLen; k++ )
	{
	    sum += a_s[threadIdx.y][k] * b_s[k][threadIdx.x];
	}
	__syncthreads();
    }

    // all done; store computed element in matrix c
    if ( col < n && row < n ) c[IDX(row,col,n)] = sum;
}

//----------------------------------------------------------------------------

void initializeMatrix( FLOAT* a, int m, int n, double k )
{
    for ( int i = 0; i < m; i++ )
    {
	for ( int j = 0; j < n; j++ )
	{
	    a[IDX(i,j,n)] = k * ( -1.0 * i + j );// / ( n * m );
	}
    }
}

//----------------------------------------------------------------------------

void dumpMatrix( FLOAT* a, int m, int n )
{
    for ( int i = 0; i < m; i++ )
    {
	printf( "[" );
	for ( int j = 0; j < n; j++ )
	{
	    printf( " %8.2f", a[IDX(i,j,n)] );
	}
	printf( "]\n" );
    }
}

//----------------------------------------------------------------------------

int main( int argc, char* argv[] )
{
    double t0, t1;

    int n = 4;
    if ( argc > 1 ) n = atoi( argv[1] );
    if ( n <= 0 ) n = 4; // safety check
    printf( "matrix-matrix product with %dx%d matrices.\n", n, n );

    // Declare and allocate memory for matrices
    FLOAT* a = new FLOAT [n * n];
    FLOAT* b = new FLOAT [n * n];
    FLOAT* c = new FLOAT [n * n];  // C = A * B

    // Initialize and display matrices (if small enough)
    initializeMatrix( a, n, n, 0.1 );
    initializeMatrix( b, n, n, 0.01 );
    if ( n <= MaxSizeToDisplay )
    {
	printf( "A =\n" );
	dumpMatrix( a, n, n );
	printf( "\nB =\n" );
	dumpMatrix( b, n, n );
    }

    // Declare and allocate memory for matrices on device
    size_t matrixSize = n * n * sizeof( FLOAT );
    FLOAT* a_d;  // device memory for first factor
    FLOAT* b_d;  // device memory for second factor
    FLOAT* c_d;  // device memory for product
    cudaChkErr( cudaMalloc( (void**) &a_d, matrixSize ) );
    cudaChkErr( cudaMalloc( (void**) &b_d, matrixSize ) );
    cudaChkErr( cudaMalloc( (void**) &c_d, matrixSize ) );

    // Initialize matrices on device
    cudaChkErr( cudaMemcpy( a_d, a, matrixSize, cudaMemcpyHostToDevice ) );
    cudaChkErr( cudaMemcpy( b_d, b, matrixSize, cudaMemcpyHostToDevice ) );

    // Set up CUDA events for timing
    cudaEvent_t event0, event1;
    cudaChkErr( cudaEventCreate( &event0 ) );
    cudaChkErr( cudaEventCreate( &event1 ) );

    // Prepare for kernel launches: use 2D grid
    dim3 blockDim( BlockSize, BlockSize );
    dim3 gridDim( ( n + blockDim.x - 1 ) / blockDim.x,
		  ( n + blockDim.y - 1 ) / blockDim.y ); 

    // Compute product using global-memory-only kernel
    t0 = wtime();
    cudaEventRecord( event0, 0 );
    matmulGlobal<<<gridDim, blockDim>>>( c_d, a_d, b_d, n );
    cudaChkErr( cudaDeviceSynchronize() );
    cudaChkErr( cudaGetLastError() );
    cudaChkErr( cudaEventRecord( event1, 0 ) );
    cudaChkErr( cudaEventSynchronize( event1 ) );// wait for event 1 to finish
    t1 = wtime();

    cudaChkErr( cudaMemcpy( c, c_d, matrixSize, cudaMemcpyDeviceToHost ) );
    if ( n <= MaxSizeToDisplay )
    {
	printf( "\n(Global Memory Only) A*B =\n" );
	dumpMatrix( c, n, n );
    }

    // Report times
    float global_time_ms;
    cudaChkErr( cudaEventElapsedTime( &global_time_ms, event0, event1 ) );
    double global_wall_time = t1 - t0;
    printf( "Global kernel time = %e sec, elapsed wall time = %e sec\n",
	    global_time_ms / 1000.0, global_wall_time );

    // Compute product using shared-memory kernel
    t0 = wtime();
    cudaChkErr( cudaEventRecord( event0, 0 ) );
    matmulShared<<<gridDim, blockDim>>>( c_d, a_d, b_d, n );
    cudaChkErr( cudaDeviceSynchronize() );
    cudaChkErr( cudaGetLastError() );
    cudaChkErr( cudaEventRecord( event1, 0 ) );
    cudaChkErr( cudaEventSynchronize( event1 ) );// wait for event 1 to finish
    t1 = wtime();

    cudaChkErr( cudaMemcpy( c, c_d, matrixSize, cudaMemcpyDeviceToHost ) );
    if ( n <= MaxSizeToDisplay )
    {
	printf( "\n(with shared memory) A*B =\n" );
	dumpMatrix( c, n, n );
    }

    double sum = 0.0;
    for ( int i = 0; i<n*n; i++ )
        sum += c[i] / double( n * n );
    printf( "sum = %f\n", sum );
    
    // Report times and speedup
    float shared_time_ms;
    cudaChkErr( cudaEventElapsedTime( &shared_time_ms, event0, event1 ) );
    double shared_wall_time = t1 - t0;
    printf( "Shared kernel time = %e sec, elapsed wall time = %e sec\n",
	    shared_time_ms / 1000.0, shared_wall_time );
    printf( "Device speedup = %6.2f, Wall clock speedup = %6.2f\n",
	    global_time_ms / shared_time_ms,
	    global_wall_time / shared_wall_time );

    // all done; "let my people go!"
    cudaChkErr( cudaFree( a_d ) );
    cudaChkErr( cudaFree( b_d ) );
    cudaChkErr( cudaFree( c_d ) );
    delete [] a;
    delete [] b;
    delete [] c;

    return 0;
}

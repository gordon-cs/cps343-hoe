// $Smake: nvcc -Xptxas -v -arch=sm_30 -O2 -o %F %f wtime.c

#include <cstdio>
#include <cuda.h>
#include "wtime.h"

#define IDX(i,j,n) ((i)*(n)+j)

#if !defined(BS)
const int BlockSize = 16;
#else
const int BlockSize = BS;  // normally 64 or less
#endif

const int MaxSizeToDisplay = 25;

typedef float FLOAT;
//typedef double FLOAT;

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

void initializeMatrix( FLOAT* a, int m, int n )
{
    for ( int i = 0; i < m; i++ )
    {
	for ( int j = 0; j < n; j++ )
	{
	    a[IDX(i,j,n)] = -1.0 * i + j;
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
    initializeMatrix( a, n, n );
    initializeMatrix( b, n, n );
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
    cudaMalloc( (void**) &a_d, matrixSize );
    cudaMalloc( (void**) &b_d, matrixSize );
    cudaMalloc( (void**) &c_d, matrixSize );

    // Initialize matrices on device
    cudaMemcpy( a_d, a, matrixSize, cudaMemcpyHostToDevice );
    cudaMemcpy( b_d, b, matrixSize, cudaMemcpyHostToDevice );

    // Set up CUDA events for timing
    cudaEvent_t event0, event1;
    cudaEventCreate( &event0 );
    cudaEventCreate( &event1 );

    // Prepare for kernel launches: use 2D grid
    dim3 blockDim( BlockSize, BlockSize );
    dim3 gridDim( ( n + blockDim.x - 1 ) / blockDim.x,
		  ( n + blockDim.y - 1 ) / blockDim.y ); 

    // Compute product using global-memory-only kernel
    t0 = wtime();
    cudaEventRecord( event0, 0 );
    matmulGlobal<<<gridDim, blockDim>>>( c_d, a_d, b_d, n );
    cudaEventRecord( event1, 0 );
    cudaEventSynchronize( event1 );  // wait for event 1 to complete
    t1 = wtime();

    cudaMemcpy( c, c_d, matrixSize, cudaMemcpyDeviceToHost );
    if ( n <= MaxSizeToDisplay )
    {
	printf( "\n(Global Memory Only) A*B =\n" );
	dumpMatrix( c, n, n );
    }

    // Report times
    float global_time_ms;
    cudaEventElapsedTime( &global_time_ms, event0, event1 );
    double global_wall_time = t1 - t0;
    printf( "Global kernel time = %e sec, elapsed wall time = %e sec\n",
	    global_time_ms / 1000.0, global_wall_time );

    // Compute product using shared-memory kernel
    t0 = wtime();
    cudaEventRecord( event0, 0 );
    matmulShared<<<gridDim, blockDim>>>( c_d, a_d, b_d, n );
    cudaEventRecord( event1, 0 );
    cudaEventSynchronize( event1 );  // wait for event 1 to complete
    t1 = wtime();

    cudaMemcpy( c, c_d, matrixSize, cudaMemcpyDeviceToHost );
    if ( n <= MaxSizeToDisplay )
    {
	printf( "\n(with shared memory) A*B =\n" );
	dumpMatrix( c, n, n );
    }

    // Report times and speedup
    float shared_time_ms;
    cudaEventElapsedTime( &shared_time_ms, event0, event1 );
    double shared_wall_time = t1 - t0;
    printf( "Shared kernel time = %e sec, elapsed wall time = %e sec\n",
	    shared_time_ms / 1000.0, shared_wall_time );
    printf( "Device speedup = %6.2f, Wall clock speedup = %6.2f\n",
	    global_time_ms / shared_time_ms,
	    global_wall_time / shared_wall_time );

    // all done; "let my people go!"
    cudaFree( a_d );
    cudaFree( b_d );
    cudaFree( c_d );
    delete [] a;
    delete [] b;
    delete [] c;

    return 0;
}

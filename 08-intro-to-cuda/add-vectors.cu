// $Smake: nvcc -O2 -o %F %f
//
// add-vectors.cu - addition of two arrays on GPU device
//
// This program follows a very standard pattern:
//  1) allocate memory on host
//  2) allocate memory on device
//  3) initialize memory on host
//  4) copy memory from host to device
//  5) execute kernel(s) on device
//  6) copy result(s) from device to host
//
// Note: it may be possible to initialize memory directly on the device,
// in which case steps 3 and 4 are not necessary, and step 1 is only
// necessary to allocate memory to hold results.

#include <stdio.h>
#include <cuda.h>

//-----------------------------------------------------------------------------
// Kernel that executes on CUDA device

__global__ void add_vectors(
    float *c,      // out - pointer to result vector c
    float *a,      // in  - pointer to summand vector a
    float *b,      // in  - pointer to summand vector b
    int n          // in  - vector length
    )
{
    // Assume single block grid and 1-D block
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // Only do calculation if we have real data to work with
    if (idx < n) c[idx] = a[idx] + b[idx];
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
// Main program executes on host device

int main(int argc, char* argv[])
{
    // determine vector length
    int n = 10;      // set default length
    if (argc > 1)
    {
        n = atoi(argv[1]);  // override default length
        if (n <= 0)
        {
            fprintf(stderr, "Vector length must be positive\n");
            return EXIT_FAILURE;
        }
    }

    // determine vector size in bytes
    const size_t vector_size = n * sizeof(float);

    // declare pointers to vectors in host memory and allocate memory
    float *a, *b, *c;
    a = (float*) malloc(vector_size);
    b = (float*) malloc(vector_size);
    c = (float*) malloc(vector_size);

    // declare pointers to vectors in device memory and allocate memory
    float *d_a, *d_b, *d_c;
    cudaMalloc((void**) &d_a, vector_size);
    cudaMalloc((void**) &d_b, vector_size);
    cudaMalloc((void**) &d_c, vector_size);

    // initialize vectors and copy them to device
    for (int i = 0; i < n; i++)
    {
        a[i] =   1.0 * i;
        b[i] = 100.0 * i;
    }
    cudaMemcpy(d_a, a, vector_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b, vector_size, cudaMemcpyHostToDevice);

    // do calculation on device
    int block_size = 16;
    int num_blocks = (n - 1 + block_size) / block_size;
    add_vectors<<<num_blocks, block_size>>>(d_c, d_a, d_b, n);

    // retrieve result from device and store on host
    cudaMemcpy(c, d_c, vector_size, cudaMemcpyDeviceToHost);

    // print results for vectors up to length 100
    if (n <= 100)
    {
        for (int i = 0; i < n; i++)
        {
            printf("%8.2f + %8.2f = %8.2f\n", a[i], b[i], c[i]);
        }
    }

    // cleanup and quit
    cudaFree(d_a);
    cudaFree(d_b);
    cudaFree(d_c);
    free(a);
    free(b);
    free(c);
  
    return 0;
}

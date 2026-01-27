/*
 * $Smake: g++ -Wall -O3 -o %F %f
 *
 * Matrix-matrix product
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <unistd.h>
#include <time.h>

// Macro to index matrices in column-major (Fortran) order
#define IDX(i,j,stride) ((i)+(j)*(stride))  // column major

//----------------------------------------------------------------------------
// Usage

void usage(char* program_name)
{
    fprintf(stderr, "Usage: %s [-v]\n", program_name);
}

//----------------------------------------------------------------------------
// Dump Matrix

void dumpMatrix(double* a, int m, int n, int stride)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf(" %8.2f", a[IDX(i,j,stride)]);
        }
        printf("\n");
    }
    printf("\n");
}

//----------------------------------------------------------------------------
// create Matrix based on supplied name - either "A" or "B"

void createMatrix(const char* name, double** a, int& rows, int& cols)
{
    if (strcmp(name, "A") == 0)
    {
        rows = 4;
        cols = 2;
        // matrix data in column order
        double a_tmp[rows*cols] = {4.0, 2.0, -2.0, 1.0, -4.0, -1.0, -3.0, 4.0};
        *a = new double [rows * cols];
        memcpy(*a, a_tmp, sizeof(*a[0])*rows*cols);
    }
    else if (strcmp(name, "B") == 0)
    {
        rows = 2;
        cols = 3;
        // matrix data in column order
        double a_tmp[rows*cols] = {5.0, -3.0, -4.0, 1.0, 2.0, -3.0};
        *a = new double [rows * cols];
        memcpy(*a, a_tmp, sizeof(*a[0])*rows*cols);
    }
    else
    {
        rows = cols = 0;
        *a = NULL;
    }
}

//----------------------------------------------------------------------------
// Form matrix product C = AB

void matmat_jki(double* c, double* a, int nrow_a, int ncol_a,
                double* b, int ncol_b)
{
    const int nrow_b = ncol_a;
    const int nrow_c = nrow_a;
    for (int j = 0; j < ncol_b; j++)
    {
        for (int i = 0; i < nrow_a; i++) c[IDX(i,j,nrow_c)] = 0.0;
        for (int k = 0; k < ncol_a; k++)
        {
            for (int i = 0; i < nrow_a; i++)
            {
                c[IDX(i,j,nrow_c)] += a[IDX(i,k,nrow_a)] * b[IDX(k,j,nrow_b)];
            }
        }
    }
}

//----------------------------------------------------------------------------
// Main program

int main(int argc, char* argv[])
{
    double* a;               // left matrix
    double* b;               // right matrix
    double* c;               // product C = AB
    int nrow_a, ncol_a;      // dimensions of left matrix
    int nrow_b, ncol_b;      // dimensions of right matrix
    int nrow_c, ncol_c;      // dimensiosn of product matrix
    int verbose = 0;         // nonzero for extra output

    // Process command line
    int ch;
    while ((ch = getopt(argc, argv, "v")) != -1)
    {
        switch (ch)
        {
            case 'v':
                verbose++;
                break;
            default:
                usage(argv[0]);
                return EXIT_FAILURE;
        }
    }
    argv[optind - 1] = argv[0];
    argv += (optind - 1);
    argc -= (optind - 1);

    // Make sure there are no additional arguments
    if (argc != 1)
    {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    // Create matrix data and optionally display it
    createMatrix("A", &a, nrow_a, ncol_a);
    createMatrix("B", &b, nrow_b, ncol_b);

    if (ncol_a != nrow_b)
    {
        fprintf(stderr, "Error: matrix dimensions are not compatible\n");
        return EXIT_FAILURE;
    }

    if (verbose)
    {
        printf("Matrix A:\n");
        dumpMatrix(a, nrow_a, ncol_a, nrow_a);
        printf("Matrix B:\n");
        dumpMatrix(b, nrow_b, ncol_b, nrow_b);
    }

    // Compute matrix product C = AB and display it
    nrow_c = nrow_a;
    ncol_c = ncol_b;
    c = new double [nrow_c * ncol_c];

    matmat_jki(c, a, nrow_a, ncol_a, b, ncol_b);

    printf("Matrix C = AB:\n");
    dumpMatrix(c, nrow_c, ncol_c, nrow_c);

    // Clean up and quit
    delete [] a;
    delete [] b;
    delete [] c;
    return 0;
}

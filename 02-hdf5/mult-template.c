/*
 * $Smake: gcc -Wall -O3 -o %F %f
 *
 * Computes a matrix-matrix product
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

/* Macro to index matrices in column-major (Fortran) order */
#define IDX(i,j,stride) ((i)+(j)*(stride))  /* column major */

/*----------------------------------------------------------------------------
 * Display string showing how to run program from command line
 *
 * Input:
 *   char* program_name (in)  name of executable
 * Output:
 *   writes to stderr
 * Returns:
 *   nothing
 */
void usage(char* program_name)
{
    fprintf(stderr, "Usage: %s [-v]\n", program_name);
}

/*----------------------------------------------------------------------------
 * Dump Matrix
 *
 * Parameters:
 *   double* a          (in)  pointer to matrix data
 *   int rows           (in)  number of rows in matrix
 *   int cols           (in)  number of columns in matrix
 *   int stride         (in)  =rows if column major or =cols if row major
* Returns:
 *   nothing
 */
void dumpMatrix(double* a, int rows, int cols, int stride)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            printf(" %8.2f", a[IDX(i,j,stride)]);
        }
        printf("\n");
    }
    printf("\n");
}

/*----------------------------------------------------------------------------
 * Create Matrix based on supplied name
 *
 * Parameters:
 *   char* name         (in)  name of matrix ("A" or "B")
 *   double** a         (out) pointer to pointer to matrix data
 *   int* rows          (out) pointer to number of rows
 *   int* cols          (out) pointer to number of cols
 * Returns:
 *   nothing
 */
void createMatrix(char* name, double** a, int* rows, int* cols)
{
    if (strcmp(name, "A") == 0)
    {
        *rows = 4;
        *cols = 2;
        *a = (double*) malloc(*rows * *cols * sizeof(double));
        (*a)[IDX(0,0,*rows)] =  4.0;
        (*a)[IDX(1,0,*rows)] =  2.0;
        (*a)[IDX(2,0,*rows)] = -2.0;
        (*a)[IDX(3,0,*rows)] =  1.0;
        (*a)[IDX(0,1,*rows)] = -4.0;
        (*a)[IDX(1,1,*rows)] = -1.0;
        (*a)[IDX(2,1,*rows)] = -3.0;
        (*a)[IDX(3,1,*rows)] =  4.0;
    }
    else if (strcmp(name, "B") == 0)
    {
        *rows = 2;
        *cols = 3;
        *a = (double*) malloc(*rows * *cols * sizeof(double));
        (*a)[IDX(0,0,*rows)] =  5.0;
        (*a)[IDX(1,0,*rows)] = -3.0;
        (*a)[IDX(0,1,*rows)] = -4.0;
        (*a)[IDX(1,1,*rows)] =  1.0;
        (*a)[IDX(0,2,*rows)] =  2.0;
        (*a)[IDX(1,2,*rows)] = -3.0;
    }
}

/*----------------------------------------------------------------------------
 * Form matrix product C = AB
 *
 * Parameters:
 *   double* c          (out) pointer to result matrix (nrows_a x ncols_b)
 *   double* a          (in)  pointer to left matrix
 *   int nrow_a         (in)  rows in left matrix
 *   int ncol_a         (in)  cols in left matrix (rows in right matrix)
 *   double* b          (in)  pointer to right matrix
 *   int ncol_b         (in)  cols in right matrix
 * Returns:
 *   nothing
 */
void matmat_jki(double* c, double* a, int nrow_a, int ncol_a,
                 double* b, int ncol_b)
{
    const int nrow_b = ncol_a;
    const int nrow_c = nrow_a;
    for (int j = 0; j < ncol_b; j++)
    {
        for (int i = 0; i < nrow_a; i++) c[IDX(i,j,nrow_c)] = 0.0;
        for (int k = 0; k < ncol_a; k++)
            for (int i = 0; i < nrow_a; i++)
                c[IDX(i,j,nrow_c)] += a[IDX(i,k,nrow_a)] * b[IDX(k,j,nrow_b)];
    }
}

/*----------------------------------------------------------------------------
 * Main program
 */
int main(int argc, char* argv[])
{
    double* a;             /* left matrix */
    double* b;             /* right matrix */
    double* c;             /* product C = AB */
    int nrow_a, ncol_a;    /* dimensions of left matrix */
    int nrow_b, ncol_b;    /* dimensions of right matrix */
    int nrow_c, ncol_c;    /* dimensions of product matrix */
    int verbose = 0;       /* nonzero for extra output */

    /* Process command line */
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

    /* Make sure there are no additional arguments */
    if (argc != 1)
    {
        usage(argv[0]);
        return EXIT_FAILURE;
    }

    /* Create matrix data and optionally display it */
    createMatrix("A", &a, &nrow_a, &ncol_a);
    createMatrix("B", &b, &nrow_b, &ncol_b);

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

    /* Compute matrix product C = AB and display it */
    nrow_c = nrow_a;
    ncol_c = ncol_b;
    c = (double*) malloc(nrow_c * ncol_c * sizeof(double));

    matmat_jki(c, a, nrow_a, ncol_a, b, ncol_b);

    printf("Matrix C = AB:\n");
    dumpMatrix(c, nrow_c, ncol_c, nrow_c);

    /* Clean up and quit */
    free(a);
    free(b);
    free(c);
    return 0;
}

/*
 * $Smake: gcc -O3 -Wall -pedantic -o %F %f -lf77blas -lsatlas -lgfortran
 *
 * Jonathan Senning <jonathan.senning@gordon.edu>
 * Department of Mathematics and Computer Science
 * Gordon College, 255 Grapevine Road, Wenham MA 01984-1899
 *
 * This program is released into the public domain.
 *
 * Reports times for computing C : = A * B using various methods where
 *     A is MxP
 *     B is PxN
 *     C is MxN
 * There are four versions:
 *     Naive            - standard "textbook" matrix product.
 *     Standard tiling  - product using predetermined tiles
 *     Recursive tiling - product using recursively determined matrix tiles
 *     BLAS             - uses BLAS routine DGEMM
 *     CBLAS            - uses CBLAS routine cblas_dgemm
 * All versions except BLAS and CBLAS can be carried out with a particular
 * one of the "ijk" loop-order forms.
 *
 * Running the naive version:
 *     matrix_prod -n <ijk> <M> <P> <N>
 *   Example:
 *     matrix_prod -n ikj 1000 500 1000
 *
 * Running the standard tiling version:
 *     matrix_prod -t <ijk> <M> <P> <N> <TILE_M> <TILE_P> <TILE_N>
 *   Example:
 *     matrix_prod -t ijk 1000 500 1000 64 32 64
 *
 * Running the recursive tiling version:
 *     matrix_prod -r <ijk> <M> <P> <N> <MAX_TILESIZE>
 *     (Note that MAX_TILESIZE is the total number of elements in the
 *     block, not the block dimension.)
 *   Example:
 *     matrix_prod -t jki 1000 500 1000 12000
 *
 * Running the BLAS version:
 *     matrix_prod -b <M> <P> <N>
 *   Example:
 *     matrix_prod -b 1000 500 1000
 *
 * Running the CBLAS version:
 *     matrix_prod -c <M> <P> <N>
 *   Example:
 *     matrix_prod -c 1000 500 1000
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>

enum {
    MODE_UNKNOWN, MODE_BLAS, MODE_CBLAS, MODE_NAIVE, MODE_RECURSIVE_TILE,
    MODE_STD_TILE
};

void (*do_product)(double**, double**, double**, int, int,
                   int, int, int, int) = NULL;

#if defined(HAVE_BLAS)
/* Prototype of Fortran BLAS routine DGEMM */
void dgemm_(char* transa, char* transb, int* m, int* n, int* k,
            double* alpha, double* a, int* lda, double* b, int* ldb,
            double* beta, double* c, int* ldc);

/*---------------------------------------------------------------------------
 *
 * C interface to Fortran BLAS routine DGEMM.  See DGEMM man page for
 * parameter description.
 */
void dgemm(char transa, char transb, int m, int n, int k, double alpha,
           double* a, int lda, double* b, int ldb, double beta, double* c,
           int ldc)
{
    dgemm_(&transa, &transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta,
           c, &ldc);
}
#endif

#if defined(HAVE_CBLAS)
#include <cblas.h>
#endif

#if defined(HAVE_GSL_CBLAS)
#include <gsl/gsl_cblas.h>
#endif

/*---------------------------------------------------------------------------
 *
 * Returns the number of seconds since some fixed arbitrary time in the past
 */
double wtime(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double) (ts.tv_sec +ts.tv_nsec / 1.0e9);
}

/*---------------------------------------------------------------------------
 *
 * Allocate matrix
 *
 * Input:
 *   int rows, cols   - rows and columns in matrix
 *
 * Output:
 *   double** a       - pointer to array of pointers to rows of matrix
 */
double** allocateMatrix(int rows, int cols)
{
    int i;
    double** a = (double**) malloc(rows * sizeof(double*));
    a[0] = (double*) malloc(rows * cols * sizeof(double));
    for (i = 1; i < rows; i++)
    {
        a[i] = &a[0][i * cols];
    }
    return a;
}

/*---------------------------------------------------------------------------
 *
 * Deallocate matrix memory
 *
 * Input:
 *   double** a       - matrix
 *
 * Output:
 *   none
 */
void deallocateMatrix(double** a)
{
    free(a[0]);
    free(a);
}

/*---------------------------------------------------------------------------
 *
 * Computes sum of matrix elements
 *
 * Input:
 *   double** a       - matrix
 *   int rows, cols   - rows and columns in matrix
 *
 * Output:
 *   double sum       - sum of matrix elements
 */
double checksum(double** a, int rows, int cols)
{
    int i, j;
    double sum = 0.0;
    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            sum += a[i][j];
    return sum;
}

/*---------------------------------------------------------------------------
 *
 * Initialize matrices A and B with random numbers and set C to zero
 *
 * Input
 *   double** a      - first matrix used in product
 *   double** b      - second matrix used in product
 *   double** c      - matrix to be updated by product
 *   int rows        - number of rows in A and C
 *   int cols        - number of columns in B and C
 *   int mids        - number of columns in A and rows in B
 *   int verbosity   - program verification: verbosity > 1 --> use fixed seed
 *
 * Output
 *   double** a      - filled with uniform random numbers on [0,1)
 *   double** b      - filled with uniform random numbers on [0,1)
 *   double** c      - filled with zero
 */
void initialize_matrices(double** a, double** b, double** c,
                         int rows, int cols, int mids, int verbosity)
{
    int i, j;
    const int REF_SEED = 42;

    srandom((unsigned int) (verbosity > 1 ? REF_SEED : time(NULL)));
    for (i = 0; i < rows; i++)
        for (j = 0; j < mids; j++)
            a[i][j] = (double) random() / RAND_MAX;

    for (i = 0; i < mids; i++)
        for (j = 0; j < cols; j++)
            b[i][j] = (double) random() / RAND_MAX;

    for (i = 0; i < rows; i++)
        for (j = 0; j < cols; j++)
            c[i][j] = 0.0;
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix operation C = A * B + C using ijk loop order
 * 
 * Input:
 *   double** a              - first factor of product
 *   double** b              - second factor of product
 *   double** c              - matrix to update
 *   int row_start, row_end  - first and last row of block of A and C
 *   int col_start, col_end  - first and last column of block of B and C
 *   int mid_start, mid_end  - first and last col of A and row of B.
 *
 * Output:
 *   double** c       - updated matrix
 */
void do_ijk_product(double** a, double** b, double** c,
                    int row_start, int row_end,
                    int col_start, int col_end,
                    int mid_start, int mid_end)
{
    int i, j, k;
    for (i = row_start; i <= row_end; i++)
        for (j = col_start; j <= col_end; j++)
            for (k = mid_start; k <= mid_end; k++)
                c[i][j] += a[i][k] * b[k][j];
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix operation C = A * B + C using jik loop order
 * 
 * Input:
 *   double** a              - first factor of product
 *   double** b              - second factor of product
 *   double** c              - matrix to update
 *   int row_start, row_end  - first and last row of block of A and C
 *   int col_start, col_end  - first and last column of block of B and C
 *   int mid_start, mid_end  - first and last col of A and row of B.
 *
 * Output:
 *   double** c       - updated matrix
 */
void do_jik_product(double** a, double** b, double** c,
                    int row_start, int row_end,
                    int col_start, int col_end,
                    int mid_start, int mid_end)
{
    int i, j, k;
    for (j = col_start; j <= col_end; j++)
        for (i = row_start; i <= row_end; i++)
            for (k = mid_start; k <= mid_end; k++)
                c[i][j] += a[i][k] * b[k][j];
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix operation C = A * B + C using ikj loop order
 * 
 * Input:
 *   double** a              - first factor of product
 *   double** b              - second factor of product
 *   double** c              - matrix to update
 *   int row_start, row_end  - first and last row of block of A and C
 *   int col_start, col_end  - first and last column of block of B and C
 *   int mid_start, mid_end  - first and last col of A and row of B.
 *
 * Output:
 *   double** c       - updated matrix
 */
void do_ikj_product(double** a, double** b, double** c,
                    int row_start, int row_end,
                    int col_start, int col_end,
                    int mid_start, int mid_end)
{
    int i, j, k;
    for (i = row_start; i <= row_end; i++)
        for (k = mid_start; k <= mid_end; k++)
            for (j = col_start; j <= col_end; j++)
                c[i][j] += a[i][k] * b[k][j];
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix operation C = A * B + C using jki loop order
 * 
 * Input:
 *   double** a              - first factor of product
 *   double** b              - second factor of product
 *   double** c              - matrix to update
 *   int row_start, row_end  - first and last row of block of A and C
 *   int col_start, col_end  - first and last column of block of B and C
 *   int mid_start, mid_end  - first and last col of A and row of B.
 *
 * Output:
 *   double** c       - updated matrix
 */
void do_jki_product(double** a, double** b, double** c,
                    int row_start, int row_end,
                    int col_start, int col_end,
                    int mid_start, int mid_end)
{
    int i, j, k;
    for (j = col_start; j <= col_end; j++)
        for (k = mid_start; k <= mid_end; k++)
            for (i = row_start; i <= row_end; i++)
                c[i][j] += a[i][k] * b[k][j];
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix operation C = A * B + C using kij loop order
 * 
 * Input:
 *   double** a              - first factor of product
 *   double** b              - second factor of product
 *   double** c              - matrix to update
 *   int row_start, row_end  - first and last row of block of A and C
 *   int col_start, col_end  - first and last column of block of B and C
 *   int mid_start, mid_end  - first and last col of A and row of B.
 *
 * Output:
 *   double** c       - updated matrix
 */
void do_kij_product(double** a, double** b, double** c,
                    int row_start, int row_end,
                    int col_start, int col_end,
                    int mid_start, int mid_end)
{
    int i, j, k;
    for (k = mid_start; k <= mid_end; k++)
        for (i = row_start; i <= row_end; i++)
            for (j = col_start; j <= col_end; j++)
                c[i][j] += a[i][k] * b[k][j];
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix operation C = A * B + C using kji loop order
 * 
 * Input:
 *   double** a              - first factor of product
 *   double** b              - second factor of product
 *   double** c              - matrix to update
 *   int row_start, row_end  - first and last row of block of A and C
 *   int col_start, col_end  - first and last column of block of B and C
 *   int mid_start, mid_end  - first and last col of A and row of B.
 *
 * Output:
 *   double** c       - updated matrix
 */
void do_kji_product(double** a, double** b, double** c,
                    int row_start, int row_end,
                    int col_start, int col_end,
                    int mid_start, int mid_end)
{
    int i, j, k;
    for (k = mid_start; k <= mid_end; k++)
        for (j = col_start; j <= col_end; j++)
            for (i = row_start; i <= row_end; i++)
                c[i][j] += a[i][k] * b[k][j];
}

/*---------------------------------------------------------------------------
 *
 * Computes block-oriented matrix-matrix product recursively.
 *
 * Input:
 *   double** c       - matrix product C = A * B
 *   double** a       - first factor of product
 *   double** b       - second factor of product
 *   int crow, ccol   - starting row and column of block of C
 *   int arow, acol   - starting row and column of block of A
 *   int brow, bcol   - starting row and column of block of B
 *   int l, m, n      - dims of blocks: A is l x m, B is m x n, C is l x n
 *   int N            - full row length (column dimension) of matrix B
 *   int threshold    - B blocks larger than this are partitioned
 *
 * Output:
 *   double** c       - matrix product C = A * B
 *
 * Algorithm based on one presented on page 276 of "Parallel Programming in
 * C with MPI and OpenMP", Michael J. Quinn, McGraw-Hill, 2004.
 *
 * **** NOTE ****: There is a typo in Quinn's code in the recursive call
 * to mm_rec().  The 5th parameter should be "ccol + nhalf[j]" and not
 * use "mhalf" as shown in the text.  The error only shows up when the
 * dimensions of the matrices are not uniform.
 */
void mm_rec(double** c, double** a, double** b,
            int crow, int ccol, int arow, int acol, int brow, int bcol,
            int l, int m, int n, int N, int threshold)
{
    int lhalf[3], mhalf[3], nhalf[3];
    int i, j, k;
    
    if (m * n > threshold)
    {
        lhalf[0] = 0; lhalf[1] = l/2; lhalf[2] = l - lhalf[1];
        mhalf[0] = 0; mhalf[1] = m/2; mhalf[2] = m - mhalf[1];
        nhalf[0] = 0; nhalf[1] = n/2; nhalf[2] = n - nhalf[1];
        for (i = 0; i < 2; i++)
        {
            for (j = 0; j < 2; j++)
            {
                for (k = 0; k < 2; k++)
                {
                    mm_rec(c, a, b,
                           crow + lhalf[i], ccol + nhalf[j],
                           arow + lhalf[i], acol + mhalf[k],
                           brow + mhalf[k], bcol + nhalf[j],
                           lhalf[i + 1], mhalf[k + 1], nhalf[j + 1],
                           N, threshold);
                }
            }
        }
    }
    else
    {
        do_product(a, b, c,
                   arow, arow + l - 1,
                   bcol, bcol + n - 1,
                   acol, acol + m - 1);
    }
}

#if defined(HAVE_BLAS)
/*---------------------------------------------------------------------------
 *
 * Compute matrix product using BLAS routine DGEMM.
 *
 * Input
 *   int argc        - length of argv[] array
 *   char* argv[]    - pointer to command line parameter array
 *   int verbosity   - program verification: verbosity > 0 gives more output
 *
 * Output
 *   double          - elapsed time for product computation
 */
double multiply_by_blas(int argc, char* argv[], int verbosity)
{
    int rows, cols, mids;
    double **a, **b, **c;
    double t1, t2;
    double sec;
    double gflop_count;

    /*
     * process command line arguments
     */
    rows = atoi(argv[0]);
    mids = atoi(argv[1]);
    cols = atoi(argv[2]);
    gflop_count = 2.0 * rows * mids * cols / 1.0e9;

    if (verbosity > 0)
    {
        printf("BLAS: rows = %d, mids = %d, columns = %d\n",
                rows, mids, cols);
    }

    /*
     * allocate and initialize matrices
     */
    a = (double**) allocateMatrix(rows, mids);
    b = (double**) allocateMatrix(mids, cols);
    c = (double**) allocateMatrix(rows, cols);
    initialize_matrices(a, b, c, rows, cols, mids, verbosity);

    /*
     * compute product: There is an implicit matrix transpose when
     * passing from Fortran to C and vice-versa.  To compute C :=
     * alpha * A * B + beta * C we use dgemm() to compute C' := alpha
     * * B' * A' + beta * C'.  The first two arguments to dgemm() are
     * 'N' indicating we don't want a transpose in addition to the
     * implicit one.  The matrices A and B are passed in reverse order
     * so dgemm() receives (after the implicit transpose) B' and A'.
     * Arguments 3 and 4 are the dimensions of C' and argument 5 is
     * the column dimension of B' (and the row dimension of A').
     */
    t1 = wtime();
    dgemm('N', 'N', cols, rows, mids, 1.0, &b[0][0], cols, &a[0][0], mids, 
          0.0, &c[0][0], cols);
    t2 = wtime();
    sec = t2 - t1;

    if (verbosity > 1)
        printf("checksum = %f\n", checksum(c, rows, cols));

    printf("BLAS:        %6.3f secs %6.3f GFlop/s (%5d x %5d x %5d)\n",
           sec, gflop_count / sec, rows, mids, cols);

    /*
     * clean up
     */
    deallocateMatrix(a);
    deallocateMatrix(b);
    deallocateMatrix(c);

    return t2 - t1;
}
#endif

#if defined(HAVE_CBLAS) || defined(HAVE_GSL_CBLAS)
/*---------------------------------------------------------------------------
 *
 * Compute matrix product using BLAS routine DGEMM.
 *
 * Input
 *   int argc        - length of argv[] array
 *   char* argv[]    - pointer to command line parameter array
 *   int verbosity   - program verification: verbosity > 0 gives more output
 *
 * Output
 *   double          - elapsed time for product computation
 */
double multiply_by_cblas(int argc, char* argv[], int verbosity)
{
    int rows, cols, mids;
    double **a, **b, **c;
    double t1, t2;
    double sec;
    double gflop_count;

    /*
     * process command line arguments
     */
    rows = atoi(argv[0]);
    mids = atoi(argv[1]);
    cols = atoi(argv[2]);
    gflop_count = 2.0 * rows * mids * cols / 1.0e9;

    if (verbosity > 0)
    {
        printf("CBLAS: rows = %d, mids = %d, columns = %d\n",
               rows, mids, cols);
    }

    /*
     * allocate and initialize matrices
     */
    a = (double**) allocateMatrix(rows, mids);
    b = (double**) allocateMatrix(mids, cols);
    c = (double**) allocateMatrix(rows, cols);
    initialize_matrices(a, b, c, rows, cols, mids, verbosity);

    /*
     * compute product
     */
    t1 = wtime();
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                rows, cols, mids, 1.0, &a[0][0], mids, &b[0][0], cols,
                0.0, &c[0][0], cols);
    t2 = wtime();
    sec = t2 - t1;

    if (verbosity > 1)
        printf("checksum = %f\n", checksum(c, rows, cols));

    printf("CBLAS:       %6.3f secs %6.3f GFlop/s (%5d x %5d x %5d)\n",
           sec, gflop_count / sec, rows, mids, cols);

    /*
     * clean up
     */
    deallocateMatrix(a);
    deallocateMatrix(b);
    deallocateMatrix(c);

    return t2 - t1;
}
#endif

/*---------------------------------------------------------------------------
 *
 * Compute matrix product using naive triple-nested loops.
 *
 * Input
 *   int argc        - length of argv[] array
 *   char* argv[]    - pointer to command line parameter array
 *   int verbosity   - program verification: verbosity > 0 gives more output
 *   char* order     - string indicating loop order, e.g., "ijk" or "jki"
 *
 * Output
 *   double          - elapsed time for product computation
 */
double multiply_by_naive_product(int argc, char* argv[], int verbosity,
                                 char* order)
{
    int rows, cols, mids;
    double **a, **b, **c;
    double t1, t2;
    double sec;
    double gflop_count;

    /*
     * process command line arguments
     */
    rows = atoi(argv[0]);
    mids = atoi(argv[1]);
    cols = atoi(argv[2]);
    gflop_count = 2.0 * rows * mids * cols / 1.0e9;

    if (verbosity > 0)
    {
        printf("naive(%3s): rows = %d, mids = %d, columns = %d\n",
               order, rows, mids, cols);
    }

    /*
     * allocate and initialize matrices
     */
    a = (double**) allocateMatrix(rows, mids);
    b = (double**) allocateMatrix(mids, cols);
    c = (double**) allocateMatrix(rows, cols);
    initialize_matrices(a, b, c, rows, cols, mids, verbosity);

    /*
     * compute product
     */
    t1 = wtime();
    do_product(a, b, c, 0, rows - 1, 0, cols - 1, 0, mids - 1);
    t2 = wtime();
    sec = t2 - t1;

    if (verbosity > 1)
        printf("checksum = %f\n", checksum(c, rows, cols));

    printf("naive(%s):  %6.3f secs %6.3f GFlop/s (%5d x %5d x %5d)\n",
           order, sec, gflop_count / sec, rows, mids, cols);

    /*
     * clean up
     */
    deallocateMatrix(a);
    deallocateMatrix(b);
    deallocateMatrix(c);

    return t2 - t1;
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix product using recursive tiling.
 *
 * Input
 *   int argc        - length of argv[] array
 *   char* argv[]    - pointer to command line parameter array
 *   int verbosity   - program verification: verbosity > 0 gives more output
 *   char* order     - string indicating loop order, e.g., "ijk" or "jki"
 *
 * Output
 *   double          - elapsed time for product computation
 */
double multiply_by_recursive_blocks(int argc, char* argv[], int verbosity,
                                    char* order)
{
    int rows, cols, mids, block_size;
    double **a, **b, **c;
    double t1, t2;
    double sec;
    double gflop_count;

    /*
     * process command line arguments
     */
    rows = atoi(argv[0]);
    mids = atoi(argv[1]);
    cols = atoi(argv[2]);
    block_size = atoi(argv[3]);
    gflop_count = 2.0 * rows * mids * cols / 1.0e9;

    if (verbosity > 0)
    {
        printf("Recursive blocks(%3s): rows = %d, mids = %d, columns = %d\n",
               order, rows, mids, cols);
        printf("block size = %d\n", block_size);
    }

    /*
     * allocate and initialize matrices
     */
    a = (double**) allocateMatrix(rows, mids);
    b = (double**) allocateMatrix(mids, cols);
    c = (double**) allocateMatrix(rows, cols);
    initialize_matrices(a, b, c, rows, cols, mids, verbosity);

    /*
     * compute product
     */
    t1 = wtime();
    mm_rec(c, a, b, 0, 0, 0, 0, 0, 0, rows, mids, cols, cols, block_size);
    t2 = wtime();
    sec = t2 - t1;

    if (verbosity > 1)
        printf("checksum = %f\n", checksum(c, rows, cols));

    printf("blocks(%3s): %6.3f secs %6.3f GFlop/s ",
           order, sec, gflop_count / sec);
    printf("(%5d x %5d x %5d) (%6d)\n", rows, mids, cols, block_size);

    /*
     * clean up
     */
    deallocateMatrix(a);
    deallocateMatrix(b);
    deallocateMatrix(c);

    return t2 - t1;
}

/*---------------------------------------------------------------------------
 *
 * Compute matrix product using tiling.  The loop order used for the tile
 * products is specified in string variable "mode".
 *
 * Input
 *   int argc        - length of argv[] array
 *   char* argv[]    - pointer to command line parameter array
 *   int verbosity   - program verification: verbosity > 0 gives more output
 *   char* order     - string indicating loop order, e.g., "ijk" or "jki"
 *
 * Output
 *   double          - elapsed time for product computation
 */
double multiply_by_tiles(int argc, char* argv[], int verbosity, char* order)
{
    int rows, cols, mids;
    int rows_per_tile, cols_per_tile, mids_per_tile;
    int row_start, row_end;
    int col_start, col_end;
    int mid_start, mid_end;
    double **a, **b, **c;
    double t1, t2;
    double sec;
    double gflop_count;

    /*
     * process command line arguments
     */
    rows = atoi(argv[0]);
    mids = atoi(argv[1]);
    cols = atoi(argv[2]);
    rows_per_tile = atoi(argv[3]);
    mids_per_tile = atoi(argv[4]);
    cols_per_tile = atoi(argv[5]);
    gflop_count = 2.0 * rows * mids * cols / 1.0e9;

    if (verbosity > 0)
    {
        printf("Tiles(%3s): rows = %d, mids = %d, columns = %d\n",
               order, rows, mids, cols);
        printf("block rows = %d, mids = %d, columns = %d\n",
               rows_per_tile, mids_per_tile, cols_per_tile);
    }

    /*
     * allocate and initialize matrices
     */
    a = (double**) allocateMatrix(rows, mids);
    b = (double**) allocateMatrix(mids, cols);
    c = (double**) allocateMatrix(rows, cols);
    initialize_matrices(a, b, c, rows, cols, mids, verbosity);

    /*
     * compute product
     */
    t1 = wtime();
    for (row_start = 0; row_start < rows; row_start += rows_per_tile)
    {
        row_end = row_start + rows_per_tile - 1;
        if (row_end >= rows) row_end = rows - 1;
        for (col_start = 0; col_start < cols; col_start += cols_per_tile)
        {
            col_end = col_start + cols_per_tile - 1;
            if (col_end >= cols) col_end = cols - 1;
            for (mid_start = 0; mid_start < mids; mid_start += mids_per_tile)
            {
                mid_end = mid_start + mids_per_tile - 1;
                if (mid_end >= mids) mid_end = mids - 1;
                do_product(a, b, c, row_start, row_end, col_start,
                           col_end, mid_start, mid_end);
            }
        }
    }
    t2 = wtime();
    sec = t2 - t1;

    if (verbosity > 1)
        printf("checksum = %f\n", checksum(c, rows, cols));

    printf("tiles(%3s):  %6.3f secs %6.3f GFlop/s ",
            order, sec, gflop_count / sec);
    printf("(%5d x %5d x %5d) (%4d x %4d x %4d)\n",
            rows, mids, cols, rows_per_tile, mids_per_tile,
            cols_per_tile);

    /*
     * clean up
     */
    deallocateMatrix(a);
    deallocateMatrix(b);
    deallocateMatrix(c);

    return t2 - t1;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

int main(int argc, char* argv[])
{
    int c;
    int mode = MODE_UNKNOWN;
    int verbosity = 0;
    char* pname = argv[0];
    char order[4];

    /*
     * process command line to determine solution method.  Begin by
     * constructing option string for getopt() depending on available
     * libraries.
     */
    const char* optstring =
#if defined(HAVE_BLAS)
        "b"
#endif
#if defined(HAVE_CBLAS) || defined(HAVE_GSL_CBLAS)
        "c"
#endif
        "bn:r:t:v";

    while ((c = getopt(argc, argv, optstring)) != -1)
    {
        switch (c)
        {
            case 'b':
                mode = MODE_BLAS;
                break;
            case 'c':
                mode = MODE_CBLAS;
                break;
            case 'n':
                mode = MODE_NAIVE;
                strncpy(order, optarg, 4);
                order[3] = 0; /* just to make sure... */
                break;
            case 'r':
                mode = MODE_RECURSIVE_TILE;
                strncpy(order, optarg, 4);
                order[3] = 0; /* just to make sure... */
                break;
            case 't':
                mode = MODE_STD_TILE;
                strncpy(order, optarg, 4);
                order[3] = 0; /* just to make sure... */
                break;
            case 'v':
                verbosity++;
                break;
            default:
                break;
        }
    }
    argc -= optind;
    argv += optind;

    if (mode == MODE_UNKNOWN)
    {
        fprintf(stderr, "Usage: one of:\n\n");
        fprintf(stderr, "  %s -n ijk M P N              (naive)\n",
                pname);
        fprintf(stderr, "  %s -t ijk M P N TM TP TN     (standard tiling)\n",
                pname);
        fprintf(stderr, "  %s -r ijk M P N TS           (recursive tiling)\n",
                pname);
#if defined(HAVE_BLAS)
        fprintf(stderr, "  %s -b M P N                  (BLAS)\n",
                pname);
#endif
#if defined(HAVE_CBLAS) || defined(HAVE_GSL_CBLAS)
        fprintf(stderr, "  %s -c M P N                  (CBLAS)\n",
                pname);
#endif
        fprintf(stderr, "\nThe product is MxN, factors are MxP and PxN, ");
        fprintf(stderr, "and 'ijk' may be replace by\n");
        fprintf(stderr, "any permution of 'i', 'j', and 'k'.\n\n");
        fprintf(stderr, "TM, TP, and TN are corresponding tile dimensions ");
        fprintf(stderr, "and TS is the maximum\n");
        fprintf(stderr, "number of elements in a tile.\n");
        exit(EXIT_FAILURE);
    }

    /*
     * unless we're using the BLAS we need to set do_function() to be
     * the function that uses the desired ordering.
     */
    if (mode != MODE_BLAS && mode != MODE_CBLAS)
    {
        if (strncmp(order, "ijk", 3) == 0) do_product = do_ijk_product;
        if (strncmp(order, "ikj", 3) == 0) do_product = do_ikj_product;
        if (strncmp(order, "jik", 3) == 0) do_product = do_jik_product;
        if (strncmp(order, "jki", 3) == 0) do_product = do_jki_product;
        if (strncmp(order, "kij", 3) == 0) do_product = do_kij_product;
        if (strncmp(order, "kji", 3) == 0) do_product = do_kji_product;
        if (do_product == NULL)
        {
            fprintf(stderr, "loop order string not recognized: ");
            fprintf(stderr, "must be one of: ijk, ikj, jik, jki, kij, kji\n");
            exit(EXIT_FAILURE);
        }
    }

    /*
     * compute the product!
     */
    switch (mode)
    {
#if defined(HAVE_BLAS)
        case MODE_BLAS:
            if (argc < 3)
            {
                fprintf(stderr, "usage: %s -b ROWS MIDS COLS\n", pname);
                exit(EXIT_FAILURE);
            }
            multiply_by_blas(argc, argv, verbosity);
            break;
#endif
#if defined(HAVE_CBLAS) || defined(HAVE_GSL_CBLAS)
        case MODE_CBLAS:
            if (argc < 3)
            {
                fprintf(stderr, "usage: %s -c ROWS MIDS COLS\n", pname);
                exit(EXIT_FAILURE);
            }
            multiply_by_cblas(argc, argv, verbosity);
            break;
#endif
        case MODE_NAIVE:
            if (argc < 3)
            {
                fprintf(stderr, "usage: %s -n ijk ROWS MIDS COLS\n", pname);
                exit(EXIT_FAILURE);
            }
            multiply_by_naive_product(argc, argv, verbosity, order);
            break;
        case MODE_STD_TILE:
            if (argc < 6)
            {
                fprintf(stderr, "usage: %s -t ijk ROWS MIDS COLS", pname);
                fprintf(stderr, " TILE_ROWS TILE_MIDS TILE_COLS\n");
                exit(EXIT_FAILURE);
            }
            multiply_by_tiles(argc, argv, verbosity, order);
            break;
        case MODE_RECURSIVE_TILE:
            if (argc < 4)
            {
                fprintf(stderr, "usage: %s -r ROWS MIDS COLS TILE_SIZE\n",
                        pname);
                exit(EXIT_FAILURE);
            }
            multiply_by_recursive_blocks(argc, argv, verbosity, order);
            break;
    }

    return 0;
}

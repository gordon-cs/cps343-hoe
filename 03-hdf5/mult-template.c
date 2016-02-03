/*
 * $Smake: gcc -Wall -O3 -o %F %f
 *
 * Matrix-matrix product
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

/* Macro to index matrices in column-major (Fortran) order */
#define IDX(i,j,stride) ((i)+(j)*(stride))  /* column major */

/*----------------------------------------------------------------------------
 * Usage
 */
void usage( char* program_name )
{
    fprintf( stderr, "Usage: %s [-v]\n", program_name );
}

/*----------------------------------------------------------------------------
 * Dump Matrix
 */
void dumpMatrix( double* a, int m, int n, int stride )
{
    int i, j;
    for ( i = 0; i < m; i++ )
    {
        for ( j = 0; j < n; j++ )
        {
            printf( " %8.2f", a[IDX(i,j,stride)] );
        }
        printf( "\n" );
    }
    printf( "\n" );
}

/*----------------------------------------------------------------------------
 * create Matrix based on supplied name
 */
void createMatrix( char* name, double** a, int* rows, int* cols )
{
    if ( strcmp( name, "A" ) == 0 )
    {
        *rows = 4;
        *cols = 2;
        *a = (double*) malloc( *rows * *cols * sizeof( double ) );
        (*a)[IDX(0,0,*rows)] =  4.0;
        (*a)[IDX(1,0,*rows)] =  2.0;
        (*a)[IDX(2,0,*rows)] = -2.0;
        (*a)[IDX(3,0,*rows)] =  1.0;
        (*a)[IDX(0,1,*rows)] = -4.0;
        (*a)[IDX(1,1,*rows)] = -1.0;
        (*a)[IDX(2,1,*rows)] = -3.0;
        (*a)[IDX(3,1,*rows)] =  4.0;
    }
    else if ( strcmp( name, "B" ) == 0 )
    {
        *rows = 2;
        *cols = 3;
        *a = (double*) malloc( *rows * *cols * sizeof( double ) );
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
 */
void matmat_jki( double* c, double* a, int nrow_a, int ncol_a,
                 double* b, int ncol_b )
{
    const int nrow_b = ncol_a;
    const int nrow_c = nrow_a;
    int i, j, k;
    for ( j = 0; j < ncol_b; j++ )
    {
        for ( i = 0; i < nrow_a; i++ ) c[IDX(i,j,nrow_c)] = 0.0;
        for ( k = 0; k < ncol_a; k++ )
        {
            for ( i = 0; i < nrow_a; i++ )
            {
                c[IDX(i,j,nrow_c)] += a[IDX(i,k,nrow_a)] * b[IDX(k,j,nrow_b)];
            }
        }
    }
}

/*----------------------------------------------------------------------------
 * Main program
 */
int main( int argc, char* argv[] )
{
    double* a;
    double* b;
    double* c;
    int nrow_a, ncol_a;
    int nrow_b, ncol_b;
    int nrow_c, ncol_c;
    int verbose = 0;

    /* Process command line */
    int ch;
    while ( ( ch = getopt( argc, argv, "v" ) ) != -1 )
    {
        switch ( ch )
        {
            case 'v':
                verbose++;
                break;
            default:
                usage( argv[0] );
                return EXIT_FAILURE;
        }
    }
    argv[optind - 1] = argv[0];
    argv += ( optind - 1 );
    argc -= ( optind - 1 );

    if ( argc != 1 )
    {
        usage( argv[0] );
        return EXIT_FAILURE;
    }

    /* create matrix data and optionally display it */
    createMatrix( "A", &a, &nrow_a, &ncol_a );
    createMatrix( "B", &b, &nrow_b, &ncol_b );

    if ( ncol_a != nrow_b )
    {
        fprintf( stderr, "Error: matrix dimensions are not compatible\n" );
        return EXIT_FAILURE;
    }

    if ( verbose )
    {
        printf( "Matrix A:\n" );
        dumpMatrix( a, nrow_a, ncol_a, nrow_a );
        printf( "Matrix B:\n" );
        dumpMatrix( b, nrow_b, ncol_b, nrow_b );
    }

    /* Compute matrix product C = AB and optionally display it */
    nrow_c = nrow_a;
    ncol_c = ncol_b;
    c = (double*) malloc( nrow_c * ncol_c * sizeof( double ) );

    matmat_jki( c, a, nrow_a, ncol_a, b, ncol_b );

    printf( "Matrix C = AB:\n" );
    dumpMatrix( c, nrow_c, ncol_c, nrow_c );

    /* Clean up and quit */
    free( a );
    free( b );
    free( c );
    return 0;
}

/*
 * $Smake: mpic++ -O3 -o %F %f
 *
 * Copyright (c) 2010 Jonathan Senning <jonathan.senning@gordon.edu.
 * Department of Mathematics and Computer Science
 * Gordon College, Wenham MA 01984-1899
 * Revised March 2018 - display communication and computation timing data
 * Revised March 2024 - Added command-line switch for timing information
 *
 * Implements Parallel Sorting with Regular Sampling (PSRS) using the
 * C library qsort() function for local sorting.
 */

#include <cstdio>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>

//----------------------------------------------------------------------------
//
// Computes the starting and ending displacements for the ith
// subinterval in an n-element array given that there are m
// subintervals of approximately equal size.
//
// Input:
//    int n    - length of array (array indexed [0]..[n-1])
//    int m    - number of subintervals
//    int i    - subinterval number
//
// Output:
//    int* s   - location to store subinterval starting index
//    int* e   - location to store subinterval ending index
//
// Suppose we want to partition a 100-element array into 3
// subintervals of roughly the same size.  The following three
// pairs of calls find the starting and ending indices of each
// subinterval:
//   decompose1d(100, 3, 0, &s, &e);  (now s =  0, e = 33)
//   decompose1d(100, 3, 1, &s, &e);  (now s = 34, e = 66)
//   decompose1d(100, 3, 2, &s, &e);  (now s = 67, e = 99)
//
// The subinterval length can be computed with e - s + 1.
//
// Based on the FORTRAN subroutine MPE_DECOMP1D in the file
// UsingMPI/intermediate/decomp.f supplied with the book
// "Using MPI" by Gropp et al.  It has been adapted to use
// 0-based indexing.

static void decompose1d(int n, int m, int i, int* s, int* e)
{
    const int length  = n / m;
    const int deficit = n % m;
    *s =  i * length + (i < deficit ? i : deficit);
    *e = *s + length - (i < deficit ? 0 : 1);
    if ((*e >= n) || (i == m - 1)) *e = n - 1;
}

//-----------------------------------------------------------------------------
// Compare two numeric values and return integer indicating their ordinal
// ranking
//
// Input:
//    const void* a, b   - pointers to elements of comparible numerical type
// Output:
//    none
// Returns:
//    positive int if *a > *b, negative int if *a < *b, and 0 if *a == *b

int cmp(const void* a, const void* b)
{
    return *((int*) a) - *((int*) b);
}

//-----------------------------------------------------------------------------
// Performs parallel sort with regular sampling
//
// Input:
//    int* inList        - pointer to first element of initial list
//    int  inLen         - length of initial list
//    int  master        - rank of MPI master process (determines splitters)
//    MPI_Comm comm      - MPI communicator common to all sorting processes
//    bool doTiming      - true for timing information
// Output:
//    int* outList       - pointer to first element of output (sorted) list
//    int  outLen        - length of output (sorted) list
// Returns:
//    0  if success; -1 to indicate error

int psrsSort(int* inList, int inLen, int* outList, int* outLen, int master,
             MPI_Comm comm, bool doTiming)
{
    enum {SPLITTER, BIN};
    int rank;              // process id
    int size;              // number of processes
    int* sample  = NULL;   // list of samples from local sorted list
    int* samples = NULL;   // samples from all processes (rank 0 only)
    int* splitter = NULL;  // list of splitters
    int* binSize  = NULL;  // list of partial destination bin sizes
    int* binSizes = NULL;  // list of all partial destination bin sizes
    int  n, m;             // temporary variables array indices

    // get communicator parameters
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    if (doTiming)
    {
        if (rank == 0)
        {
            printf("\n");
            printf("                broad-     all-            ");
            printf("    pre-   split-    post-      comm\n");
            printf("rank   gather     cast   gather sendrecv   ");
            printf("    sort     list     sort     fract\n");
            printf("---- -------- -------- -------- --------   ");
            printf("-------- -------- --------   -------\n");
        }
        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    // sort original list
    double t0 = MPI_Wtime();
    qsort(inList, inLen, sizeof(int), cmp);
    double time_sort1 = MPI_Wtime() - t0;

    // collect set of regularly-spaced samples from sorted list and
    // gather the samples from all processes on master process in array
    sample = new int [size];
    samples = (rank == master ? new int [size * size] : NULL);
    for (int i = 0; i < size; i++) sample[i] = inList[(i * inLen) / size];
    t0 = MPI_Wtime();
    MPI_Gather(sample, size, MPI_INT, samples, size, MPI_INT, master, comm);
    double time_gather = MPI_Wtime() - t0;

    // master process sorts list of samples to determine splitter values
    // and broadcasts lists of splitters to all processes
    splitter = new int [size - 1];
    if (rank == master)
    {
        qsort(samples, size * size, sizeof(int), cmp);
        for (int i = 1; i < size; i++) splitter[i - 1] = samples[size * i];
    }
    t0 = MPI_Wtime();
    MPI_Bcast(splitter, size - 1, MPI_INT, master, comm);
    double time_broadcast = MPI_Wtime() - t0;
    delete [] samples;
    delete [] sample;

    // splitters will be used to partition the inList; here we compute the
    // length of each i_th portion of the list and store it in binSize[i]
    binSize = new int [size];
    t0 = MPI_Wtime();
    n = m = 0;
    for (int i = 0; i < size - 1; i++)
    {
        while (m < inLen && inList[m] < splitter[i]) m++;
        binSize[i] = m - n;
        n = m;
    }
    binSize[size - 1] = inLen - n;
    double time_split = MPI_Wtime() - t0;
    delete [] splitter;    

    // distribute all bin size information to all processes and
    // determine length of list needed to hold all data this process
    // will receive; binSizes[i*size+j] is the number of elements
    // that process i will send to process j.
    binSizes = new int [size * size];
    t0 = MPI_Wtime();
    MPI_Allgather(binSize, size, MPI_INT, binSizes, size, MPI_INT, comm);
    double time_allgather = MPI_Wtime() - t0;
    n = 0;
    for (int i = 0; i < size; i++) n += binSizes[i * size + rank];
    if (*outLen < n) return -1;
    *outLen = n;

    // exchange partial bins with other processes.
    n = m = 0;
    t0 = MPI_Wtime();
    for (int i = 0; i < size; i++)
    {
        int j = i * size + rank;
        MPI_Sendrecv(&inList[n],  binSize[i],  MPI_INT, i, BIN,
                     &outList[m], binSizes[j], MPI_INT, i, BIN,
                     comm, MPI_STATUS_IGNORE);
        n += binSize[i];
        m += binSizes[j];
    }
    double time_sendrecv = MPI_Wtime() - t0;
    delete [] binSizes;

    // finally, sort received data
    t0 = MPI_Wtime();
    qsort(outList, *outLen, sizeof(int), cmp);
    double time_sort2 = MPI_Wtime() - t0;

    // compute elapsed communication and sorting time
    double tcomm = time_gather+time_broadcast+time_allgather+time_sendrecv;
    double tcomp = time_sort1 + time_sort2 + time_split;
    double ratio = tcomm / (tcomm + tcomp);

    // display detailed timing information if requested
    if (doTiming)
    {
        printf("%4d %8.5f %8.5f %8.5f %8.5f | %8.5f %8.5f %8.5f | %7.4f\n",
               rank,
               time_gather, time_broadcast, time_allgather, time_sendrecv,
               time_sort1, time_split, time_sort2, ratio);
    }

    return 0;
}

//-----------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    bool doTiming = false; // show detailed timing information
    int rank;              // process id
    int size;              // number of processes
    int nelem;             // number of elements in global list
    int n = 0;             // number of elements in local list
    int outLen;            // number of elements in local sorted list
    int start;             // starting value for range
    int end;               // ending value for range
    int* list = NULL;      // pointer to local portion of unsorted list
    int* outList = NULL;   // pointer to local portion of sorted list;
    double t0, t1;         // timing variables
    int status;            // sorting status (0 means okay)

    int opt;
    while ((opt = getopt(argc, argv, "t")) != -1)
    {
        switch (opt)
        {
            case 't':
                doTiming = true;
                break;
            default:
                fprintf(stderr, "usage: %s [-t] N\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    argv[optind - 1] = argv[0];
    argc -= (optind - 1);
    argv += (optind - 1);

    // initalize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // initalize random number generator
    srandom((unsigned int) (time(NULL) * rank));

    // get total list length from command line
    if (argc != 2)
    {
        fprintf(stderr, "usage: %s [-t] N\n", argv[0]);
        exit(EXIT_FAILURE);
    }
    nelem = atol(argv[1]);
    if (nelem <= 0)
    {
        fprintf(stderr, "Error: list length N must be positive\n");
        exit(EXIT_FAILURE);
    }

    // determine the length of the local portion of the list
    decompose1d(nelem, size, rank, &start, &end);
    n = end - start + 1;

    // generate list of random numbers
    list = new int [n];
    for (int i = 0; i < n; i++) list[i] = random();

    // allocate memory for sorted list
    outLen = 2 * n;
    outList = new int [outLen];

    // sort the list
    t0 = MPI_Wtime();
    status = psrsSort(list, n, outList, &outLen, 0, MPI_COMM_WORLD,
                      doTiming);
    t1 = MPI_Wtime();

    if (status < 0) fprintf(stderr, "Error in psrsSort()\n");

    // display sorted list if it is fairly short
    if (nelem <= 100)
    {
        for (int r = 0; r < size; r++)
        {
            if (rank == r)
            {
                for (int i = 0; i < outLen; i++)
                {
                    printf("rank %2d: list[%3d] = %12d\n",
                           rank, i, outList[i]);
                }
                fflush(stdout);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }

    MPI_Finalize();

    // report sort time
    if (rank == 0)
    {
        printf("\n%12d element list sorted in %12.8f seconds\n", 
               nelem, t1 - t0);
    }

    // all done; clean up
    delete [] list;
    delete [] outList;

    return 0;
}

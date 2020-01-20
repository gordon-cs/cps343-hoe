Compiling the programs
======================

All the programs in each subdirectory directory can be compiled by
typing `make` at the command prompt.

Note that in order to compile `pi_mpi.cc` the MPI development headers
and libraries must be available.  In a system using Lmod environment
modules this can be accomplished by typing a command like `module load
openmpi` or `module load mpich`.

Description and instructions
============================

Subdirectory `compute_pi`
-----------------------

* `pi_serial.cc`

  This program computes an approximation of pi using a serial
  algorithm.  The actual computation is done by approximating the
  integral of 4/(1+x^2) on the interval [0,1] with the midpoint rule
  with 400,000,000 subintervals.

* `pi_omp.cc`

  Same as above except OpenMP is used to provide multithreading if
  supported by the hardware.  In particular, a parallel thread is
  created for each processor or processor core on the machine.  This
  can be overridden by setting the environment variable
  `OMP_NUM_THREADS`.  For example `OMP_NUM_THREADS=1 ./pi_omp` will
  force the program to use a single thread regardless of of the
  hardware.

* `pi_mpi.cc`

  This version uses MPI and is suitable for distributed memory
  system like a cluster.  To run on a single multicore machine, use
  ```shell
  mpiexec -n <N> ./pi_mpi
  ```
  with <N> replaced by the number of processes to be used.  To run on a
  cluster try
  ```shell
  mpiexec -machinefile <file> -n <N> ./pi_mpi
  ```

  where `<file>` is a file containing the names of the individual machines
  that comprise the cluster.  If the cluster has the resource manager SLURM
  then try
  ```shell
  salloc -Q -n <N> mpiexec ./pi_mpi
  ```

Subdirectory `matrix_matrix_multiply`
-------------------------------------


* `matmat_c1.c` - Matrices are statically allocated 2-Dimensional arrays

  Since the value of N is known at compile-time, the matrices are
  declared with statements like `double a[N][N]`.  They are global so
  that they will not be allocated on the stack.  The (i,j) array
  element is referred to with `a[i][j]`.

* `matmat_c2.c` - Matrices are dynamically allocated 1-Dimensional arrays
  (row-major)

  The C library routine `malloc()` is used to allocate each matrix as
  an (N^2)-element one-dimensional array.  A macro `IDX()` is defined
  to facilitate indexing into the array using row and column values and
  here assumes a *row-major* ordering, i.e., elements along a row in
  the matrix are stored consecutively in memory.  The (i,j) array element
  is referred to with `a[IDX(i,j,stride)]` where `stride` is the matrix
  row length.

* `matmat_c3.c` - Matrices are dynamically allocated 1-Dimensional arrays
  (column-major)

  Same as `matmat_c2.c` except *column-major* ordering is used.  The only
  change to the code is the definition of the `IDX()` macro.  This also
  means that `stride` is the matrix column length.

* `matmat_c4.c` - Matrices are dynamically allocated using an array of pointers

  To allow dynamic memory allocation but to facilitate normal C-style
  two-dimensional array indexing even when matrices are passed as
  function parameters, this program allocates an (N^2)-element
  one-dimensional array to hold the matrix data but additionally
  allocates an N-element one-dimensional array of pointers.  Each
  element of this latter array is set to point to the first element of
  a matrix row.  The (i,j) array element is referred to with `a[i][j]`.

* `matmat_f77.f` - Written in FORTRAN 77

* `matmat_ijk.c` - Template with *ijk* and *jki* loop-orderings; user writes
  for *ikj*, *jik*, *kij*, and *kji* orderings.

* `matmat_ijk_macro.c` - same as `matmat_ijk.c` except this version uses
  1D arrays for the matrices and a macro for indexing.

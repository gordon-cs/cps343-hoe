Compiling the programs
======================

All the programs in this directory can be compiled by typing `make` at
the command prompt.

Description of files
====================

Matrix-Matrix multiplication
----------------------------

* matmat_c1.c - Matrices are statically allocated 2-Dimensional arrays

  Since the value of N is known at compile-time, the matrices are
  declared with statements like `double a[N][N]`.  They are global so
  that they will not be allocated on the stack.  The (i,j) array
  element is referred to with `a[i][j]`.
  
* matmat_c2.c - Matrices are dynamically allocated 1-Dimensional arrays

  The C library routine `malloc()` is used to allocate each matrix as
  an (N^2)-element one-dimensional array.  A macro `IDX()` is defined
  to facilitate indexing into the array using row and column values.
  The (i,j) array element is referred to with `a[IDX(i,j)]`.
  
* matmat_c3.c - Matrices are dynamically allocated using an array of pointers

  To allow dynamic memory allocation but to facilitate normal C-style
  two-dimensional array indexing even when matrices are passed as
  function parameters, this program allocates an (N^2)-element
  one-dimensional array to hold the matrix data but additionally
  allocates an N-element one-dimensional array of pointers.  Each
  element of this latter array is set to point to the first element of
  a matrix row.  The (i,j) array element is referred to with `a[i][j]`.

* matmat_f77.f - Uses FORTRAN 77

* matmat_f95.f - Uses FORTRAN 95

* matmat_ijk.c - Template with ijk and jki loop-orderings; user writes for ikj, jik, kij, and kji orderings.

* timing.c - Timing test for matrix multiplication with ijk loop ordering using two different matrix storage schemes.
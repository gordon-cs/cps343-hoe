Compiling the programs
======================

All the programs in this directory can be compiled by typing 
```shell
make
```
at the command prompt.

Description and instructions
============================

pi_serial.cc	This program computes an approximation of pi using a serial
		algorithm.  The actual computation is done by approximating
		the integral of 4/(1+x^2) on the interval [0,1] with the
		midpoint rule with 400,000,000 subintervals.

pi_omp.cc	Same as above except OpenMP is used to provide multithreading
		if supported by the hardware.  In particular, a parallel
		thread is created for each processor or processor core on the
		machine.  This can be overridden by setting the environment
		variable OMP_NUM_THREADS.  For example
			OMP_NUM_THREADS=1 ./pi_omp
		will force the program to use a single thread regardless of
		of the hardware.

pi_omp_dyn.cc	This version also uses OpenMP but with dynamic load balancing.
		Each thread computes a small task and then asks for another
		one.  Since the tasks in this program are of uniform size this
		will not provide any improvement, but is very useful if tasks
		require differing amounts of time to carry out.

pi_mpi.cc	This version uses MPI and is suitable for distributed memory
		system like a cluster.  To run on a single machine, use
			mpirun -n <N> pi_mpi
		with <N> replaced by the number of processes to be used.  To
		run on a cluster try
			mpirun -machinefile <file> -n <N> pi_mpi
		where <file> is a file containing the names of the individual
		machines that comprise the cluster.  If you're using SLURM
		then the command would be
			salloc -n<N> mpirun pi_mpi

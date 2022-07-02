# Gauss-Seidel Method

<!-- Compile instructions -->
make sequential
make mpi
make omp

## Or you can use:
1. For sequential: 
	gcc -std=c99  gseidel_sequential.c -o sequential.out -lm

2. For OMP:
	gcc  gseidel_mp_schedule.c -std=c99 -lm -fopenmp -o openmp.out

3. For MPI:
	mpicc gseidel_MPI_by_cols.c -std=c99 -lm -o openmpi.out 

<!-- Windows instructions -->
    Remove -lm attribute

<!-- OSX instructions -->
    To avoid message: "A system call failed during shared memory initialization that should
    not have.  It is likely that your MPI job will now either abort or
    experience performance degradation."

    Please type "export TMPDIR=/tmp" (without quotes) at your terminal

    At compile, Remove -lm attribute.

<!-- Linux instructions -->
    You must compile using -lm

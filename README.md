# Gauss-Seidel Method

## Compile instructions

Use:

1. make sequential

or

2. make mpi

or

3. make omp

Or you can use:

4. For sequential:

	gcc -std=c99  gseidel_sequential.c -o sequential.out -lm

5. For OMP:

	gcc  gseidel_mp_schedule.c -std=c99 -lm -fopenmp -o openmp.out

6. For MPI:

	mpicc gseidel_MPI_by_cols.c -std=c99 -lm -o openmpi.out 

## Windows instructions

    To compile, Remove -lm attribute


## OSX instructions

    To avoid message: "A system call failed during shared memory initialization that should
    not have.  It is likely that your MPI job will now either abort or
    experience performance degradation."

    Please type "export TMPDIR=/tmp" (without quotes) at your terminal

    To compile, Remove -lm attribute.

## Linux instructions
    You must compile using -lm
    
# Running at all operation systems

For sequential: ./sequential.out [precision] [filename]

For OMP: ./omp.out [precision] [filename] [num_threads]

For MPI: mpirun -np [num_process] ./mpi.out [precision] [filename]

# DISCLAIMER:
until MPI 0.1.0 version you only calculate matrix multiply processes numbers.

e.g. 4 process, matrix 16, 20, 24, ...

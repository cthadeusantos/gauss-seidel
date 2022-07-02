sequential: 
	clear
	gcc -std=c99  gseidel_sequential.c -o sequential.out -lm

omp:
	clear
	gcc  gseidel_mp_schedule.c -std=c99 -lm -fopenmp -o openmp.out

mpi:
	clear
	mpicc gseidel_MPI_by_cols.c -std=c99 -lm -o openmpi.out 

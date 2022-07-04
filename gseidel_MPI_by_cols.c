/* Program Gauss-Seidel
   Solution of a system of linear equations by Gauss-Seidel's
   iteration method. Assume that the coefficient matrix satisfies
   the condition of convergence.

    A.x = B

    Variables:
    nvar            dimension of the matrix
    tolerance       precision error
    matrix          matrix A (coefficient matrix)
    vector          matrix B (constant matrix)
    solution        matrix x (variable matrix)
    flag            tolerance reached

    OSX INSTRUCTIONS
    To avoid message: "A system call failed during shared memory initialization that should
    not have.  It is likely that your MPI job will now either abort or
    experience performance degradation."

    Please type "export TMPDIR=/tmp" (without quotes) at your terminal

    LINUX INSTRUCTIONS
    You must compile using -lm

    References:
    https://www.rookiehpc.com/mpi/docs/mpi_scatterv.php
    https://stackoverflow.com/questions/9269399/sending-blocks-of-2d-array-in-c-using-mpi/9271753#9271753
   */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdbool.h>
#include<string.h>
#include<time.h>
#include<mpi.h>
#include "functools.h"

// #define max(a, b) ((a) > (b) ? (a) : (b)) // replace by get_max at functools.h

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    MPI_Datatype subarray_type;

    double start, end; // To compute running time

    int rank, size;
    MPI_Comm_size (MPI_COMM_WORLD, &size);  // running rank number
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);  // Processes quantity

    if (argc < 3){
        printf("Usage: output [num_decimals] [filename] [0-1]\n");
        MPI_Finalize();
        exit(1);
    }

    clock_t t; //variável para armazenar tempo

	long double tolerance = 1.000, sum=0, temp=0, global_delta=0;
    long double local_sum, delta_x, solucao;

    char *filename;
    FILE *infile;
	
    int offset, flag=0, i, j, k, n, aux, index_mtx=0;

    int root=0, i_local=0, nlocal=0;

    int dimension = 0, ncols_submatrix=0, nelements_submatrix=0;
    
    /*
     * DEFINE DIMENSIONS FROM FIRST LINE AT FILE
     */

    if (rank == 0){
        // Read command-line arguments
        i = (-1) * atoi(argv[1]);       
        filename = argv[2];
        if (argv[3] != NULL)
            if (atoi(argv[3]) == 1)
                flag = atoi(argv[3]);
        tolerance = tolerance * pow(10, i);
        dimension = read_dimension(filename);
        
        // Set submatrix' size
        ncols_submatrix = (int)dimension / size;
        // Adjusting if round it down
        if (ncols_submatrix * size < dimension) ncols_submatrix++;
        nelements_submatrix = ncols_submatrix * dimension;
    }

    MPI_Bcast( &dimension , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &tolerance , 1 , MPI_LONG_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // Define parameters to create submatrix and 
    nelements_submatrix = get_size_submatrix(rank, size, ncols_submatrix, dimension);
    ncols_submatrix = get_num_lines_submtx(rank, size, dimension);

    /*
    * Create submatrix and vectors for all processes
    */
    long double *vector = calloc(ncols_submatrix , sizeof(long double));
    long double *submatrix = calloc(nelements_submatrix, sizeof(long double));
    long double *solution = calloc(ncols_submatrix , sizeof(long double));

    // Init solution vector 
    for(i=0;i<ncols_submatrix;i++)
        solution[i]=0;

    // READ MATRIX AND VECTOR PROCEDURE
    if (rank == 0){
        printf("2.Montando matrizes, rank %i DIMENSAO %i!\n", rank, dimension);

        long double *matrix = calloc(dimension * dimension , sizeof(long double));
        long double *vector_aux = calloc(dimension, sizeof(long double));

        // READ MATRIX & VECTOR FILE
        read_matrix(filename, matrix, vector_aux);

        // SPLIT MATRIX AND SEND TO PROCESSES
        int linha, coluna;
        for (i = 0;i < size; i++){ // SPLIT MATRIX IN COLUMNS
            offset = matrix_offset(i, size, dimension);
            long double *auxiliary = calloc(offset * dimension, sizeof(long double));
            int count = 0;
            int lower_limit, upper_limit;

            for (linha=0; linha < dimension; linha++){  // Process all lines
                lower_limit=offset * i;
                upper_limit=offset * ( i + 1 );

                // Process columns from submatrix limits right and left 
                for (coluna=lower_limit; coluna < upper_limit; coluna++){
                    int index = linha * dimension + coluna;
                    if (i!=0)
                        auxiliary[count] = matrix[index];
                    else
                        submatrix[count] = matrix[index];
                    count++;
                }
            }
            if (i != 0){ // SEND SUBMATRIX AND VECTOR(PARCIALLY) TO RANKS
                MPI_Isend(auxiliary, dimension * offset, MPI_LONG_DOUBLE, i, 4000+i, MPI_COMM_WORLD, &request);
                int index = offset * i;
                int nelements = get_num_lines_submtx(i, size, dimension);
                MPI_Send(&vector_aux[index], nelements, MPI_LONG_DOUBLE, i, 3000+i, MPI_COMM_WORLD);
            } else {    // COPY VECTOR_AUX TO VECTOR (RANK 0)
                int nelements = get_num_lines_submtx(i, size, dimension);
                for (linha = 0; linha < nelements; linha++){
                    vector[linha] = vector_aux[linha];
                }
            }
            free(auxiliary);
        }
        free(vector_aux);
        free(matrix);
    } 

    if (rank != 0){     // RECEIVE SUBMATRIX & VECTOR
        MPI_Recv(submatrix, nelements_submatrix, MPI_LONG_DOUBLE, 0, 4000+rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int nelements = get_num_lines_submtx(i, size, dimension);
        MPI_Recv(&vector[0], nelements, MPI_LONG_DOUBLE, 0, 3000+rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD); // Sync before start operation

    if (rank == 0){
        printf("****Waiting! MPI calculating! \n" );
    }

    nlocal = get_num_lines_submtx(rank, size, dimension);
    solucao=0.0;
    delta_x = 0.0;
    start = MPI_Wtime();
    // int iterations = 0;
    // int global_iter = 0;
    
    long double global_delta_previous = 0;  
    delta_x = 0.0;
    
    do {
        delta_x = 0.0;
        global_delta_previous = global_delta;
        for (i = 0; i < dimension; i++){    // Process all lines
            local_sum = 0.0;
            
            for (j = 0; j < nlocal; j++){   // Process submatrix columns    
                if (j + rank * nlocal != i){    // NON DIAGONAL ELEMENT
                    index_mtx = i * nlocal + j; // Indice da matriz
                    local_sum = local_sum + submatrix[index_mtx] * solution[j];
                }
            }

            root = i / nlocal;      // Define rank que a submatrix pertence
            i_local = i % nlocal;   // Define column submatrix
            MPI_Reduce(&local_sum, &solution[i_local], 1, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            
            if (rank == root){      // DIAGONAL ELEMENT
                i_local = i % nlocal;
                index_mtx = i * nlocal + i_local; // Indice da matriz
                solucao = (vector[i_local] - solution[i_local]) / submatrix[index_mtx];
                delta_x = get_max(delta_x, fabsl(solution[i_local] - solucao));
                solution[i_local] = solucao;
            }
        }
        // iterations++;
        MPI_Allreduce(&delta_x , &global_delta, 1 , MPI_LONG_DOUBLE , MPI_MAX , MPI_COMM_WORLD);
        // MPI_Allreduce(&iterations , &global_iter, 1 , MPI_INT , MPI_MAX , MPI_COMM_WORLD);
    
    // STOPPED CONDITION
    } while (fabsl(global_delta - global_delta_previous) > tolerance );
    
    end = MPI_Wtime();

    if (rank==0){
        printf("time elapsed during the job: %.2f miliseconds.\n", (end - start) * 1000);
        if (flag==1){
            for(int i = 0; i < dimension; i++){
                printf("X%i = %LF\n", i+1, solution[i]);
            }
        }
        filename = "results.log";
        FILE *outfile = fopen(filename, "a");
        if (outfile == NULL){
            printf("Error opening file!\n");
            exit(1);
        }
        fprintf(outfile, "Matriz de tamanho %i Tolerance %50.50Lf \n", dimension,tolerance);
        fprintf(outfile, "Tempo de execucao: %lf milisegundos\n", (end - start) * 1000);
        fprintf(outfile, "\n");
    }
    MPI_Finalize();
    free(submatrix);
    free(solution);
    free(vector);
}
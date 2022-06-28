/* Program Gauss-Seidel
   Solution of a system of linear equations by Gauss-Seidel's
   iteration method. Assume that the coefficient matrix satisfies
   the condition of convergence.

   INSTRUCAO PARA OSX
Para evitar a msg:
A system call failed during shared memory initialization that should
not have.  It is likely that your MPI job will now either abort or
experience performance degradation.

digite export TMPDIR=/tmp no seu terminal

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
#include "matrixtools.h"

// #define max(a, b) ((a) > (b) ? (a) : (b))

int main(int argc, char *argv[])
{
    system("clear");
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    MPI_Datatype subarray_type;

    double start, end; // To compute running time

    int rank, size;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank); 

    if (argc < 3){
        printf("Usage: output [num_decimals] [filename] [0-1]\n");
        MPI_Finalize();
        exit(1);
    }

    clock_t t; //variável para armazenar tempo

	long double tolerance = 1.000, sum, temp, auxiliary, global_delta;
    long double local_sum, delta_x, solucao;

    char *filename;
    FILE *infile;
	
    int offset, flag=0, i, j, k, n, aux, index_mtx=0;

    int root=0, i_local=0, nlocal=0;

    int dimension = 0, nline_submatrix=0, nelements_submatrix=0;
    
    if (rank == 0){
        // printf("1. Preparando variáveis.\n");
        // Read command-line arguments
        i = (-1) * atoi(argv[1]);       
        filename = argv[2];
        if (argv[3] != NULL)
            if (atoi(argv[3]) == 1)
                flag = atoi(argv[3]);
        tolerance = tolerance * pow(10, i);

        dimension = read_dimension(filename);
        
        // Set submatrix' size
        nline_submatrix = (int)dimension / size;
        // Adjusting if round it down
        if (nline_submatrix * size < dimension) nline_submatrix++;
        nelements_submatrix = nline_submatrix * dimension;
    }

    MPI_Bcast( &dimension , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &tolerance , 1 , MPI_LONG_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    nelements_submatrix = get_size_submatrix(rank, size, nline_submatrix, dimension);
    nline_submatrix = get_num_lines_submtx(rank, size, dimension);

    long double *vector = calloc(nline_submatrix , sizeof(long double));
    long double *submatrix = calloc(nelements_submatrix, sizeof(long double));
    long double *solution = calloc(nline_submatrix , sizeof(long double));

    if (rank == 0){
        printf("2.Montando matrizes, rank %i DIMENSAO %i!\n", rank, dimension);

        long double *matrix = calloc(dimension * dimension , sizeof(long double));
        long double *vector_aux = calloc(dimension, sizeof(long double));

        read_matrix(filename, matrix, vector_aux);

        // printf("VETOR CONSTANTE LIDO\n");
        // for(int hh=0; hh < dimension; hh++)
        //     printf("indice %i {{%Lf}} ", hh, vector_aux[hh]);
        int linha, coluna;
        for (i = 0;i < size; i++){
            offset = matrix_offset(i, size, dimension);

            long double *auxiliary = calloc(offset * dimension, sizeof(long double));
            int count = 0;
            int lower_limit, upper_limit;
            for (linha=0; linha < dimension; linha++){
                lower_limit=offset * i;
                upper_limit=offset * ( i + 1);
                for (coluna=lower_limit; coluna < upper_limit; coluna++){
                    int index = linha * dimension + coluna;
                    if (i!=0)
                        auxiliary[count] = matrix[index];
                    else
                        submatrix[count] = matrix[index];
                    count++;
                }
            }
            if (i != 0){
                MPI_Isend(auxiliary, dimension * offset, MPI_LONG_DOUBLE, i, 3000, MPI_COMM_WORLD, &request);
                int index = offset * i;
                int nelements = get_num_lines_submtx(i, size, dimension);
                MPI_Send(&vector_aux[index], nelements, MPI_LONG_DOUBLE, i, 3001, MPI_COMM_WORLD);
            } else {
                for (linha = 0; linha < nline_submatrix; linha++){
                    vector[linha] = vector_aux[linha];
                }
            }
            free(auxiliary);
        }
        free(vector_aux);
        free(matrix);
    } 
    if (rank != 0){
        MPI_Recv(submatrix, nelements_submatrix, MPI_LONG_DOUBLE, 0, 3000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        int nelements = get_num_lines_submtx(i, size, dimension);
        MPI_Recv(&vector[0], nelements, MPI_LONG_DOUBLE, 0, 3001, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0){
        printf("****Waiting! Calculating! \n" );
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    nlocal = matrix_offset(rank, size, dimension);
    printf("Rank %i size %i offset %i\n", rank, size, nlocal);
    solucao=0.0;
    delta_x = 0.0;
    start = MPI_Wtime();
    int iterations = 0;
    int global_iter =0;
    delta_x = 0.0;
    do {
        delta_x = 0.0;
        for (i = 0; i < dimension; i++){
            local_sum = 0.0;
            for (j = 0; j < nlocal; j++){
                if (j + rank * nlocal != i){
                    index_mtx = i * nlocal + j; // Indice da matriz
                    local_sum = local_sum + submatrix[index_mtx] * solution[j];
                }
            }            
            root = i / nlocal;      // Define submatrix number 
            i_local = i % nlocal;   // Define column submatrix
            MPI_Reduce(&local_sum, &solution[i_local], 1, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            if (rank == root){
                index_mtx = i * nlocal + i_local; // Indice da matriz
                solucao = (vector[i_local] - solution[i_local]) / submatrix[index_mtx];
                delta_x = get_max(delta_x, fabsl(solution[i_local] - solucao));
                solution[i_local] = solucao;
            }
        }
        iterations++;
        MPI_Allreduce(&delta_x , &global_delta, 1 , MPI_LONG_DOUBLE , MPI_MAX , MPI_COMM_WORLD);
        MPI_Allreduce(&iterations , &global_iter, 1 , MPI_INT , MPI_MAX , MPI_COMM_WORLD);
    } while (global_delta > tolerance && global_iter < 10000);
    
    end = MPI_Wtime();
    printf("time elapsed during the job: %.2f miliseconds.\n", (end - start) * 1000);

    if (rank==0){
        if (flag==1){
            for(int i = 0; i < dimension; i++){
                printf("X%i = %LF\n", i+1, solution[i]);
            }
        }
        filename = "time.slt";
        FILE *outfile = fopen(filename, "w");
        if (outfile == NULL){
            printf("Error opening file!\n");

            exit(1);
        }
        fprintf(outfile, "Matriz de tamanho %i Tolerance %50.50Lf \n", dimension,tolerance);
        fprintf(outfile, "Tempo de execucao: %lf milisegundos\n", (end - start) * 1000);
    }
    MPI_Finalize();
    free(submatrix);
    free(solution);
    free(vector);
}
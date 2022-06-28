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

#define max(a,b) ((a) > (b) ? (a) : (b))


// PARA USLEEP
#include <unistd.h>

// #include "determinante.c"
// float determinante(double **A, int nvar);

#include "functools.h"
#include "matrixtools.h"

int main(int argc, char *argv[])
{
    system("clear");
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Request request;
    MPI_Datatype subarray_type;


    double start, end;

    int rank, size;
    MPI_Comm_size (MPI_COMM_WORLD, &size);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank); 

    if (argc < 3){
        printf("Usage: output [num_decimals] [filename] [0-1]\n");
        MPI_Finalize();
        exit(1);
    }

    clock_t t; //variável para armazenar tempo

	long double tolerance, sum, temp, auxiliary, global_delta;
    long double local_sum, delta_x, s_k, solucao;

    char *filename;
    FILE *infile;
	
    int offset, flag=0, i, j, k, n, aux, index_mtx=0;

    int root=0, i_local=0, nlocal=0;

    int nvar = 0, nline_submatrix=0, nelements_submatrix=0;
    
    if (rank == 0){
        // printf("1. Preparando variáveis.\n");
        // Read command-line arguments
        i = (-1) * atoi(argv[1]);       
        filename = argv[2];
        if (argv[3] != NULL)
            if (atoi(argv[3]) == 1)
                flag = atoi(argv[3]);
        tolerance = tolerance * pow(10, i);

        nvar = read_dimension(filename);
        
        // Set submatrix' size
        nline_submatrix = (int)nvar / size;
        // Adjusting if round it down
        if (nline_submatrix * size < nvar) nline_submatrix++;
        nelements_submatrix = nline_submatrix * nvar;
    }

    MPI_Bcast( &nvar , 1 , MPI_INT , 0 , MPI_COMM_WORLD);
    MPI_Bcast( &tolerance , 1 , MPI_LONG_DOUBLE , 0 , MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    
    nline_submatrix = get_num_lines_submtx(rank, size, nvar);
    nelements_submatrix = get_size_submatrix(rank, size, nline_submatrix, nvar);
    // printf("rank %i Elementos da submatrix %i\n", rank, nelements_submatrix);
    
    long double *vector = calloc(nvar , sizeof(long double));
    long double *submatrix = calloc(nelements_submatrix, sizeof(long double));
    long double *solution = calloc(nvar , sizeof(long double));

    if (rank == 0){
        printf("2.Montando matrizes, rank %i DIMENSAO %i!\n", rank, nvar);

        long double *matrix = calloc(nvar * nvar , sizeof(long double));
        read_matrix(filename, matrix, vector);

        // for (i = 1; i < size; i++){
        //     MPI_Send(vector, nvar, MPI_LONG_DOUBLE, i, 123, MPI_COMM_WORLD);
        // }
        MPI_Bcast(&vector, nvar, MPI_LONG_DOUBLE, 0, MPI_COMM_WORLD);
        int linha, coluna;
        for (i = 0;i < size; i++){
            // int dimensions_full_array[2] = {nvar, nvar};
            // int dimensions_subarray[2] = {nvar, o};
            // int start_coordinates[2] = {0, offset};

            // offset = (int)nvar/size;
            offset = matrix_offset(i, size, nvar);

            long double *auxiliary = calloc(offset * nvar, sizeof(long double));
            int count = 0;
            int lower_limit, upper_limit;
            for (linha=0; linha < nvar; linha++){
                lower_limit=offset * i;
                upper_limit=offset * ( i + 1);
                for (coluna=lower_limit; coluna < upper_limit; coluna++){
                    int index = linha * nvar + coluna;
                    if (i!=0)
                        auxiliary[count] = matrix[index];
                    else
                        submatrix[count] = matrix[index];
                    // printf("(rank %i) valor %Lf\n", i, submatrix[count]);
                    count++;
                }
            }
            if (i != 0){
                printf("ENVIANDO do %i para %i\n", rank, i);
                printf("Enviando %i elementos\n", nvar * offset);
                MPI_Isend(auxiliary, nvar * offset, MPI_LONG_DOUBLE, i, 3000, MPI_COMM_WORLD, &request);
            }
            free(auxiliary);
        }
        free(matrix);
    } else {
    // if (rank != 0){
        // MPI_Recv(vector, nvar, MPI_LONG_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("RECEBENDO %i\n", rank);
        printf("recebendo %i elementos\n", nelements_submatrix);
        MPI_Recv(submatrix, nelements_submatrix, MPI_LONG_DOUBLE, 0, 3000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("RECEBIDO %i\n", rank);
    }
    if (rank == 0){
        printf("****Waiting! Calculating! \n" );
    }
        printf("SUBMATRIZ ------------%i\n", rank);
    for (i=0; i<nelements_submatrix;i++)
        printf("(%i)[%i]%Lff ", rank,i,submatrix[i]);
    // MPI_Barrier(MPI_COMM_WORLD);
    nlocal = matrix_offset(rank, size, nvar);
    // nlocal = nvar/size ;
    printf("Rank %i size %i offset %i\n", rank, size, nlocal);
    // delta_x = 0.0;
    // s_k =0.0;
    solucao=0.0;

    start = MPI_Wtime();
    int iterations = 0;
    int global_iter =0;
    do {
        printf("%Lf ", solution[0]);
        // delta_x = 0.0;
        for (i = 0; i < nvar; i++){
            s_k = 0.0;
            for (j = 0; j < nlocal; j++){
                if (j + rank * nlocal != i){
                    index_mtx = i * nlocal + j; // Indice da matriz
                    s_k = s_k + submatrix[index_mtx] * solution[j];
                }
            }            
            root = i / nlocal;
            i_local = i % nlocal;
            MPI_Reduce(&s_k, &solution[i_local], 1, MPI_LONG_DOUBLE, MPI_SUM, root, MPI_COMM_WORLD);
            printf("SOLUTION %i After reduction %Lf global_delta %Lf global_iter%i\n", i_local, solution[i_local], global_delta, global_iter);
            if (rank == root){
                index_mtx = i_local * nlocal + i_local; // Indice da matriz
                solucao = (vector[i_local] - solution[i_local]) / submatrix[index_mtx];
                delta_x = max(delta_x, fabsl(solution[i_local] - solucao));
                solution[i_local] = solucao;
            }
        }
        iterations++;
        MPI_Allreduce(&delta_x , &global_delta, 1 , MPI_LONG_DOUBLE , MPI_MAX , MPI_COMM_WORLD);
        MPI_Allreduce(&iterations , &global_iter, 1 , MPI_INT , MPI_MAX , MPI_COMM_WORLD);
    } while (global_delta > tolerance && global_iter < 12000);
    
    end = MPI_Wtime();
    // end = clock();
    printf("time elapsed during the job: %.2f miliseconds.\n", (end - start) * 1000);


    if (rank==0){
        if (flag==1){
            for(int i = 0; i < nvar; i++){
                printf("X%i = %LF\n", i+1, solution[i]);
            }
        }
        filename = "time.slt";
        FILE *outfile = fopen(filename, "w");
        if (outfile == NULL){
            printf("Error opening file!\n");

            exit(1);
        }
        fprintf(outfile, "Matriz de tamanho %i Tolerance %50.50Lf \n", nvar,tolerance);
        fprintf(outfile, "Tempo de execucao: %lf milisegundos\n", (end - start) * 1000);
    }
        MPI_Finalize();
    free(submatrix);
    free(solution);
    free(vector);
}
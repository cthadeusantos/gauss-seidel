/* Program Gauss-Seidel
   Solution of a system of linear equations by Gauss-Seidel's
   iteration method. Assume that the coefficient matrix satisfies
   the condition of convergence.
 
    nvar            dimension of the matrix
    num_thread      the number of threads will be used
    size_per_thread the number of elements calculates by thread
    tolerance       precision error
    matrix          matrix A (coefficient matrix)
    vector          matrix B (constant matrix)
    solution        matrix x (variable matrix)
    previous        matrix with previous values from solution matrix (initial values are zero)
*/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    if (argc != 4){
        printf("Usage: output [num_decimals] [filename] [threads]\n");
        exit(1);
    }

    clock_t t; //variável para armazenar tempo

	long double tolerance=1, sum, auxiliary;

	int i, j, flag, count, nvar = 0;

    char *filename;
    FILE *infile;
    
    i = (-1) * atoi(argv[1]);       
    filename = argv[2];
    
    int num_threads;
    num_threads = atoi(argv[3]);

    tolerance = tolerance * pow(10, i);
    
    infile=fopen(filename,"rt");
    if (NULL == infile){
        printf("file can't be opened \n");
    }
    

    // READ MATRIX FILE
    fscanf(infile,"%d",&nvar);    // Read matrix size

    printf("Fill matrix and arrays!\n");

    long double **matrix = calloc(nvar , sizeof(long double*)); 
    long double *vector = calloc(nvar , sizeof(long double));   // Declare an array to store independent set
    long double *previous = calloc(nvar , sizeof(long double)); // Declare an auxiliary array
    long double *solution = calloc(nvar , sizeof(long double)); // Declare an array to store variables' values


    // build matrix columns
    for (i = 0; i < nvar; i++){
        matrix[i] = calloc(nvar, sizeof(long double));
    }

    printf("Reading matrix file!\n");
    count = 0;
    long double element[3];
    while (fscanf(infile,"%Lf", &auxiliary) == 1)
    {
        element[count] = auxiliary;
        count++;
        if (count == 3){
            int row = (int)element[0] - 1;
            int col = (int)element[1] - 1;
            
            if (col < nvar){
                matrix[row][col] = element[2] ;
            } else {
                vector[row] = element[2];
            }
            count = 0;
        }
    }
    fclose(infile);


    // INITIALIZE MATRIXs
	for(i = 0; i < nvar; i++){
		previous[i]=0.0; //initialize
        solution[i] = 0.0;
    }

    // Gauss-seidel parallel method 
    printf("Waiting! Calculating!\n");
    
    omp_set_num_threads(num_threads);

    t = clock();
    do {
        flag = 0;
        for (i = 0; i < nvar; i++){
            // sum = vector[i];
            sum = 0;
            {
                #pragma omp parallel for schedule(guided) num_threads(num_threads) shared(i, matrix, previous) reduction ( +: sum )
                    for (j = 0; j < nvar; j++){
                        if (j != i){
                            sum = sum + matrix[i][j] * previous[j];
                        }
                    }
            }
            sum = vector[i] - sum;
            solution[i] = sum / matrix[i][i];
            if (fabsl(solution[i] - previous[i])/solution[i] > tolerance){
                flag = 1;
            }
            previous[i] = solution[i];
        }
    } while(flag == 1);

    t = clock() - t;

    printf("Tempo de execucao: %lf milisegundos\n", ((double)t)/((CLOCKS_PER_SEC/1000))); //conversão para double
    printf("Results saved at solution.slt file!\n");
    for (i = 0; i < nvar; i = i + (int)(nvar / 10) + 1){
        printf("X%i = %30.30Lf \n", i+1, solution[i]);
    }
    filename = "solution.slt";
    FILE *outfile = fopen(filename, "w");
    if (outfile == NULL){
        printf("Error opening file!\n");
        exit(1);
    }
    fprintf(outfile, "Matriz de tamanho %i Tolerance %50.50Lf \n", nvar,tolerance);
    fprintf(outfile, "Tempo de execucao: %lf milisegundos\n", ((double)t)/((CLOCKS_PER_SEC/1000)));
    fprintf(outfile, "Solution values to checking:");
    for (i = 0; i < nvar; i = i + (int)(nvar / 10) + 1){
        fprintf(outfile, "X%i = %30.30Lf \n", i+1, solution[i]);
    }

    fclose(outfile);

    free(matrix);
    free(vector);
    free(solution);
    free(previous);
}

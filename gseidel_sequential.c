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
    previous        matrix with previous values from solution matrix (initial values are zero)
    flag            tolerance reached
   */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <stdbool.h>
#include <time.h>

int main(int argc, char *argv[])
{
    if (argc != 3){
        printf("Usage: output [num_decimals] [filename]\n");
        exit(1);
    }

    clock_t t; //variável para armazenar tempo

    /*
     *  Define variables
     */
	long double tolerance=1, sum, auxiliary;
	int flag, i, j, k, n;
    char *filename;
    int nvar = 0;
    FILE *infile;
    
    i = (-1) * atoi(argv[1]);
    filename = argv[2];
    
    tolerance = tolerance * pow(10, i);
    
    /*
    * READ MATRIX DIMENSION
    */
    infile=fopen(filename,"rt");
    if (NULL == infile){
        printf("file can't be opened \n");
    }
    fscanf(infile,"%d",&nvar);    // Read matrix dimension

    /*
    * Create matrix and vectors
    */
    long double **matrix = calloc(nvar , sizeof(long double*)); // Declare matrix A
    long double *vector = calloc(nvar , sizeof(long double));   // Declare an array to store independent set
    long double *previous = calloc(nvar , sizeof(long double)); // Declare an auxiliary array with previous solutions
    long double *solution = calloc(nvar , sizeof(long double)); // Declare an array to store variables' values

    // build matrix columns
    for (i = 0; i < nvar; i++){
        matrix[i] = calloc(nvar, sizeof(long double));
    }

    int count = 0;          // auxiliary counter
    long double element[3]; // elements: row col value
    
    /*
    * READ MATRIX
    */
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

	for(i = 0; i < nvar; i++){
		previous[i]=0.0;    //initialize
        solution[i] = 0.0;  // initialize
    }

    printf("Waiting! Sequential calculating!\n");

    /*
     * Gauss-Seidel method
     */
    t = clock();    // Start time
	do {
        flag = 0;   // tolerance reached
        for (i = 0; i < nvar; i++){
            sum = vector[i];
            for (j = 0; j < nvar; j++){
                if (j != i){    // Element is not diagonal
                    sum = sum - matrix[i][j] * previous[j]; // sum
                }
            }

            solution[i] = sum / matrix[i][i]; // Calculate new solution like GS recommendation
            
            if (fabsl(solution[i] - previous[i])/solution[i] > tolerance){
                flag = 1;   // If any solution not reached tolerance
            }
            previous[i] = solution[i];
        }
	} while(flag == 1);
    t = clock() - t;    // End time

    printf("Tempo de execucao: %lf milisegundos\n", ((double)t)/((CLOCKS_PER_SEC/1000))); //conversão para double

    filename = "solution.slt";
    FILE *outfile = fopen(filename, "w");
    if (outfile == NULL){
        printf("Error opening file!\n");
        exit(1);
    }
    fprintf(outfile, "Matriz de tamanho %i Tolerance %50.50Lf \n", nvar,tolerance);
    fprintf(outfile, "Tempo de execucao: %lf milisegundos\n", ((double)t)/((CLOCKS_PER_SEC/1000)));
    // fprintf(outfile, "Solution values to checking:");
    // for (i = 0; i < nvar; i = i + (int)(nvar / 10) + 1){
    //     fprintf(outfile, "X%i = %30.30Lf \n", i+1, solution[i]);
    // }

    fclose(outfile);

    free(matrix);
    free(vector);
    free(solution);
    free(previous);
}

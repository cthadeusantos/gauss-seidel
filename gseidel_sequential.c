/* Program Gauss-Seidel
   Solution of a system of linear equations by Gauss-Seidel's
   iteration method. Assume that the coefficient matrix satisfies
   the condition of convergence.
   */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include <stdbool.h>
#include <time.h>
// #include "determinante.c"
// float determinante(double **A, int nvar);

int main(int argc, char *argv[])
{
    if (argc != 3){
        printf("Usage: output [num_decimals] [filename]\n");
        exit(1);
    }

    clock_t t; //variável para armazenar tempo

	long double tolerance=1, sum, temp, auxiliary;
	// double x[10], solution[10], tolerance=1.0, sum;

	int flag, i, j, k, n;

    char *filename;
    int nvar = 0;
    FILE *infile;
    
    i = (-1) * atoi(argv[1]);       
    filename = argv[2];
    
    tolerance = tolerance * pow(10, i);
    printf("tolerance %i decimals, value %20.15Lf\n", i*(-1), tolerance);
    
    infile=fopen(filename,"rt");
    if (NULL == infile){
        printf("file can't be opened \n");
    }
    
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
    int count = 0;
    long double element[3];
    while (fscanf(infile,"%Lf", &auxiliary) == 1)
    {
        // printf("%Lf\n", auxiliary);
        element[count] = auxiliary;
        count++;
        if (count == 3){
            int row = (int)element[0] - 1;
            int col = (int)element[1] - 1;
            
            if (col < nvar){
                // printf("<%i %i>\n", row, col);
                matrix[row][col] = element[2] ;
            } else {
                // printf("<%i> %i\n", row, col);
                vector[row] = element[2];
            }
            count = 0;
        }
    }
    fclose(infile);

	for(i = 0; i < nvar; i++){
		previous[i]=0.0; //initialize
        solution[i] = 0.0;
    }

    printf("Waiting! Calculating!\n");

    t = clock();
	do {
        flag = 0;

        for (i = 0; i < nvar; i++){
            sum = vector[i];
            for (j = 0; j < nvar; j++){
                if (j != i){
                    sum = sum - matrix[i][j] * previous[j];
                }
            }
            solution[i] = sum / matrix[i][i];
            if (fabsl(solution[i] - previous[i])/solution[i] > tolerance){
                flag = 1;
            }
            previous[i] = solution[i];
        }

	} while(flag == 1);
    t = clock() - t;

	// printf("Solution is \n");
	// for(i = 0;i < nvar ; i++)
	// 	printf("%30.30Lf %30.30Lf\n",solution[i], previous[i]);
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

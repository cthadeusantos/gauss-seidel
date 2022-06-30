#include<stdio.h>

int read_dimension(char *filename){
    int nvar;
    FILE *infile;
    
    infile=fopen(filename,"rt");

    if (NULL == infile)
        printf("file can't be opened \n");
    
    fscanf(infile,"%d",&nvar);    // Read matrix size
    fclose(infile);
    
    infile = NULL;
    return nvar;
}

int read_matrix(char *filename, long double *matrix, long double *vector){
        int dimension;
        long double auxiliary;
        FILE *infile;
        infile=fopen(filename,"rt");
        if (NULL == infile){
            printf("file can't be opened \n");
        }
        
        fscanf(infile,"%d",&dimension);    // Read matrix size

        int count = 0;
        int block = 0;
        long double element[3];
        while (fscanf(infile,"%Lf", &auxiliary) == 1) {
            element[count] = auxiliary;
            count++;
            if (count == 3){
                int row = (int)element[0] - 1;
                int col = (int)element[1] - 1;
                
                if (col < dimension){
                    matrix[dimension * row + col] = element[2] ;
                    // printf("elemento %Lf \n", element[2]);
                } else {
                    vector[row] = element[2];
                }
                count = 0;
            }
        }
        fclose(infile);
        infile = NULL;
    return 0;
}
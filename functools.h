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

int get_size_submatrix(int rank, int size, int nline_submatrix, int dimension){
    if (rank != (size - 1)){
        return nline_submatrix * dimension;
    } else {
            // Set submatrix' size
        int offset = (int)dimension / size;
        // Adjusting if round it down
        if (offset * size < dimension) offset++;
        return (dimension * dimension) - ((size  - 1) * dimension * (offset));
    }
}

int num_elements_matrix(int dimension){
    return dimension * dimension;
}

int get_num_lines_submtx(int rank, int size, int dimension){
    // Set submatrix' size
    int nline_submatrix = (int)dimension / size;
    // Adjusting if round it down
    if (nline_submatrix * size < dimension) nline_submatrix++;
    if (rank != (size - 1)){
        return nline_submatrix;
        // nelements_submatrix = nline_submatrix * nvar;
    } else {
        return dimension - (nline_submatrix * rank);
    }
}

int matrix_offset(int rank, int size, int dimension){
    // Set submatrix' size
    int nline_submatrix = (int)dimension / size;
    // Adjusting if round it down
    if (nline_submatrix * size < dimension) nline_submatrix++;
    return nline_submatrix;
}

int get_size_subvector(int rank, int size, int dimension){
    int vector_size = (int)dimension / size;
    // Adjusting if round it down
    if (vector_size * size < dimension) vector_size++;
    if (rank != (size - 1)){
        return vector_size;
        // nelements_submatrix = nline_submatrix * nvar;
    } else {
        return dimension - (vector_size * rank);
    }
}

int get_offset_subvector(int rank, int size, int dimension){
    int offset = (int)dimension / size;
    // Adjusting if round it down
    if (offset * size < dimension) offset++;
    // if (rank != (size - 1)){
        return offset;
        // nelements_submatrix = nline_submatrix * nvar;
    // } else {
        // return dimension - (offset * rank);
    // }  
}

long double get_max(long double a, long double b){
    if (a > b)
        return a;
    return b;
}

long double get_min(long double a, long double b){
    if (a < b)
        return a;
    return b;
}

float determinante(double **A, int nvar){
    int matriz[nvar][nvar];
    int fdr = nvar, i , j, k;
    for (i = 0; i < nvar; i++){
        for (j = 0; j < nvar; j++){
            matriz[i][j] = A[i][j];
        }
    }
    float mult;
    float deter=1;
    for(i=0;i<fdr;i++)
    {
        for(j=0;j<fdr;j++)
        {
            mult=matriz[j][i]/matriz[i][i];
            for(k=0;k<fdr;k++)
            {
                if(i==j) break;
                    matriz[j][k]=matriz[j][k]-matriz[i][k]*mult;
            }
        }
    }
    for(i=0;i<fdr;i++)
        {
            deter=deter*matriz[i][i];
        }
    return deter;
}


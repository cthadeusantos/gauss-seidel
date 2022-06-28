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
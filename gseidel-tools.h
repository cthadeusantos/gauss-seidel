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
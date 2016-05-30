#include "LUdecomp.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

LUdecomp *LUdecompose(int N, double **A){ // Performs the decomposition

  LUdecomp *decomp = (LUdecomp *) malloc(sizeof(LUdecomp)); // Allocate memory
  decomp->N = N; // Store N
  decomp->LU = (double **) malloc(N*sizeof(double*)); // Allocate memories
  for(int i = 0; i < N; i++){
    decomp->LU[i] = (double *) malloc(N*sizeof(double));
  }
  decomp->mutate = (short *) malloc(N*sizeof(short));
  decomp->d = 1; // initialize row swap parity value

  short *mutate = decomp->mutate; // Alias for LUdecomp->mutate
  double **LU = decomp->LU; // Alias for LUdecomp->LU  
  for(int r = 0; r < N; r++) // Copy A into LUdecomp->LU
    for(int c= 0; c < N; c++)
      LU[r][c] = A[r][c];

  // Code from LUdecomp pdf
  
  for(int i = 0; i <= N-1; i++) // Initialize row permutation array
    mutate[i] = i;
  
  for(int j = 0; j <= N-1; j++){ // Replace A with LU column by column
    for(int i = 0; i <= j; i++){ // Compute aij <- Bij on and above diagonal
      double sum = 0.0;
      for (int k = 0; k <= i-1; k++)
	sum += LU[i][k] * LU[k][j];
      LU[i][j] = LU[i][j] - sum;
    }
    
    double p = fabs(LU[j][j]); // p = initial pivot value
    int n = j; // n = initial pivot row
    
    for(int i = j+1; i <= N-1; i++){  // Compute aij <- aij below diagonal
      double sum = 0.0;
      for(int k = 0; k <= j-1; k++)
	sum += LU[i][k] * LU[k][j];
      LU[i][j] = LU[i][j] - sum;
      if(fabs(LU[i][j]) > p){ // If better pivot found
	p = fabs(LU[i][j]); // ...then record new pivot
	n = i;
      }
    }
    if(p == 0.0) // If p == 0 abort
      break;
    if(n != j){ // If best pivot off diagonal
      double* rowTemp = LU[n]; // swap rows n and j
      LU[n] = LU[j];
      LU[j] = rowTemp;
      
      short mutateTemp = mutate[n]; // swap mutate[n] and mutate[j]
      mutate[n] = mutate[j];
      mutate[j] = mutateTemp;
      
      decomp->d = -(decomp->d); //flip parity
    }
    for(int i = j+1; i <= N-1; i++){ // perfrom divisions below the diagonal
      LU[i][j] = LU[i][j] / LU[j][j];
    }
  }
  decomp->mutate = mutate;
  decomp->LU = LU;
  return decomp;
}

void LUdestroy(LUdecomp *decomp){ //deallocates the data that was allocated
  const int N = decomp->N;
  free(decomp->mutate);
  for(int r = N-1; r >= 0; r--)
    free(decomp->LU[r]);
  free(decomp->LU);
  free(decomp);
}

void LUsolve(LUdecomp *decomp, double *b, double *x){ //solves the system Ax - b for x

  int N = decomp->N;
  double y[N]; // initialize y
  
  y[0] = b[decomp->mutate[0]]; // forward substitution
  for(int i = 1; i <= N-1; i++){
    double sum = 0.0;
    for(int j = 0; j <= i-1; j++)
      sum += decomp->LU[i][j] * y[j];
    y[i] = b[decomp->mutate[i]] - sum;
  }
  
  x[N-1] = y[N-1] / decomp->LU[N-1][N-1]; // back substitution
  for(int i = N-2; i >= 0; i--) {
    double sum = 0.0;
    for(int j = i+1; j <= N-1; j++)
      sum += decomp->LU[i][j] * x[j];
    x[i] = (y[i] - sum)/decomp->LU[i][i];
  }
}

#include "LUdecomp.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void solveAndPrint(int N, double **AMatrix, double *bMatrix, double *xMatrix){
  int j = 0;
  LUdecomp *LUdecomp = malloc(sizeof(LUdecomp));
  LUdecomp = LUdecompose(2*N, AMatrix); // perform decomp
  LUsolve(LUdecomp, bMatrix, xMatrix); // solve x
  for(int i = 0; i < 2*N; i++){ // print to stdout
    printf("%.10f ", xMatrix[i]); //print
    j++;
    if(j == 3){ // 3 values, move to next row
      j = 0;
      printf("\n");
    }
  }
  printf("1.0000000000"); // print 1 at the end
}

int main (int argc, char *argv[]){
  int N, j = 0;
  double **AMatrix, *bMatrix, *xMatrix, **xy, **xyPrime, x, y, xPrime, yPrime;
  
  scanf("%d", &N); // read first line of file from stdin and assign to N

  /* allocating space for various matrices */

  AMatrix = malloc(2*N*sizeof(double*)); // allocate space for matrix A
  for(int i = 0; i < 2*N; i++)
    AMatrix[i] = malloc(8*sizeof(double));
  bMatrix = malloc(2*N*sizeof(double)); // allocate space for b
  xMatrix = malloc(2*N*sizeof(double)); //allocate space for x
  xy = malloc(2*sizeof(double*)); // allocate space for the x and y values
  xy[0] = malloc(N*sizeof(double));
  xy[1] = malloc(N*sizeof(double));
  xyPrime = malloc(2*sizeof(double*)); // allocate space for the x' and y' values
  xyPrime[0] = malloc(N*sizeof(double));
  xyPrime[1] = malloc(N*sizeof(double));
  
  for(int i = 0; i <= N-1; i++) // store values for x and y for equation 7
    scanf("%lf %lf", &xy[0][i], &xy[1][i]);
  for(int i = 0; i <= N-1; i++) // store values for x' and y' for equation 7
    scanf("%lf %lf", &xyPrime[0][i], &xyPrime[1][i]);

  for(int i = 0; i < 2*N; i++){ // equation 7
    x = xy[0][j];
    y = xy[1][j];
    xPrime = xyPrime[0][j];
    yPrime = xyPrime[1][j];
    
    AMatrix[i][0] = x;
    AMatrix[i][1] = y;
    AMatrix[i][2] = 1;
    AMatrix[i][3] = 0;
    AMatrix[i][4] = 0;
    AMatrix[i][5] = 0;
    AMatrix[i][6] = -x*xPrime;
    AMatrix[i][7] = -y*xPrime;
    bMatrix[i] = xPrime;
    i++; // next row
    
    AMatrix[i][0] = 0;
    AMatrix[i][1] = 0;
    AMatrix[i][2] = 0;
    AMatrix[i][3] = x;
    AMatrix[i][4] = y;
    AMatrix[i][5] = 1;
    AMatrix[i][6] = -x*yPrime;
    AMatrix[i][7] = -y*yPrime;
    bMatrix[i] = yPrime;
    j++; // x, y, x', y', changes every 2 rows
  }
  solveAndPrint(N, AMatrix, bMatrix, xMatrix);
}

#ifndef DFT_CPP
#define DFT_CPP

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include "globals.h"

using namespace std;


double computeNorm(gsl_vector_complex* a, gsl_vector_complex* b, int& N )
{
  double norm, A, B = 0;
  for(int i = 0; i<N; i++)
  {
    A = pow(GSL_REAL(gsl_vector_complex_get(a,i)),2) + pow(GSL_IMAG(gsl_vector_complex_get(a,i)),2);
    B = pow(GSL_REAL(gsl_vector_complex_get(b,i)),2) + pow(GSL_IMAG(gsl_vector_complex_get(b,i)),2);
    A = sqrt(A);
    B = sqrt(B);
    norm = max(norm, abs(A, B));
  }
  
  return norm;
}  

int iterativeSolver(int N)
{
  bool ok =false;
  int success = 0;
  gsl_vector_complex * X = gsl_vector_complex_alloc(N);
  gsl_vector_complex * XO = gsl_vector_complex_calloc(N); //initialize vector of size and make all elements 0
  gsl_complex result; gsl_complex mult;
  
  //STEP 1
  int k = 1;
  
  //STEP 2
   while(k<=N)
   {
   //STEP 3
     for(int i=1; i<=N; i++)
     {
       result = gsl_complex_rect(0, 0);
       for(int j=1; j<=n; j++)
       {
         if(j!=i)
         {
           mult = gsl_complex_mul(gsl_matrix_complex_get(Z, i-1, j-1), gsl_vector_complex_get(XO, j-1));
           result = gsl_complex_add(result, mult);
         }
       }
       result = gsl_complex_mul_real(result, -1.0);
       result = gsl_complex_add(result, gsl_vector_complex_get(GC, i-1); 
       result = gsl_complex_div(result, gsl_vector_complex_get(Z, i-1, i-1)
       gsl_vector_complex_set(X, i-1, result);
     }
     
     //STEP 4
     if(computeNorm(X, XO, N) < TOL)
     {
       ok = true;
       gsl_vector_complex_memcpy(X, Y);
       break;
     }
     
     //Step 5
     k++;
     
     //Step 6
     gsl_vector_complex_memcpy(XO, X);
   } 
   if(ok ==false)
   {
     success = 1;
     inverseSolver();
   }
   return success;
}

void directSolver(gsl_matrix_complex* Z, gsl_vector_complex* GC, gsl_vector_complex* y, const int& n)
{
     
    int s; 
    gsl_permutation * p = gsl_permutation_alloc(n);
    gsl_linalg_complex_LU_decomp(Z, p, &s);
    gsl_linalg_complex_LU_solve(Z, p, GC, y);
    gsl_permutation_free(p);
}


void inverseSolver(gsl_matrix_complex* Z, gsl_vector_complex* c, gsl_vector_complex* y, const int& n)
{
        //compute Z inverse by finding the complex conjugate of all the entries in Z
       gsl_matrix_complex * ZI = gsl_matrix_complex_alloc(n, n);
       gsl_complex tempZI;
       for(int j = 0; j<n; j++)
      {
        for(int k = 0; k<n; k++)
        {
          tempZI = gsl_complex_conjugate(gsl_matrix_complex_get(Z, j, k));
          gsl_matrix_complex_set(ZI, j, k, tempZI); 
        }
    }
    
    //retrieve the filtered y points by multiplying Zinverse by c
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, ZI, c, GSL_COMPLEX_ZERO, y);
    
    gsl_matrix_complex_free(ZI);
}


void initializeGSL()
{
  //instantiate the vectors and matrices
    int n = myPoints.size();
    gsl_vector_complex * y = gsl_vector_complex_alloc(n); //parameter of size
    gsl_vector_complex * c = gsl_vector_complex_alloc(n); //parameter of size
    gsl_matrix_complex * Z = gsl_matrix_complex_alloc(n, n); //parameter of size
    gsl_matrix_complex * G = gsl_matrix_complex_alloc(n, n); //parameter of size
    

    //load the y vector with the y values
    gsl_complex tempY;
    for(int i = 0; i < n; i++)
    {  
      tempY = gsl_complex_rect(myPoints[i].y, 0);
      gsl_vector_complex_set(y, i, tempY);
    }
    
    //load the Z matrix using equation 2 from project prompt

    double myNum = -1*(2*M_PI)/n;
    gsl_complex tempZ;
    double otherTemp = 0;
    for(int j = 0; j<n; j++)
    {
      for(int k = 0; k<n; k++)
      {
        otherTemp = myNum*j*k;
        tempZ = gsl_complex_rect(cos(otherTemp), sin(otherTemp));
        tempZ = gsl_complex_div_real(tempZ, sqrt(n));
        gsl_matrix_complex_set(Z, j, k, tempZ); 
      }
    }

    //find the coefficients for vector c, by multiplying Z*y
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, Z, y, GSL_COMPLEX_ZERO, c);
    
    //load the G matrix using equation 4 from project prompt

    double tempG = 0;
    gsl_complex compG;
    for(int i=0; i<n; i++)
     {
       for(int j=0; j<n; j++)
       {
         tempG = (-4*log(2)*i*j)/(pow(n, 1.5));
         if(i==j)
         {
           tempG = exp(tempG);
           compG = gsl_complex_rect(tempG, 0);
         }
         else{
           compG = GSL_COMPLEX_ZERO;
         }
         gsl_matrix_complex_set(G, i, j, compG);
       }
     } 

    //multiply G by c to get new coefficients
    gsl_vector_complex * GC = gsl_vector_complex_alloc(n);
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, G, c, GSL_COMPLEX_ZERO, GC);
    
    switch(myAnalysis.sizeFilter)
    {
      //apply inverse solver
      case 0:
        inverseSolver(Z, GC, y, n);
        break;
      case 1:
        directSolver(Z, GC, y, n);
        break;
      case 2:
        iterativeSolver();
        break;
    }
    
    //load the y points in the gsl vector back into the regular vector
    //when filtering is complete
    /*for(int i = 0; i < n; i++)
    {
      myPoints[i] = gsl_vector_get(y, i);
    }*/
    
    //free the vectors and matrices after being done with them
    gsl_vector_complex_free(y);
    gsl_vector_complex_free(GC);
    gsl_vector_complex_free(c);
    gsl_matrix_complex_free(G);
    gsl_matrix_complex_free(Z);
    
}


#endif
#ifndef DFT_CPP
#define DFT_CPP

#include <iostream>
#include <vector>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>
#include "prototypes.h"
#include "globals.h"

using namespace std;

void initializeGSL()
{
  //THESE MAY ALL NEED TO BE COMPLEX??
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
    double myNum = (2*M_PI)/n;
    gsl_complex myComp = gsl_complex_rect(cos(myNum), -sin(myNum));
    gsl_complex tempZ;
    for(int j = 0; j<n; j++)
    {
      for(int k = 0; k<n; k++)
      {
        tempZ = gsl_complex_pow_real(myComp, (j*k));
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
         tempG = (-4*log(2)*i*j)/(pow(n, (3/2)));
         tempG = exp(tempG);
         if(i==j)
         {
           tempG = tempG*1;
         }
         else{
           tempG = tempG*0;
         }
         compG = gsl_complex_rect(tempG, 0);
         gsl_matrix_complex_set(G, i, j, compG);
       }
     } 
    
    //multiply G by c to get new coefficients
    gsl_blas_zgemv(CblasNoTrans, GSL_COMPLEX_ONE, G, c, GSL_COMPLEX_ZERO, c);
    
    //find inverse of Z by replacing each element by its complex conjugate
    
    /* EQUATION 5 RETRIEVAL*/
    
    //compute Z inverse by finding the complex conjugate of all the entries in Z
    /* gsl_matrix_complex * ZI = gsl_matrix_complex_alloc(n, n);
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
    
    gsl_matrix_complex_free(ZI);*/
    
    /*DIRECT SOLVER--LU Decomposition*/
    int s; 
    gsl_permutation * p = gsl_permutation_alloc(n);
    gsl_linalg_complex_LU_decomp(Z, p, &s);
    gsl_linalg_complex_LU_solve(Z, p, c, y);
    
    /*ITERATIVE SOLVER--JACOBI*/

    
    //load the y points in the gsl vector back into the regular vector
    //when filtering is complete
    /*for(int i = 0; i < n; i++)
    {
      myPoints[i] = gsl_vector_get(y, i);
    }*/
    
    //free the vectors and matrices after being done with them
    gsl_vector_complex_free(y);
    gsl_vector_complex_free(c);
    gsl_matrix_complex_free(G);
    gsl_matrix_complex_free(Z);
    
}

#endif
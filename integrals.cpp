#ifndef INTEGRALS_CPP
#define INTEGRALS_CPP

//this file will eventually hold all of the integral methods
//composite simpsons, romberg, adaptive quadrature, guassian??
//TMS peak not included in rest of peaks for area
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "globals.h"
using namespace std;

#define PI 3.14159265

//return the cubic spline function
double myFunc(double x, int i)
{
    return (spline.a[i] + spline.b[i]*x + spline.c[i]*pow(x,2) + spline.d[i]*pow(x,3));
}

//composite simpson's method taken from book algorithm 4.1
double simpson(double a, double b, int n, int x)
{
    vector<double> XI;
    double X = 0;
    double result = 0;
    //STEP 1
    double h = (b-a)/n;

    //STEP 2
    XI.push_back(myFunc(a, x) + myFunc(b, x)); //XI0 
    XI.push_back(0); //XI1
    XI.push_back(0); //XI2
    //STEP 3
    for(int i=1; i<n; i++)
    {
        //STEP 4
        X = a + (i*h);

        //STEP 5
        if(i%2 == 0)
        {
            XI[2] = XI[2] + myFunc(X, x);
        }
        else{
            XI[1] = XI[1] + myFunc(X, x);
        }
    }

    //STEP 6
    result = h*((XI[0]) + (2*XI[2]) + (4*XI[1]))/3;

    //STEP 7
    return result;
}

double summation(int i, double a, double h, int x)
{
    double result = 0;
    for(int k=1; k<pow(2,i-2)+1; k++)
    {
        result = result + (myFunc(a+(k-0.5)*h, x));
    }
    
    return result;
}


//Romberg Integration method from algorithm 4.2 in the book
double romberg(double a, double b, int n, int x)
{
    //vector with n rows and n columns
    vector<vector<double> > R(3, vector<double>(n+1)); //vector with n rows and n columns
    //STEP 1
    double h = b-a;
    R[1][1] = (h/2)*(myFunc(a, x) + myFunc(b, x));


    //STEP 3
    for(int i=2; i<n+1;i++)
    {
        //STEP 4
        R[2][1] = (0.5)*(R[1][1] + (h*summation(i, a, h, x)));

        //STEP 5
        for(int j=2; j<i+1; j++)
        {
            R[2][j] = R[2][j-1]+ ((R[2][j-1]-R[1][j-1])/(pow(4, j-1)-1));
        }

        //STEP 7
        h=h/2;

        //STEP 8
        for(int j=1; j<i+1; j++)
        {
            R[1][j] = R[2][j];
        }
    }
    return R[1][n-1];

}

double adQuad(double a, double b, double tol, int n, int x)
{
    vector<double> TOL(n+1, 0), A(n+1, 0), H(n+1, 0), FA(n+1, 0), FC(n+1, 0), FB(n+1, 0), S(n+1, 0), L(n+1, 0), v(9,0);
    double FD = 0;
    double FE = 0;
    double S1 = 0;
    double S2 = 0;
    //STEP 1
    double APP = 0;
    int i = 1;
    TOL[i] = 10*tol;
    A[i] = a;
    H[i] = (b-a)/2;
    FA[i] = myFunc(a, x);
    FC[i] = myFunc(a+H[i], x);
    FB[i] = myFunc(b, x);
    S[i] = H[i]*(FA[i]+4*FC[i]+FB[i])/3;
    L[i]=1;
   
    //STEP 2
   while(i>0)
   {
       //STEP 3
       FD = myFunc(A[i]+(H[i]/2), x);
       FE = myFunc(A[i] + (3*H[i]/2), x);
       S1 = H[i]*(FA[i]+4*FD + FC[i])/6;
       S2 = H[i]*(FC[i]+4*FE + FB[i])/6;
       v[1] = A[i];
       v[2] = FA[i];
       v[3] = FC[i];
       v[4] = FB[i];
       v[5] = H[i];
       v[6] = TOL[i];
       v[7] = S[i];
       v[8] = L[i];

       //STEP 4
       i = i-1;

       //STEP 5
       if(abs(S1+S2-v[7])<v[6])
       {
           APP = APP + (S1+S2);
       }
       else{
           if(v[8]>=n)
           {
               break;
           }
           else{
               i = i + 1;
               A[i] = v[1]+v[5];
               FA[i] = v[3];
               FC[i] = FE;
               FB[i] = v[4];
               H[i] = v[5]/2;
               TOL[i] = v[6]/2;
               S[i] = S2;
               L[i] = v[8]+1;

                i = i+1;
                A[i] = v[1];
                FA[i] = v[2];
                FC[i] = FD;
                FB[i] = v[3];
                H[i] = H[i-1];
                TOL[i] = TOL[i-1];
                S[i] = S1;
                L[i] = L[i-1];
           }
       }
   }

   //STEP 6

   return APP;

}

 

double gaussLeg(double a, double b, int x)
{
    vector<double> r {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,  0.1834346424956498, 0.5255324099163290,  0.7966664774136267, 9.602898564975363};
    vector<double> c {0.1012285362903706, 0.2223810344533744, 0.3137066458778874, 0.3626837833783621, 0.3626837833783621, 0.3137066458778874, 0.2223810344533744, 0.1012285362903706};
    double result = 0;
    double funcMath = 0;
    for(int i=0; i<r.size();i++)
    {
        funcMath = (((b-a)*r[i])+(b+a))/2;
        result = result + c[i]*myFunc(funcMath, x)*((b-a)/2.0);
    }
    return result;
    
}

double performIntegration(int key)
{


  double integral = 0;
  double minIntegral = 100000; 
  double position = 0;
  

  switch(key){ 
  
    //Adaptive Quad
    case 0:
    //MAY NEED TO CHANGE THIS TO GO FROM END OF ROOTS TO BEGINNING
    //iterate through the roots, add the integral to the areas vector, find the minimum integral for the hydrogen
       for(int i = roots.size()-1; i > 0; i=i-2)
      {
          position = posRoots[i]/10; 
          integral = adQuad(roots[i].x-myPoints[position].x, splinePoints[posRoots[i]-1].x-myPoints[position].x, myAnalysis.tolerance, 10, position);
          for (int j=posRoots[i]-1;j>posRoots[i-1]+1;j--)
          {
            position = j/10;
            integral = integral + adQuad(splinePoints[j].x-myPoints[position].x,splinePoints[j-1].x-myPoints[position].x, myAnalysis.tolerance, 10, position);  
          }
          position = (posRoots[i-1]+1)/10;
          integral = integral + adQuad(splinePoints[posRoots[i-1]+1].x-myPoints[position].x,roots[i-1].x-myPoints[position].x, myAnalysis.tolerance, 10, position);
          areas.push_back(integral);
          if(integral<minIntegral)
          {
            minIntegral = integral;
          }
      }
        break;
        
    //Romberg Method
    case 1:
       for(int i = roots.size()-1; i > 0; i=i-2)
      {
          position = posRoots[i]/10; 
          integral = romberg(roots[i].x-myPoints[position].x, splinePoints[posRoots[i]-1].x-myPoints[position].x, 5, position);
          for (int j=posRoots[i]-1;j>posRoots[i-1]+1;j--)
          {
            position = j/10;
            integral = integral + romberg(splinePoints[j].x-myPoints[position].x,splinePoints[j-1].x-myPoints[position].x, 5, position);  
          }
          position = (posRoots[i-1]+1)/10;
          integral = integral + romberg(splinePoints[posRoots[i-1]+1].x-myPoints[position].x,roots[i-1].x-myPoints[position].x, 5, position);
          areas.push_back(integral);
          if(integral<minIntegral)
          {
            minIntegral = integral;
          }
      }
      break;
      
    //Simpson's Method
    case 2:
       for(int i = roots.size()-1; i > 0; i=i-2)
      {
          position = posRoots[i]/10; 
          integral = simpson(roots[i].x-myPoints[position].x, splinePoints[posRoots[i]-1].x-myPoints[position].x, 10, position);
          for (int j=posRoots[i]-1;j>posRoots[i-1]+1;j--)
          {
            position = j/10;
            integral = integral + simpson(splinePoints[j].x-myPoints[position].x,splinePoints[j-1].x-myPoints[position].x, 10, position);  
          }
          position = (posRoots[i-1]+1)/10;
          integral = integral + simpson(splinePoints[posRoots[i-1]+1].x-myPoints[position].x,roots[i-1].x-myPoints[position].x, 10, position);
          areas.push_back(integral);
          if(integral<minIntegral)
          {
            minIntegral = integral;
          }
      }
      break;
      
    //Gaussian Quadrature
    case 3:
     for(int i = roots.size()-1; i > 0; i=i-2)
      {
          position = posRoots[i]/10; 
          integral = gaussLeg(roots[i].x-myPoints[position].x, splinePoints[posRoots[i]-1].x-myPoints[position].x, position);
          for (int j=posRoots[i]-1;j>posRoots[i-1]+1;j--)
          {
            position = j/10;
            integral = integral + gaussLeg(splinePoints[j].x-myPoints[position].x,splinePoints[j-1].x-myPoints[position].x, position);  
          }
          position = (posRoots[i-1]+1)/10;
          integral = integral + gaussLeg(splinePoints[posRoots[i-1]+1].x-myPoints[position].x,roots[i-1].x-myPoints[position].x, position);
          areas.push_back(integral);
          if(integral<minIntegral)
          {
            minIntegral = integral;
          }
      }
      break;
  }
  return minIntegral;
}
#endif

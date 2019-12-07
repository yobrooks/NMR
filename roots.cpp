#ifndef ROOTS_CPP
#define ROOTS_CPP
//this file will hold all the functions dealing with the roots
#include <iostream>
#include <cmath>
#include <vector>
#include "structs.h"
#include "prototypes.h"
#include "globals.h"
using namespace std;

double bisection(double a, double b, int N, int i){
   double FA = myFunc(a, i);
  int k = 1;
   double p = 0; double FP = 0;
   while(k<=N)
   {
       p = a+(b-a)/2;
       FP = myFunc(p, i);

       if(FP == 0 || (b-a)/2<myAnalysis.tolerance)
       {
           break;
       }
       else{
           k = k+1;
           if((FA-myAnalysis.baseline)*(FP-myAnalysis.baseline)>0)
           {
               a = p;
               FA = FP;
           }
           else{
               b = p;
           }
       }
   }
   return p;
}

void findIntersect()
{
    double A = 0, B = 0, rootX = 0; Point tempRoot;
    //find the intersection between the points and the baseline
    for(int i = 0; i<splinePoints.size()-1; i++)
    {
        A = splinePoints[i].y-myAnalysis.baseline;
        B = splinePoints[i+1].y-myAnalysis.baseline;
        if(A*B<0)
        {
            
            rootX = bisection(splinePoints[i].x, splinePoints[i+1].x, 10, (i/10)); //use the bisection method to find the root
            tempRoot.x = rootX; //also just evaluate it at myFunc
            tempRoot.y = myAnalysis.baseline; //if a root is found then the y value should be at the baseline //also just evaluate it at myFunc
            roots.push_back(tempRoot);
            posRoots.push_back(i);
            
        }
    }
}  

double findTMS()
{
    double rightMost = -500000;

    //iterate through the x values to find the right most point which will be the TMS peak
    for(int i=0; i<myPoints.size(); i++)
    {
        if(myPoints[i].y>myAnalysis.baseline && myPoints[i].x>rightMost)
        {
            rightMost = myPoints[i].x;
        }
    }

    //iterate through all of the points and shift them so they are all on the 
    //left side of the TMS peak; should be still same interval away as before
    //will end up making TMS peak centered on 0 and everything else less than 0.
    for(int i=0; i<myPoints.size();i++)
    {
        myPoints[i].x = myPoints[i].x - rightMost;
    }
    return rightMost;
}

#endif

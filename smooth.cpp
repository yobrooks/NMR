#ifndef SMOOTH_CPP
#define SMOOTH_CPP
//this will hold functions for filtering and cubic spline
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include "globals.h"
#include "prototypes.h"
#include "structs.h"
using namespace std;


void prepVec(int n)
{
    //wrap the data at the beginning and end so the number of points doesn't change
    //add the beginning to the end
    for(int i=0; i<n; i++)
    {
        yFilter.push_back(yFilter[i]);
    }

    //add end to the beginning
    for(int i=0; i<n; i++)
    {
        yFilter.insert(yFilter.begin(), yFilter[yFilter.size()-(n-1)-i]);
    }
}

void boxcarFilter()
{
    
    int half = myAnalysis.sizeFilter/2; //needed to know how many points go before and after the one you're averaging for
    double sum = 0;
    
    //add y points to new vector used for filtering the y values
    for(int i = 0; i< myPoints.size(); i++)
    {
        yFilter.push_back(myPoints[i].y);
    }

     prepVec(half);


    //apply filtering starting at the middle of the filter size
    for(int i=half; i<yFilter.size()-half;i++)
    {
        for(int j = i-half; j<half+1+i; j++)
        {
            sum = sum + yFilter[j];
        }
        
        sum = sum/myAnalysis.sizeFilter;
        yFilter[i] = sum;
        sum = 0;
    }
    
    int count = 0;
     for(int i=half; i<yFilter.size()-half; i++)
    {
      myPoints[count].y = yFilter[i];
      count++;
    }
    
    yFilter.clear();

}

void sgFilter(int c, vector<int> typeFilter){
    //add y points to new vector used for filtering the y values
    for(int i = 0; i< myPoints.size(); i++)
    {
        yFilter.push_back(myPoints[i].y);
    }

    int half = myAnalysis.sizeFilter/2; //needed to know how many points go before and after the one you're averaging for
    double sum = 0;


    prepVec(half);
    vector<double> oldY = yFilter;
    for(int i=half; i<yFilter.size()-half; i++)
    {
        for(int j = i-half; j<i+half+1; j++)
        {
            sum = sum + (oldY[j]*typeFilter[j-i+half]);
        }
        yFilter[i] = sum/c;
        sum = 0;
    }
    int count  = 0;
    //once done reload points back into myPoints vector
    for(int i=half; i<yFilter.size()-half; i++)
    {
      myPoints[count].y = yFilter[i];
      count++;
    }
    
    yFilter.clear();
}

//function to make cubic spline
void cubicSpline(int n){
    spline.a.resize(n,0);
    spline.b.resize(n,0);
    spline.c.resize(n+1,0);
    spline.d.resize(n,0);

    //STEP 1
    vector<double> h;
    h.resize(n,0);
    for(int i=0; i<=n-1; i++)
    {
        h[i] = myPoints[i+1].x - myPoints[i].x;
        spline.a[i] = myPoints[i].y; //fill up the a coefficients vector 
    }

    //STEP 2
    vector<double> alpha;
    alpha.resize(n, 0);
    for(int i = 1; i<=n-1; i++)
    {
        alpha[i] = ((3/h[i])*(spline.a[i+1]-spline.a[i]))-((3/h[i-1])*(spline.a[i]-spline.a[i-1]));
    }

    //STEP 3
    vector<double>l, mu, z;
    l.resize(n+1,0); mu.resize(n+1,0); z.resize(n+1,0);
    l[0] = 1; mu[0] = 0; z[0] = 0;

    //STEP 4
    for(int i=1; i<=n-1; i++)
    {
        l[i] = 2*(myPoints[i+1].x-myPoints[i-1].x)-(h[i-1]*mu[i-1]);
        mu[i] = h[i]/l[i];
        z[i] = (alpha[i] - (h[i-1]*z[i-1]))/l[i];
    }

    //STEP 5
    l[n] = 1;
    z[n] = 0;
    spline.c[n] = 0;

    //STEP 6
     int j = 0;
    for(int i = 0; i<=n-1; i++)
    {
        j=n-1-i;
        spline.c[j] = z[j] - mu[j]*spline.c[j+1];
        spline.b[j] = (spline.a[j+1]-spline.a[j])/h[j]-h[j]*(spline.c[j+1]+2.0*spline.c[j])/3.0; 
        spline.d[j] = (spline.c[j+1]-spline.c[j])/(3.0*h[j]);
    }


    //after finding the coefficients for the spline
    //find the new points from the piecewise functions
    double numPoints = 10; //number of points between each point created by the spline
    Point tempPoint; double step = 0;
    for(int i = 0; i<spline.a.size(); i++)
    {
      step = (myPoints[i+1].x - myPoints[i].x)/numPoints;
      for(int j = 0; j <numPoints; j++)
      {
        tempPoint.x = myPoints[i].x + (j+1)*step;
        tempPoint.y = myFunc((j+1)*step, i);
        splinePoints.push_back(tempPoint);
      }
    }
    
}
#endif
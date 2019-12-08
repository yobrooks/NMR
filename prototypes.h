#ifndef PROTOTYPES_H
#define PROTOTYPES_H

#include <iostream>
#include <vector>
#include "structs.h"
using namespace std; 

double myFunc(double x, int i);
double bisection(double a, double b, int N, int i);
void findIntersect();
double findTMS();
void prepVec(int n);
void boxcarFilter();
void sgFilter(int c, vector<int> typeFilter);
void cubicSpline(int n);
double simpson(double a, double b, int n, int x);
double summation(int i, double a, double h, int x);
double romberg(double a, double b, int n, int x);
double adQuad(double a, double b, double tol, int n, int x);
double guassLeg(double a, double b, int x);
void readInNMR(string filename);
void readInData();
int initializeGSL();
double performIntegration(int key); //returns the minimum integral
#endif

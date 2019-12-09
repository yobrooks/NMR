#ifndef GLOBALS_H
#define GLOBALS_H

#include "structs.h"
#include <vector>
using namespace std;

extern Analysis myAnalysis;
extern Coefficients spline;
extern vector <Point> myPoints; //vector to hold the points
extern vector<double> yFilter; //vector that will be changing throughout filtering
extern vector<Point> splinePoints; //points after making the cubic spline
extern vector<Point> roots; //vector to hold the roots
extern vector<int> posRoots; //vector to hold what position in the splinePoints vector the roots occur at
extern vector<double> areas;
#endif
#ifndef STRUCTS_H
#define STRUCTS_H
#include <vector>
#include <iostream>
#include <string>
using namespace std;
struct Analysis{
    string inputFile;       //name of input file
    double baseline;        //baseline adjustment
    double tolerance;       //tolerance for numerical integration
    int filter;             //type of filter
    int sizeFilter;         //size of boxcar
    int numPasses;          //number of passes for filter
    int typeIntegration;    //type of numerical integration
    string outputFile;      //name of output file
};

struct Point{
    double x, y;
};

struct Coefficients{
    vector <double> a, b, c, d;
};
#endif
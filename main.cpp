//this file will run the NMR program
#ifndef MAIN_CPP
#define MAIN_CPP

#include <iostream>
#include <fstream>
#include <string>
#include "prototypes.h"
#include "structs.h"
#include "globals.h"
#include "constants.h"
using namespace std;


//coefficients used for the SG filter
vector<int> SG5 = {-3, 12, 17, 12, -3};
vector<int> SG11 = {-36, 9, 44, 69, 84, 89, 84, 69, 44, 9, -36};
vector<int> SG17 = {-21, -6, 7, 18, 27, 34, 39, 42, 43, 42, 39, 34, 27, 18, 7, -6, -21};


int main(int argc,char* argv[])
{

    string fileName;
    ofstream outputFile; 
    
    //read in the input file (nmr.in)
    readInNMR("nmr.in");
    readInData();
    
    cout << "Original Size of vector: " << myPoints.size() << endl;
    
    //open up an output file to output values to for the summary
    //output the beginning few lines of the summary
    outputFile.open(myAnalysis.outputFile);
    outputFile << "\t\t\t-=> NMR ANALYSIS <=-" << endl;
    outputFile << "Program Options" << endl;
    outputFile << "+++++++++++++++++++++++++++++++" << endl;
    outputFile << "Baseline Adjustment : " << myAnalysis.baseline << endl;
    outputFile << "Tolerance           : " << myAnalysis.tolerance << endl;
    
    double peak = findTMS();
    cout << "Peak: " << peak << endl;
    
    //switch statement to turn on the correct filtering technique.  
    switch(myAnalysis.filter)
    {
      //if 0 is inputted, then no filtering occurs
      case 0:
        outputFile << "Filtering has been turned off." << endl;
        break;
      case 1:
      //checking if filter size is odd 
        if(myAnalysis.sizeFilter%2==0)
        {
          outputFile << "Invalid Filter Size. Please enter an odd number" << endl;
        }
        else{
          outputFile << "Boxcar Filtering " << endl;
          outputFile << "Boxcar Size (Cyclic): " << myAnalysis.sizeFilter << endl;
          outputFile << "Boxcar Passes       : " << myAnalysis.numPasses << endl;
          
          //apply the boxcar filter using yFilter vector for the correct number of passes
          for(int i = 0; i <myAnalysis.numPasses; i++)
          {
            boxcarFilter();
          }
        }
        break;
        
      /* Apply the SG filter*/
        case 2: 
          
        outputFile << "Savitzy-Golay Filtering" << endl;
        outputFile << "Savitzy-Golay Size: " << myAnalysis.sizeFilter<< endl;
        outputFile << "Savitzy-Golay Passes: " << myAnalysis.numPasses << endl;
        
        //if the SG filter size is 5
        if(myAnalysis.sizeFilter == 5){
          for(int i = 0; i<myAnalysis.numPasses; i++){
            sgFilter(35, SG5);
          }
        }
        
        //if the filter size is 11
        else if(myAnalysis.sizeFilter == 11){
          for(int i = 0; i<myAnalysis.numPasses; i++){
           sgFilter(429, SG11);
           }
        }
        
        //if the filter size is 17
        else if(myAnalysis.sizeFilter == 17){
          for(int i=0; i<myAnalysis.numPasses; i++){
            sgFilter(323, SG17);
          }
        }
        else
          outputFile<< "Invalid filter size for Savitzy-Golay Filtering entered. Valid filter sizes are 5, 11, and 17." << endl;
        break;
      case 3:
        outputFile << "Discrete Fourier Transform Filter" << endl;
        initializeGSL();
        break;
    }
    
    /*
    outputFile << "\nIntegration Method" << endl;
    outputFile << "+++++++++++++++++++++++++++++++" << endl;
    //switch statement to print the appropriate integration technique
    switch(myAnalysis.typeIntegration)
    {
      case 0:
        outputFile << "Adaptive Quadrature" << endl;
        break;
      case 1: 
        outputFile << "Romberg Method" <<endl;
        break;
      case 2: 
        outputFile << "Simpson's Method" << endl;
        break;
      case 3:
        outputFile << "Gaussian-Legendre Quadrature" << endl; 
        break;
    } 
    
    outputFile << "\nPlot File Data" << endl;
    outputFile << "++++++++++++++++++++++++++++++++" << endl;
    outputFile << "File : " << myAnalysis.inputFile << endl;
    outputFile << "Plot Shifted " << peak << " ppm for this TMS calibration\n\n" << endl;
    
    outputFile << "Peak" << "Begin" << "End"<< "Location" << "Area" << "Hydrogens" << endl;
    
    //perform cubic spline
    cubicSpline(myPoints.size());
   
    //find the roots
    findIntersect();
    cout << "found intersect" << endl;
    double minInt = 0;
    minInt = performIntegration(myAnalysis.typeIntegration); //returns the smallest area
    cout << "found integration" << endl;
    cout << "AREAS IN MAIN: " << areas.size() << endl;
    for(int i = 0; i < (roots.size()/2); i=i+2)
    {
      cout << "for loop " << endl;
      cout << roots.size() << endl;
      outputFile << i+1 << roots[i*2].x << roots[i*2+1].x << roots[i*2+1].x-(roots[i*2].x/2) << areas[i] << int(areas[i]/minInt) << endl;
    }*/
    
    outputFile.close();
}
#endif

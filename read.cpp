//this file holds all of the data input functions
#ifndef READ_CPP
#define READ_CPP

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "structs.h"
#include "globals.h"

using namespace std; 

/*function to read in the nmr.in file and parse through it to get the correct values 
  to complete NMR. All of the information from here will be put into an Analysis struct*/
void readInNMR(string filename)
{
      string line;
      vector<string> values;
      ifstream myFile(filename);
    
    if(!myFile) //if file cannot be opened
    {
        cout << "Unable to open file." << endl;
        //exit(1);
    }
    else{ //if file can be opened
        string val;
        while(getline(myFile, line)) //get each line
        {
            val=line.substr(0, line.find("#")); //get value of string until # symbol
            val=val.substr(0,line.find(" ")); //cut off extra whitespace
            values.push_back(val); //put into vector which will be used to assign to 
        }
        //close the file when done
        myFile.close();
    }

    //assign values into Analysis struct
    myAnalysis.inputFile=values[0];
    myAnalysis.baseline=stod(values[1]);
    myAnalysis.tolerance=stod(values[2]);
    myAnalysis.filter=stoi(values[3]);
    myAnalysis.sizeFilter=stoi(values[4]);
    myAnalysis.numPasses=stoi(values[5]);
    myAnalysis.typeIntegration=stoi(values[6]);
    myAnalysis.outputFile=values[7];
}

/*function to read in the data points from the file specified by nmr.in. The data points will be put into a temporary Point and then pushed into the vector holding the original points*/
void readInData()
{
    string line;
    ifstream input(myAnalysis.inputFile);
    if(!input) //if file cannot be opened
    {
        cout << "Unable to open file 2." << endl;
      //  exit(1);
    }
    else{ //if file can be opened
        Point temp;
        while(getline(input, line)) //get each line
        {
            //assign values to temporary point and add it to points vector
            temp.x=stod(line.substr(0, line.find(" ")));
            temp.y=stod(line.substr(line.find(" ")+1));
            myPoints.push_back(temp);
        }
        //close the file when done
        input.close();
    }
}

#endif

#include <iostream>
#include <vector>
#include <cmath>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "prototypes.h"
#include "globals.h"

using namespace std;

gsl_vector * y; //MAKE THIS IN GLOBAL VARIABLES
gsl_matrix * Z;
gsl_vector * c;
gsl_matrix * G;

void initializeGSL()
{
    c = gsl_vector_alloc(); //parameter of size
    

    //WHEN LOADING VARIABLES ALSO ALLOW THIS TO BE 
    for(int i = 0;)   
}


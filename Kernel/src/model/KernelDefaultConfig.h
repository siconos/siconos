#ifndef __KERNELDEFAULTCONFIG__
#define __KERNELDEFAULTCONFIG__


#include "NewSiconosVector.h"


/* A matrix is saved in the XML output file if his size is not higher then MatrixMaxSize */
/*const*/ int   MATRIX_MAX_SIZE = 10;
/*const*/
int   VECTOR_MAX_SIZE = 10;
/*const*/
string  FILE_STORAGE = N_ASCII; // N_ASCII or N_BINARY

/*const*/
string  XML_SCHEMA = "/config/xmlschema/SiconosModelSchema-V1.2.xsd";

string  DefaultSolver = "default";
string  DefaultAlgoName = "default";
string  DefaultAlgoNormType = "default";
double  DefaultAlgoTolerance = -1.0;
int   DefaultAlgoMaxIter = -1;
double  DefaultAlgoSearchDirection = -1.0;


//#else
//
//extern int MATRIX_MAX_SIZE;
//extern int VECTOR_MAX_SIZE;
//extern string FILE_STORAGE;
//extern string XML_SCHEMA;
//
//extern string   DefaultSolver;
//extern string   DefaultAlgoName;
//extern string   DefaultAlgoNormType;
//extern double   DefaultAlgoTolerance;
//extern int    DefaultAlgoMaxIter;
//extern double   DefaultAlgoSearchDirection;


#endif


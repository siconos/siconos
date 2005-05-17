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
string  XML_SCHEMA = "/share/SICONOS/SiconosModelSchema-V1.2.xsd";

string  DefaultSolver = "default";
string  DefaultAlgoName = "default";
string  DefaultAlgoNormType = "default";
double  DefaultAlgoTolerance = -1.0;
int   DefaultAlgoMaxIter = -1;
double  DefaultAlgoSearchDirection = -1.0;

string DefaultComputeInput = "BasicPlugin:computeInput";
string DefaultComputeOutput = "BasicPlugin:computeOutput";

// To initialize pointers:
#ifndef NULL
extern const int NULL = 0;
#endif

#endif


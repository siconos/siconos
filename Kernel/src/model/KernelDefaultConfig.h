//$Id: KernelDefaultConfig.h,v 1.5 2005/02/10 10:35:18 jbarbier Exp $
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

//$Log: KernelDefaultConfig.h,v $
//Revision 1.5  2005/02/10 10:35:18  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.4  2005/01/25 10:33:15  jbarbier
//- modifications for test purpose
//
//Revision 1.3  2005/01/25 09:27:16  jbarbier
//- save of Solver tag in the OneStepNSProblem tag available when saving without XML input file and with partial XML input file
//
//Revision 1.2  2005/01/20 09:05:34  jbarbier
//- configuration file available and usable
//
//- save of vectors and matrices into external files (version 0.1)
//
//Revision 1.1  2005/01/18 16:46:03  jbarbier
//- file containing default parameters for the kernel added in src/model
//

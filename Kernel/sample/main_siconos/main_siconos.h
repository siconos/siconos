#ifndef MAIN_SICONOS_H
#define MAIN_SICONOS_H

#include <stdio.h>
#include <sys/time.h>
#include <iostream>
#include "Model.h"
#include "check.h"

#include "TimeStepping.h"
#include "EventDriven.h"

#include "DynamicalSystem.h"
#include "LagrangianDS.h"
#include "LinearSystemDS.h"
#include "LagrangianLinearTIDS.h"
#include "LagrangianNonLinearR.h"
#include "LCP.h"

#include "LagrangianDSIO.h"
#include "LinearTIEC.h"

#include <libxml/parser.h>
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "SiconosMatrix.h"
#include "SiconosDOMTreeTools.h"

//#include "SiconosNumerics.h"

using namespace std;

Model m;

extern "C"
{
  void cartouche();
  void essai_model();
  void essai_model_XML(char []);
  void essai_model2();
  void test_schema();
  void bench();
}

//int main(int argc, char* argv[]);



#endif //MAIN_SICONOS_H

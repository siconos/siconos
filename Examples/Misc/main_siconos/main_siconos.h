/* Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/
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

//Model m;

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

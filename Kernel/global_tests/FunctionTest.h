/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

/** \class FunctionTest
 *  \brief Tool to generate a function to be tested
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) february 28, 2006
 *
 */

#ifndef FUNCTIONTEST_H
#define FUNCTIONTEST_H

#include "SiconosException.h"
#include "RuntimeException.h"
#include "SiconosMatrix.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>
#include <set>
#include<string>

typedef bool (*pointerToTestFunction)();
typedef bool func();
class FunctionTest
{
private:

  //  SiconosMatrix * dataOut;
  /** name of the test \fn void vectorField (const double& time) */
  std::string name;

  /** true if test succeed */
  bool isOk;

  //std::string referenceFile;

  /** pointer to function that handle the simulation */
  pointerToTestFunction fp;

public:

  /** \fn void FunctionTest(const std::string&, func);
   * \brief constructor from function name
   * \param the name of the test
   * \param a function of type func
   */
  FunctionTest(const std::string&, func);

  /** \fn ~FunctionTest();
   * \brief Destructor
   */
  ~FunctionTest() {};

  /** \fn bool hasFailed()
   * \brief test if test failed or not
   */
  inline bool hasFailed() const
  {
    return !isOk;
  } ;

  /** \fn std::string getName() const {return name;} ;
   * \brief return test name
   */
  inline std::string getName() const
  {
    return name;
  } ;

  /** \fn void run()
   * \brief compute simulation
   */
  void  run();

  //  void  compare();
};

#endif

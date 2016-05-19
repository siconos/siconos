/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/** \class FunctionTest
 *  \brief Tool to generate a function to be tested
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) february 28, 2006
 *
 */

#ifndef FUNCTIONTEST_H
#define FUNCTIONTEST_H

#include "SiconosException.hpp"
#include "RuntimeException.hpp"
#include "SiconosMatrix.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include<string>

typedef bool (*pointerToTestFunction)();
typedef bool func();
class FunctionTest
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FunctionTest);


  //  SiconosMatrix * dataOut;
  /** name of the test function void vectorField (double time) */
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

  /** \fn std::string name() const {return name;} ;
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

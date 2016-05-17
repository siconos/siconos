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

/** \class GlobalTest
 *  \brief Tool to generate a list of function to be tested
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) february 28, 2006
 *
 */

#ifndef GLOBALTEST_H
#define GLOBALTEST_H

#include "SiconosException.hpp"
#include "RuntimeException.hpp"
#include"FunctionTest.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include<string>

class GlobalTest
{
private:

  /** list of tests to be run */
  std::set<FunctionTest*> testsList;

  /** output log file */
  std::ofstream out;

  /** number of tests that failed (ie isOk parameter = false) */
  unsigned int numberOfFailedTests;

public:

  /** \fn GlobalTest(const std::ofstream&)
   * \brief constructor from output file name
   * \param the name of the output file
   */
  GlobalTest(const std::string& = "none");

  /** \fn ~GlobalTest()
   * \brief Destructor
   */
  ~GlobalTest();

  /**\fn   inline unsigned int getNumberOfTests() const
   * \brief return the total number of tests
   */
  inline unsigned int getNumberOfTests() const
  {
    return testsList.size();
  }  ;

  /** \fn   inline unsigned int getNumberOfFailedTests() const
   * \brief return the number of tests that failed
   */
  inline unsigned int getNumberOfFailedTests() const
  {
    return numberOfFailedTests;
  }  ;

  /** \fn void addTest()
   * \brief add a new FunctionTest in the list
   */
  void addTest(const std::string&, func);

  /** \fn void run()
   * \brief compute all the tests
   */
  void  run();

  /** \fn void print()
   * \brief print results into a file
   */
  void print();
};

#endif

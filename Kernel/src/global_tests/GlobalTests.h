/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

/** \class GlobalTest
 *  \brief Tool to generate a list of function to be tested
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) february 28, 2006
 *
 */

#ifndef GLOBALTEST_H
#define GLOBALTEST_H

#include "SiconosException.h"
#include "RuntimeException.h"
#include"FunctionTest.h"
#include <fstream>
#include <iostream>
#include <math.h>
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

  /** \fn bool hasFailed()
   * \brief return the total number of tests
   */
  inline const unsigned int getNumberOfTests() const
  {
    return testsList.size();
  }  ;

  /** \fn bool hasFailed()
   * \brief return the number of tests that failed
   */
  inline const unsigned int getNumberOfFailedTests() const
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

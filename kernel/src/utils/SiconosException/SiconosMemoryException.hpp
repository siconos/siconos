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
/*! \file SiconosMemoryException.hpp
    \brief  SiconosMemoryException class

*/
#ifndef SICONOSMEMORYEXCEPTION_H
#define SICONOSMEMORYEXCEPTION_H

#include "SiconosException.hpp"

/** Exceptions for SiconosMemory
 *
 * \author SICONOS Development Team - copyright INRIA
 * \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 */
class SiconosMemoryException : public SiconosException
{
public:

  /**
   * constructor
   */
  SiconosMemoryException();

  /**
   * constructor
   * \param report the exception description
   */
  SiconosMemoryException(const std::string& report);

  /**
   * destructor
   */
  ~SiconosMemoryException();

  static void selfThrow() NO_RETURN;

  static void selfThrow(const std::string& report) NO_RETURN;

};

#endif //SICONOSMEMORYEXCEPTION_H

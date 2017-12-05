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
#ifndef SICONOSSHAREDLIBRARYEXCEPTION_H
#define SICONOSSHAREDLIBRARYEXCEPTION_H

#include "SiconosException.hpp"

/*! \file SiconosSharedLibraryException.hpp

*/

/** Exceptions for SiconosSharedLibrary
 *
 * \author SICONOS Development Team - copyright INRIA
 * \date (creation) 07/21/2006
 *  Matrices can be either block or Simple.
 *  See Derived classes for details.
 */
class SiconosSharedLibraryException : public SiconosException
{
public:

  /**
   * constructor
   */
  SiconosSharedLibraryException();

  /**
   * constructor
   * \param report the description of the exception
   */
  SiconosSharedLibraryException(const std::string& report);

  /**
   * destructor
   */
  ~SiconosSharedLibraryException();

  static void selfThrow() NO_RETURN;

  static void selfThrow(const std::string& report) NO_RETURN;

};

#endif //SICONOSSHAREDLIBRARYEXCEPTION_H

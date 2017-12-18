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


/*! \file SiconosMatrixException.hpp
 */
#ifndef __SiconosMatrixException__
#define __SiconosMatrixException__

#include "SiconosException.hpp"

/** Exception caused by a SiconosMatrix
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 05/25/2004
 *
 *
 *
 * SiconosMatrixException must be throws when an error is find in a SiconosMatrix
 * This exception can be caught by "catch(SiconosMatrixException)" or "catch(SiconosException)"
 *
 */
class SiconosMatrixException : public SiconosException
{
public:

  /** constructor
   */
  SiconosMatrixException();

  /** constructor with a report
   * \param report exception description
   */
  SiconosMatrixException(const std::string& report);

  /** destructor
   */
  ~SiconosMatrixException();

  /** static function which throw a SiconosMatrixException
   */
  static void selfThrow() NO_RETURN ;

  /** static function which throw a SiconosMatrixException with a report
   * \param report exception description
   * \exception SiconosMatrixException
   */
  static void selfThrow(const std::string& report) NO_RETURN;
};

#endif //__SiconosMatrixException__

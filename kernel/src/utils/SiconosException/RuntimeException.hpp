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
/*! \file RuntimeException.hpp
 */

#ifndef __RuntimeException__
#define __RuntimeException__

#include "SiconosException.hpp"

/** Runtime exceptions
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 05/25/2004
 *
 *
 * RuntimeException can be throws for example when a pointer is used but not allocated
 * This exception can be caught by "catch(RuntimeException)" or "catch(SiconosException)"
 *
 */
class RuntimeException: public SiconosException
{
public:

  /** constructor
   */
  RuntimeException();

  /** constructor with a report
   * \param report exception description
   */
  RuntimeException(const std::string& report);

  /** destructor
   */
  ~RuntimeException();

  /** static function which throw a RuntimeException
   *
   */
  static void selfThrow() NO_RETURN;

  /** static function which throw a RuntimeException with a report
   * \param report exception description
   *
   */
  static void selfThrow(const std::string& report) NO_RETURN;

};

#endif //__RuntimeException__

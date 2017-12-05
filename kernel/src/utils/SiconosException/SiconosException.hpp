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

/*! \file SiconosException.hpp
 */

#ifndef __SiconosException__
#define __SiconosException__

#include <string>
#include "SiconosSerialization.hpp"

#ifdef __clang_analyzer__
#define NO_RETURN  __attribute__((analyzer_noreturn))
#else
#define NO_RETURN
#endif

/** General Siconos Exception
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 05/25/2004
 *
 *
 * SiconosException should not be throws directly; prefer to use an inherit class
 * This exception can be caught by "catch(SiconosException)"
 *
 */
class SiconosException
{
public:
  /** constructor
   */
  SiconosException();

  /** constructor with a report
   * \param report exception description
   */
  SiconosException(const std::string& report);

  /** destructor
  */
  virtual ~SiconosException();

  /** return the report of the exception
   * \return std::string report : exception description
   */
  inline std::string report() const
  {
    return _reportMsg;
  } ;

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosException);

  /** report message which describe the exception */
  std::string _reportMsg;
};

#endif //__SiconosException__

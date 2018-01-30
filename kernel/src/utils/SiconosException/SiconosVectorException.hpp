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

/*! \file SiconosVectorException.hpp
 */

#ifndef __SiconosVectorException__
#define __SiconosVectorException__

#include "SiconosException.hpp"

/** Exception caused by a SiconosVector
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 05/25/2004
 *
 *
 *
 * SiconosVectorException must be throws when an error is find in a SiconosVector
 * This exception can be caught by "catch(SiconosVectorException)" or "catch(SiconosException)"
 *
 */
class SiconosVectorException : public SiconosException
{
public:

  /** constructor
   */
  SiconosVectorException();

  /** constructor with a report
   * \param report exception description
   */
  SiconosVectorException(const std::string& report);

  /** destructor
   */
  ~SiconosVectorException();

  /** static function which throw a SiconosVectorException
   */
  static void selfThrow() NO_RETURN;

  /** static function which throw a SiconosVectorException with a report
   * \param report exception description
   * \exception SiconosVectorException
   */
  static void selfThrow(const std::string& report) NO_RETURN;

};

#endif //__SiconosVectorException__

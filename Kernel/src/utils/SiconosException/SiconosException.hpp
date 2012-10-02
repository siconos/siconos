/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file SiconosException.hpp
 */

#ifndef __SiconosException__
#define __SiconosException__

#include <string>
#include "SiconosSerialization.hpp"

/** General Siconos Exception
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 05/25/2004
 *
 *
 * SiconosException should not be throws directly; prefer to use an inherit class
 * This exception can be catched by "catch(SiconosException)"
 *
 */
class SiconosException
{
public:
  /** constructor
   */
  SiconosException();

  /** constructor with a report
   * \param string report : exception description
   */
  SiconosException(const std::string&);

  /** destructor
  */
  virtual ~SiconosException();

  /** return the report of the exception
   * \return string report : exception description
   */
  inline std::string report() const
  {
    return reportMsg;
  } ;

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosException);

  /** report message which describe the exception */
  std::string reportMsg;
};

#endif //__SiconosException__

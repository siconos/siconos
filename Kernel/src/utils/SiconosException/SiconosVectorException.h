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

/*! \file SiconosVectorException.h
 */

#ifndef __SiconosVectorException__
#define __SiconosVectorException__

#include "SiconosException.h"

/** Exception caused by a SiconosVector
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date (Creation) 05/25/2004
 *
 *
 *
 * SiconosVectorException must be throws when an error is find in a SiconosVector
 * This exception can be catched by "catch(SiconosVectorException)" or "catch(SiconosException)"
 *
 */
class SiconosVectorException : public SiconosException
{
public:

  /** constructor
   */
  SiconosVectorException();

  /** constructor with a report
   * \param std::string report : exception description
   */
  SiconosVectorException(const std::string& report);

  /** destructor
   */
  ~SiconosVectorException();

  /** static function which throw a SiconosVectorException
   * \exception SiconosVectorException
   */
  static void selfThrow()  ;

  /** static function which throw a SiconosVectorException with a report
   * \param std::string report : exception description
   * \exception SiconosVectorException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__SiconosVectorException__

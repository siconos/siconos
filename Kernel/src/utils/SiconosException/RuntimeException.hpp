/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file RuntimeException.h
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
 * This exception can be catched by "catch(RuntimeException)" or "catch(SiconosException)"
 *
 */
class RuntimeException: public SiconosException
{
public:

  /** constructor
   */
  RuntimeException();

  /** constructor with a report
   * \param string report : exception description
   */
  RuntimeException(const std::string& report);

  /** destructor
   */
  ~RuntimeException();

  /** static function which throw a RuntimeException
   * \exception RuntimeException
   */
  static void selfThrow() ;

  /** static function which throw a RuntimeException with a report
   * \param string report : exception description
   * \exception RuntimeException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__RuntimeException__

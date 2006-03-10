/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
/** \class RuntimeException
*   \brief This class represent a runtime exeption causing by the plateforme
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.1.3.
*   \date (Creation) 05/25/2004
*
*
* RuntimeException can be throws for example when a pointer is used but not allocated
* This exception can be catched by "catch(RuntimeException)" or "catch(SiconosException)"
*
*/

#ifndef __RuntimeException__
#define __RuntimeException__

#include "SiconosException.h"

class RuntimeException: public SiconosException
{
public:

  /**
   * \fn RuntimeException()
   * \brief constructor
   */
  RuntimeException();

  /**
   * \fn RuntimeException(const string& report)
   * \brief constructor with a report
   * \param string report : exception description
   */
  RuntimeException(const std::string& report);

  /**
   * \fn ~RuntimeException()
   * \brief destructor
   */
  ~RuntimeException();

  /**
   * \fn static void selfThrow()
   * \brief static function which throw a RuntimeException
   * \exception RuntimeException
   */
  static void selfThrow() ;

  /**
   * \fn static void selfThrow(string report)
   * \brief static function which throw a RuntimeException with a report
   * \param string report : exception description
   * \exception RuntimeException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__RuntimeException__

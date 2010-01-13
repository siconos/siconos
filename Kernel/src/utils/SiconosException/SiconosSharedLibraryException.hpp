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
#ifndef SICONOSSHAREDLIBRARYEXCEPTION_H
#define SICONOSSHAREDLIBRARYEXCEPTION_H

#include "SiconosException.hpp"

/*! \file SiconosSharedLibraryException.h

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
   * @param string which describe the exception
   */
  SiconosSharedLibraryException(const std::string& report);

  /**
   * destructor
   */
  ~SiconosSharedLibraryException();

  static void selfThrow() ;

  static void selfThrow(const std::string& report) ;

};

#endif //SICONOSSHAREDLIBRARYEXCEPTION_H

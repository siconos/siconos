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

/*! \file ioObject.hpp
  Interface for input/output of SiconosAlgebra
*/

#ifndef __ioObject__
#define __ioObject__

#include <string>
#include "RuntimeException.hpp"

/** Interface for read/write objects from/to a file.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.1.0.
 *   \date (Creation) 07/21/2006
 *
 */
template <class T> class ioObject
{
protected:

  /** std::string fileName : the file to read
   */
  std::string FileName;

  /**std::string mode : ASCII or BINARY
   */
  std::string Mode;

  /** default constructor
   */
  ioObject(const std::string& m = "ascii"): FileName("NoName.dat"), Mode(m) {}

public :

  /** constructor
   *  \param string: input/output file name
   *  \param string: ascii or binary
   */
  ioObject(const std::string& file, const std::string& m): FileName(file), Mode(m) {}

  /** destructor
   */
  ~ioObject() {}

  /** read the matrix in the file "Filename" and write it into matrix A
   *  \param a SiconosMatrix
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool read(T&) const
  {
    RuntimeException::selfThrow("ioObject::read(T) - not implemented");
    return false;
  }

  /** write the matrix A in the file "Filename"
   *  \param a SiconosMatrix
   *  \param a string: type of output - See on top of file for details
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  bool write(const T&, const std::string& = "python") const
  {
    RuntimeException::selfThrow("ioObject::write(T) - not implemented");
    return false;
  }
};
#endif

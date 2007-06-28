/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

/*! \file ioMatrix.h
    \brief

*/

#ifndef __ioMatrix__
#define __ioMatrix__

#include "ioObject.h"

class SimpleMatrix;

/** Interface for read/write matrices from/to a file.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date (Creation) 07/21/2006
 *
 *
 *  Type of Output for write function:
 *    - "boost": boost way:
 *                    [row,col] ((a00,a01,...),(a10,...),...
 *    - "python"(default):
 *                    row col
 *                    a00 a01 a02 ...
 *                    a10 ...
 *    - "noDim":
 *                    a00 a01 a02 ...
 *                    a10 ...
 *
 * Reading input format is the one corresponding to "python".
 *
 */
class ioMatrix : public ioObject
{
private :

  /** default constructor with Mode = "ascii"
  */
  ioMatrix();

public :
  /**static SimpleMatrix *temporary
   * used for storing matrix written in a file when Mode = "binary"
   */
  static SimpleMatrix* temporary;

  /** static writeSimpleBinary
   * true if temporary is allocated, else false.
   */
  static bool writeSimpleBinary;

  /** constructor with FileName = file and Mode = mode
   *  \param 2 std::string
   */
  ioMatrix(const std::string&, const std::string&);

  /** destructor
   */
  ~ioMatrix(void);

  /** read the matrix in the file "Filename" and write it into matrix A
   *  \param a SiconosMatrix
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  const bool read(SiconosMatrix&) const;

  /** write the matrix A in the file "Filename"
   *  \param a SiconosMatrix
   *  \param a string: type of output - See on top of file for details
   *  \exception SiconosMatrixException
   *  \return true if no error
   */
  const bool write(const SiconosMatrix&, const std::string& = "python") const;

};
#endif

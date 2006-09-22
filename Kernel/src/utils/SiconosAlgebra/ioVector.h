/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

/** \class ioVector
*   \brief This class is an interface for read/write vectors in a file.
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.3.0.
*   \date (Creation) 07/21/2006
*
*/


#ifndef __ioVector__
#define __ioVector__

#include "ioObject.h"

class ioVector : public ioObject
{
private :

  /** \fn ioVector ()
  *  \brief default constructor with Mode = "ascii"
  */
  ioVector(void);

public :
  /** \var static MySimpleVector *temporary
  *  \brief used for storing vector read from a file when Mode = "binary"
  */
  static MySimpleVector* temporary;

  /** \var static bool writeSimpleBinary
  *  \brief true if temporary is allocated, else false.
  */
  static bool writeSimpleBinary;

  /** \fn ioVector (const std::string& file, const std::string& mode)
  *  \brief constructor with FileName = file and Mode = mode
  *  \param 2 std::string
  */
  ioVector(const std::string&, const std::string&);

  /** \fn ~ioVector ()
  *  \brief destructor
  */
  ~ioVector(void);


  /** \fn bool read(MySiconosVector &A)
  *  \brief read the vector in the file "Filename" and write it into vector A
  *  \param a MySiconosVector
  *  \exception SiconosVectorException
  *  \return true if no error
  */
  bool    read(MySiconosVector&)const;

  /** \fn bool write(MySiconosVector &A)
  *  \brief write the vector A in the file "Filename"
  *  \param a MySiconosVector
  *  \exception SiconosVectorException
  *  \return true if no error
  */
  bool    write(const MySiconosVector&);

  /** \fn bool rawWrite(MySiconosVector &A)
  *  \brief write the vector A in the file "Filename" without its dimensions
  *  \param a MySiconosVector
  *  \exception SiconosVectorException
  *  \return true if no error
  */
  bool    rawWrite(const MySiconosVector&);

  /** \fn bool read(MySiconosMatrix &A)
  *  \brief read the matrix in the file "Filename" and write it into matrix A, useless for ioVector
  *  \param a MySiconosMatrix
  *  \exception SiconosMatrixException
  *  \return true if no error
  */
  bool    read(MySiconosMatrix&)const;

  /** \fn bool write(MySiconosMatrix &A)
  *  \brief write the matrix A in the file "Filename", useless for ioVector
  *  \param a MySiconosMatrix
  *  \exception SiconosMatrixException
  *  \return true if no error
  */
  bool    write(const MySiconosMatrix&);

  /** \fn bool rawWrite(MySiconosMatrix &A)
  *  \brief write the matrix A in the file "Filename" without its dimensions, useless for ioVector
  *  \param a MySiconosMatrix
  *  \exception SiconosMatrixException
  *  \return true if no error
  */
  bool    rawWrite(const MySiconosMatrix&);


};
#endif

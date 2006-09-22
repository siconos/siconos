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

/** \class ioObject
*   \brief This class is an interface for read/write objects in a file.
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.3.0.
*   \date (Creation) 07/21/2006
*
*/

#ifndef __ioObject__
#define __ioObject__

#include <fstream>
#include <iostream>
#include <string>
#include <cassert>

#include "MyBlockMatrix.h"

class ioObject
{
public :

  /** \var std::string fileName : the file to read
  */
  std::string FileName;

  /** \var std::string mode : ASCII or BINARY
  */
  std::string Mode;

  /** \fn ~ioObject ()
  *  \brief destructor
  */
  virtual ~ioObject(void)        = 0;

  /** \fn bool read(MySiconosMatrix &A)
  *  \brief read the matrix in the file "Filename" and write it into matrix A
  *  \param a MySiconosMatrix
  *  \exception SiconosMatrixException
  *  \return true if no error
  */
  virtual bool    read(MySiconosMatrix&)const  = 0;

  /** \fn bool write(MySiconosMatrix &A)
  *  \brief write the matrix A in the file "Filename"
  *  \param a MySiconosMatrix
  *  \exception SiconosMatrixException
  *  \return true if no error
  */
  virtual bool    write(const MySiconosMatrix&) = 0;

  /** \fn bool rawWrite(MySiconosMatrix &A)
  *  \brief write the matrix A in the file "Filename" without its dimensions
  *  \param a MySiconosMatrix
  *  \exception SiconosMatrixException
  *  \return true if no error
  */
  virtual bool    rawWrite(const MySiconosMatrix&) = 0;

  /** \fn bool read(MySiconosVector &A)
  *  \brief read the vector in the file "Filename" and write it into vector A
  *  \param a MySiconosVector
  *  \exception SiconosVectorException
  *  \return true if no error
  */
  virtual bool    read(MySiconosVector&)const  = 0;

  /** \fn bool write(MySiconosVector &A)
  *  \brief write the vector A in the file "Filename"
  *  \param a MySiconosVector
  *  \exception SiconosVectorException
  *  \return true if no error
  */
  virtual bool    write(const MySiconosVector&) = 0;

  /** \fn bool rawWrite(MySiconosVector &A)
  *  \brief write the vector A in the file "Filename" without its dimensions
  *  \param a MySiconosVector
  *  \exception SiconosVectorException
  *  \return true if no error
  */
  virtual bool    rawWrite(const MySiconosVector&) = 0;

};
#endif

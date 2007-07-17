/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
/*! \file LinearDSIO.h

*/
#ifndef LINEARDSIO_H
#define LINEARDSIO_H

#include "DSInputOutput.h"
#include "LinearDSIOXML.h"

//! Linear DSInputOutput
/**         { y = A.x
 *         { R = B.lambda
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date Apr 27, 2004
 *
 *
 *
 */
class LinearDSIO : public DSInputOutput
{
public:

  /** Default constructor
   */
  LinearDSIO();

  /** constructor with XML object of the LinearDSIO
   *  \param LinearDSIOXML* : the XML object corresponding
   */
  LinearDSIO(DSInputOutputXML*);

  /** Destructor
   */
  ~LinearDSIO();

  /** copy the data of the DSInputOutput to the XML tree
   */
  void saveDSInputOutputToXML();

  /** allows to create the DSInputOutput with an xml file, or the needed data
   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
   *  \param int : the number of the DSInputOutput
   *  \param SiconosMatrix* : matrix A of the linear DSIO
   *  \param SiconosMatrix* : matrix B of the linear DSIO
   *  \exception RuntimeException
   */
  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
                           SiconosMatrix *A = NULL, SiconosMatrix *B = NULL);

protected:
  /** uses the DSInputOutputXML of the LinearDSIO to fill the fields of this DSInputOutput
   */
  void fillDSInputOutputWithDSInputOutputXML();

private:
  SiconosMatrix* A;
  SiconosMatrix* B;
};

#endif // LINEARDSIO_H

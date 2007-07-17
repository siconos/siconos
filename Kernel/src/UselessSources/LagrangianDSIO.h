/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file LagrangianDSIO.h

*/
#ifndef LAGRANGIANDSIO_H
#define LAGRANGIANDSIO_H

#include "DSInputOutput.h"
#include "LagrangianDSIOXML.h"
#include "check.h"

//! Lagrangian DSInputOutput
/**         { y = H(q, t)
 *         { R = G(lambda)
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date 17/01/2005
 *
 *
 */
class LagrangianDSIO : public DSInputOutput
{
public:

  /** default constructor
  */
  LagrangianDSIO();
  ~LagrangianDSIO();

  //   /** copy the data of the DSInputOutput to the XML tree
  //   */
  //  void saveDSInputOutputToXML();

  //   /** allows to create the DSInputOutput with an xml file, or the needed data
  //   *  \param DSInputOutputXML * : the XML object for this DSInputOutput
  //   *  \exception RuntimeException
  //   */
  //  void createDSInputOutput(DSInputOutputXML * dsioXML, int number = -1,
  //                string computeInput=DefaultComputeInput,
  //                string computeOutput=DefaultComputeOutput);

private:

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  //   /** to be defined
  //   */
  //  void (*computeJacobianPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* jacobPtr);
  //
  //  void (*computeHPtr)(int* sizeOfQ, double* qPtr, int* sizeOfY, double* yPtr);

  /** computes y
  *  \param double* xPtr : the pointer to the first element of the vector x
  *  \param double* time : the current time
  *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
  *  \param double* yPtr : the pointer to the first element of the vector y (in-out parameter)
  */
  void (*computeOutputPtr)(double* xPtr, double* time, double* lambdaPtr, double* yPtr);

  /** computes r
  *  \param double* xPtr : the pointer to the first element of the vector x
  *  \param double* time : the current time
  *  \param double* lambdaPtr : the pointer to the first element of the vector lambda
  *  \param double* rPtr : the pointer to the first element of the vector r (in-out parameter)
  */
  void (*computeInputPtr)(double* xPtr, double* time, double* lambdaPtr, double* rPtr);
};

#endif

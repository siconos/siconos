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
/*! \file LagrangianLinearDSIOXML.h

*/
#ifndef LAGRANGIANLINEARDSIOXML_H
#define LAGRANGIANLINEARDSIOXML_H

#include "LagrangianDSIOXML.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

const std::string LLDSIO_H = "H";
const std::string LLDSIO_B = "b";

//! XML management for LagrangianLinear DSInputOutput
/**  \author SICONOS Development Team - copyright INRIA
*   \version 2.1.1.
*   \date 05/25/2004
*
*
*
* LagrangianLinearDSIOXML allows to manage data of a LagrangianLinearDSIO DOM tree.
*/
class LagrangianLinearDSIOXML : public LagrangianDSIOXML
{
public:

  LagrangianLinearDSIOXML();
  virtual ~LagrangianLinearDSIOXML();

  /** Build a LagrangianLinearDSIOXML object from a DOM tree describing a DSIO with LagrangianLinear type
  *   \param LagrangianLinear : the LagrangianLinearDSIO DOM tree
  *   \exception XMLException : if a property of the LagrangianLinear DSIO lacks in the DOM tree
  */
  LagrangianLinearDSIOXML(xmlNode * dsioNode);

  //////////////////////////////////////////////

  /** Return the H of the LagrangianLinearDSIOXML
  *   \return The H SimpleMatrix of the LagrangianLinearDSIOXML
  */
  inline SimpleMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }


  /** Return b vector of the LagrangianLinearDSIOXML
  *   \return SimpleVector : b vector of the LagrangianLinearDSIOXML
  */
  inline SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->bNode);
  }

  /** Change the H matrix value (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for H matrix
  */
  void setH(SiconosMatrix *matrix);

  /** Change the b vector value (in xml file or external data file switch his origin position)
  *   \param SiconosVector vector : the new value for b vector
  */
  void setB(SiconosVector *vector);


private:
  //Nodes
  xmlNode * HNode;
  xmlNode * bNode;
};

#endif // LAGRANGIANLINEARDSIOXML_H

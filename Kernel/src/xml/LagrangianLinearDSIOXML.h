/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#ifndef LAGRANGIANLINEARDSIOXML_H
#define LAGRANGIANLINEARDSIOXML_H

#include "LagrangianDSIOXML.h"

const std::string LLDSIO_H = "H";
const std::string LLDSIO_B = "b";

/** \class LagrangianLinearDSIOXML
*   \brief This class manages LagrangianLinear DSInputOutput data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.1.4.
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

  /** \fn LagrangianLinearDSIOXML(xmlNode * dsioNode)
  *   \brief Build a LagrangianLinearDSIOXML object from a DOM tree describing a DSIO with LagrangianLinear type
  *   \param LagrangianLinear : the LagrangianLinearDSIO DOM tree
  *   \exception XMLException : if a property of the LagrangianLinear DSIO lacks in the DOM tree
  */
  LagrangianLinearDSIOXML(xmlNode * dsioNode);

  //////////////////////////////////////////////

  /** \fn SimpleMatrix getH()
  *   \brief Return the H of the LagrangianLinearDSIOXML
  *   \return The H SimpleMatrix of the LagrangianLinearDSIOXML
  */
  inline SimpleMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }


  /** \fn SimpleVector getB()
  *   \brief Return b vector of the LagrangianLinearDSIOXML
  *   \return SimpleVector : b vector of the LagrangianLinearDSIOXML
  */
  inline SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->bNode);
  }

  /** \fn void setH(SiconosMatrix *matrix)
  *   \brief Change the H matrix value (in xml file or external data file switch his origin position)
  *   \param SiconosMatrix matrix : the new value for H matrix
  */
  void setH(SiconosMatrix *matrix);

  /** \fn void setB(SiconosVector *vector)
  *   \brief Change the b vector value (in xml file or external data file switch his origin position)
  *   \param SiconosVector vector : the new value for b vector
  */
  void setB(SiconosVector *vector);


private:
  //Nodes
  xmlNode * HNode;
  xmlNode * bNode;
};

#endif // LAGRANGIANLINEARDSIOXML_H

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

/*! \file LinearBCXML.h

*/
#ifndef __LINEARBCXML__
#define __LINEARBCXML__

#include "BoundaryConditionXML.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"

//Tags
const std::string LINEARBC_OMEGA = "Omega";
const std::string LINEARBC_OMEGA0 = "Omega0";
const std::string LINEARBC_OMEGAT = "OmegaT";


//! XML management for LinearBC
/**  \author SICONOS Development Team - copyright INRIA
*   \version 2.1.1.
*   \date 05/25/2004
*
*
*
* LinearBCXML allows to manage data of a LinearBC DOM tree.
*/
class LinearBCXML : public BoundaryConditionXML
{
public:

  LinearBCXML();

  /** Build a LinearBCXML object from a DOM tree describing a LinearBC
  *   \param xmlNode * LinearBCNode : the LinearBC DOM tree
  */
  LinearBCXML(xmlNode * LinearBCNode);

  ~LinearBCXML();

  /** Return Omega of the LinearBCXML
  *   \return SimpleVector : Omega of LinearBCXML
  */
  inline /*SiconosVector*/SimpleVector getOmega()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->omegaNode);
  }

  /** Return the Omega0 of the LinearBCXML
  *   \return The Omega0 SimpleMatrix of the LinearBCXML
  */
  inline SimpleMatrix getOmega0()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->omega0Node);
  }

  /** Return the OmegaT of the LinearBCXML
  *   \return The OmegaT SimpleMatrix of the LinearBCXML
  */
  inline SimpleMatrix getOmegaT()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->omegaTNode);
  }

  /** allows to save the Omega of the LinearBCXML
  *   \param The Omega SiconosVector to save
  */
  inline void setOmega(SiconosVector *v)
  {
    if (this->omegaNode == NULL)
    {
      this->omegaNode = SiconosDOMTreeTools::createVectorNode(this->rootBCNode, LINEARBC_OMEGA, *v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->omegaNode, *v);
  }

  /** allows to save the Omega0 of the LinearBCXML
  *   \param The Omega0 SiconosMatrix to save
  */
  inline void setOmega0(SiconosMatrix *m)
  {
    if (this->omega0Node == NULL)
    {
      this->omega0Node = SiconosDOMTreeTools::createMatrixNode(this->rootBCNode, LINEARBC_OMEGA0, *m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->omega0Node, *m);
  }

  /** allows to save the OmegaT of the LinearBCXML
  *   \param The OmegaT SiconosMatrix to save
  */
  inline void setOmegaT(SiconosMatrix *m)
  {
    if (this->omegaTNode == NULL)
    {
      this->omegaTNode = SiconosDOMTreeTools::createMatrixNode(this->rootBCNode, LINEARBC_OMEGAT, *m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->omegaTNode, *m);
  }


  /** makes the operations to add a BoundaryCondition to the DynamicalSystemXML
  *  \param xmlNode* : the root node of this BoundaryCondition
  //     *  \param BoundaryCondition* : the BoundaryCondition of the DS
  */
  void updateBoundaryConditionXML(xmlNode* node/*, BoundaryCondition* bc*/);


private:

  //Nodes
  xmlNode * omegaNode;
  xmlNode * omega0Node;
  xmlNode * omegaTNode;

  //Methods
  /** load the different properties of a LinearBC
  *   \exception XMLException : if a property of the LinearBC lacks in the DOM tree
  */
  void loadLinearBCProperties();
};


#endif

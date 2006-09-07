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

/** \class LinearBCXML
*   \brief This class manages Linear BC data part
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.3.0.
*   \date 05/25/2004
*
*
*
* LinearBCXML allows to manage data of a LinearBC DOM tree.
*/

#ifndef __LINEARBCXML__
#define __LINEARBCXML__

#include "BoundaryConditionXML.h"

//Tags
const std::string LINEARBC_OMEGA = "Omega";
const std::string LINEARBC_OMEGA0 = "Omega0";
const std::string LINEARBC_OMEGAT = "OmegaT";


class LinearBCXML : public BoundaryConditionXML
{
public:

  LinearBCXML();

  /** \fn LinearBCXML(xmlNode * LinearBCNode)
  *   \brief Build a LinearBCXML object from a DOM tree describing a LinearBC
  *   \param xmlNode * LinearBCNode : the LinearBC DOM tree
  */
  LinearBCXML(xmlNode * LinearBCNode);

  ~LinearBCXML();

  /** \fn SimpleVector getOmega()
  *   \brief Return Omega of the LinearBCXML
  *   \return SimpleVector : Omega of LinearBCXML
  */
  inline /*SiconosVector*/SimpleVector getOmega()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(this->omegaNode);
  }

  /** \fn SimpleMatrix getOmega0()
  *   \brief Return the Omega0 of the LinearBCXML
  *   \return The Omega0 SimpleMatrix of the LinearBCXML
  */
  inline SimpleMatrix getOmega0()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->omega0Node);
  }

  /** \fn SimpleMatrix getOmegaT()
  *   \brief Return the OmegaT of the LinearBCXML
  *   \return The OmegaT SimpleMatrix of the LinearBCXML
  */
  inline SimpleMatrix getOmegaT()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->omegaTNode);
  }

  /** \fn void setOmega(SiconosVector *v)
  *   \brief allows to save the Omega of the LinearBCXML
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

  /** \fn void setOmega0(SiconosMatrix *m)
  *   \brief allows to save the Omega0 of the LinearBCXML
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

  /** \fn void setOmegaT(SiconosMatrix *m)
  *   \brief allows to save the OmegaT of the LinearBCXML
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


  /** \fn void updateBoundaryConditionXML( xmlNode* node)//, BoundaryCondition* bc )
   *  \brief makes the operations to add a BoundaryCondition to the DynamicalSystemXML
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
  /** \fn loadLinearBCProperties(xmlNode * LinearBCnode)
  *   \brief load the different properties of a LinearBC
  *   \exception XMLException : if a property of the LinearBC lacks in the DOM tree
  */
  void loadLinearBCProperties();
};


#endif

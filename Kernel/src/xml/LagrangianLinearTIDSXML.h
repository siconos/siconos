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
/** \class LagrangianLinearTIDSXML
 *   \brief This class manages Lagrangian TIDS data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.1.4.
 *   \date 05/11/2004
 *
 *
 *
 * LagrangianLinearTIDSXML allows to manage data of a LagrangianLinearTIDS DOM tree.
 */


#ifndef __LAGRANGIANTIDSXML__
#define __LAGRANGIANTIDSXML__

#include "LagrangianDSXML.h"

const std::string LTIDS_K = "K";
const std::string LTIDS_C = "C";


class LagrangianLinearTIDSXML : public LagrangianDSXML
{
private:

  xmlNode * KNode;
  xmlNode * CNode;

public:

  /** \fn LagrangianLinearTIDSXML()
   *   \brief default constructor
   */
  LagrangianLinearTIDSXML();

  /** \fn LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, int number)
   *   \brief Build a LagrangianLinearTIDSXML object from a DOM tree describing a LagrangianLinearTIDS
   *   \param LagrangianLinearTIDSNode : the LagrangianLinearTIDS DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP LagrangianLinearTIDS have boundary condition
   */
  LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, bool isBVP);

  /** \fn const SimpleMatrix getK() const
   *   \brief Return the K of the LagrangianLinearTIDSXML
   *   \return The K SimpleMatrix of the LagrangianLinearTIDSXML
   */
  inline const SimpleMatrix getK() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(KNode);
  }

  /** \fn void setK(SiconosMatrix *m)
   *   \brief allows to save the K of the LagrangianLinearTIDSXML
   *   \return The K SiconosMatrix to save
   */
  inline void setK(const SiconosMatrix& m)
  {
    if (KNode == NULL)
    {
      KNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LTIDS_K, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(KNode, m);
  }

  /** \fn const SimpleMatrix getC() const
   *   \brief Return the C of the LagrangianLinearTIDSXML
   *   \return The C SimpleMatrix of the LagrangianLinearTIDSXML
   */
  inline const SimpleMatrix getC() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(CNode);
  }

  /** \fn void setC(SiconosMatrix *m)
   *   \brief allows to save the C of the LagrangianLinearTIDSXML
   *   \return The C SiconosMatrix to save
   */
  inline void setC(const SiconosMatrix& m)
  {
    if (CNode == NULL)
    {
      CNode = SiconosDOMTreeTools::createMatrixNode(rootDynamicalSystemXMLNode, LTIDS_C, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(CNode, m);
  }


  /** \fn void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition*)
   *   \brief makes the operations to add a DynamicalSystem to the NonSmoothDynamicalSystemXML
   *   \param xmlNode* : the root node of this DynamicalSystem
   *   \param DynamicalSystem* : the DynamicalSystem of this DynamicalSystemXML
   *   \param BoundaryCondition* : the BoundaryCondition of the DS if the NonSmoothDynamicalSystem is BVP (optional)
   */
  void updateDynamicalSystemXML(xmlNode*, DynamicalSystem*, BoundaryCondition* bc = NULL);

  /** \fn bool hasK()
   *  \brief determines if K is defined in the DOM tree
   *  \return bool : true if K is defined, false otherwise
   */
  inline const bool hasK() const
  {
    return (KNode != NULL);
  }

  /** \fn bool hasC()
   *  \brief determines if C is defined in the DOM tree
   *  \return bool : true if C is defined, false otherwise
   */
  inline const bool hasC() const
  {
    return (CNode != NULL);
  }
};

#endif

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

/*! \file LagrangianLinearTIDSXML.h

*/

#ifndef __LAGRANGIANTIDSXML__
#define __LAGRANGIANTIDSXML__

#include "LagrangianDSXML.h"

const std::string LTIDS_K = "K";
const std::string LTIDS_C = "C";

/** XML management for LagrangianLinearTIDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.1.
 *   \date 05/11/2004
 *
 *
 *
 * LagrangianLinearTIDSXML allows to manage data of a LagrangianLinearTIDS DOM tree.
 */
class LagrangianLinearTIDSXML : public LagrangianDSXML
{
private:

  xmlNodePtr KNode;
  xmlNodePtr CNode;

  /** default constructor
   */
  LagrangianLinearTIDSXML();

public:

  /** Build a LagrangianLinearTIDSXML object from a DOM tree describing a LagrangianLinearTIDS
   *   \param LagrangianLinearTIDSNode : the LagrangianLinearTIDS DOM tree
   *   \param bool isBVP : if NonSmoothDynamicalSystem is BVP LagrangianLinearTIDS have boundary condition
   */
  LagrangianLinearTIDSXML(xmlNode * LagrangianLinearTIDSNode, bool isBVP);

  /** Destructor */
  ~LagrangianLinearTIDSXML();

  /** Return the K of the LagrangianLinearTIDSXML
  *   \return The K SimpleMatrix of the LagrangianLinearTIDSXML
  */
  inline const SimpleMatrix getK() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(KNode);
  }

  /** allows to save the K of the LagrangianLinearTIDSXML
  *   \return The K SiconosMatrix to save
  */
  inline void setK(const SiconosMatrix& m)
  {
    if (KNode == NULL)
    {
      KNode = SiconosDOMTreeTools::createMatrixNode(rootNode, LTIDS_K, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(KNode, m);
  }

  /** Return the C of the LagrangianLinearTIDSXML
  *   \return The C SimpleMatrix of the LagrangianLinearTIDSXML
  */
  inline const SimpleMatrix getC() const
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(CNode);
  }

  /** allows to save the C of the LagrangianLinearTIDSXML
  *   \return The C SiconosMatrix to save
  */
  inline void setC(const SiconosMatrix& m)
  {
    if (CNode == NULL)
    {
      CNode = SiconosDOMTreeTools::createMatrixNode(rootNode, LTIDS_C, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(CNode, m);
  }

  /** determines if K is defined in the DOM tree
  *  \return bool : true if K is defined, false otherwise
  */
  inline const bool hasK() const
  {
    return (KNode != NULL);
  }

  /** determines if C is defined in the DOM tree
  *  \return bool : true if C is defined, false otherwise
  */
  inline const bool hasC() const
  {
    return (CNode != NULL);
  }
};

#endif

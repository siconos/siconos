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

/*! \file
*/

#ifndef _LagrangianLinearECXML_
#define _LagrangianLinearECXML_

#include "LagrangianECXML.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"

const std::string LLEC_H = "H";
const std::string LLEC_B = "b";

//! XML management for LagrangianLinear Relation
/**  \author SICONOS Development Team - copyright INRIA
*   \version 2.1.1.
*   \date 05/25/2004
*
*
*
* LagrangianLinearECXML allows to manage data of a LagrangianLinearEC DOM tree.
*/
class LagrangianLinearECXML : public LagrangianECXML
{
public:

  LagrangianLinearECXML();

  /** Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LagrangianLinearECXML(xmlNode*, std::vector<int>);

  ~LagrangianLinearECXML();

  /** Return the H of the LagrangianLinearECXML
  *   \return The H SimpleMatrix of the LagrangianLinearECXML
  */
  inline SimpleMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(HNode);
  }


  /** Return b vector of the LLRelationXML
  *   \return SimpleVector : b vector of the LLRelationXML
  */
  inline SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
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

#endif

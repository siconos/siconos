/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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

/** \class LagrangianLinearECXML
*   \brief This class manages LagrangianLinear Relation data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.1.3.
*   \date 05/25/2004
*
*
*
* LagrangianLinearECXML allows to manage data of a LagrangianLinearEC DOM tree.
*/


#ifndef _LagrangianLinearECXML_
#define _LagrangianLinearECXML_

#include "LagrangianECXML.h"

const std::string LLEC_H = "H";
const std::string LLEC_B = "b";

class LagrangianLinearECXML : public LagrangianECXML
{
public:

  LagrangianLinearECXML();

  /** \fn LagrangianLinearECXML(xmlNode * , vector<int> )
  *   \brief Build a EqualityConstraintXML object from a DOM tree describing a EqualityConstraint
  *   \param xmlNode* : the EqualityConstraint DOM tree
  *   \param vector<int>  : vector of DSXML numbers to verify DS concerned by the EqualityConstraint (identified by number) exists
  */
  LagrangianLinearECXML(xmlNode*, std::vector<int>);

  ~LagrangianLinearECXML();

  /** \fn SimpleMatrix getH()
  *   \brief Return the H of the LagrangianLinearECXML
  *   \return The H SimpleMatrix of the LagrangianLinearECXML
  */
  inline SimpleMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->HNode);
  }


  /** \fn SimpleVector getB()
  *   \brief Return b vector of the LLRelationXML
  *   \return SimpleVector : b vector of the LLRelationXML
  */
  inline /*SiconosVector*/SimpleVector getB()
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

#endif

/* Siconos version 1.0, Copyright INRIA 2005.
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

/** \class LagrangianLinearRXML
 *   \brief This class manages LagrangianLinear Relation data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/25/2004
 *
 *
 *
 * LagrangianLinearRXML allows to manage data of a LLRelation DOM tree.
 */

#ifndef __LLRelationXML__
#define __LLRelationXML__

#include "LagrangianRXML.h"

const std::string  LLR_H = "H";
const std::string  LLR_B = "b";

class LagrangianLinearRXML : public LagrangianRXML
{

private:

  //Nodes
  xmlNode * HNode;
  xmlNode * bNode;
  xmlNode * DNode;

public:

  LagrangianLinearRXML();

  /** \fn LagrangianLinearRXML(xmlNode * LLRelationNode)
   *   \brief Build a LagrangianLinearRXML object from a DOM tree describing a Relation with LL type
   *   \param LagrangianLinearRXML : the LagrangianLinearR DOM tree
   *   \exception XMLException : if a property of the LagrangianLinear Relation lacks in the DOM tree
   */
  LagrangianLinearRXML(xmlNode * LLRelationNode);

  ~LagrangianLinearRXML();

  /** \fn SiconosMatrix getH()
   *   \brief Return the H of the LLRelationXML
   *   \return The H SiconosMatrix of the LLRelationXML
   */
  inline SiconosMatrix getH()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(HNode);
  }


  /** \fn SimpleVector getB()
   *   \brief Return b vector of the LLRelationXML
   *   \return SimpleVector : b vector of the LLRelationXML
   */
  inline /*SiconosVector*/SimpleVector getB()
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(bNode);
  }

  /** \fn SiconosMatrix getD()
   *   \brief Return the D of the LLRelationXML
   *   \return The D SiconosMatrix of the LLRelationXML
   */
  inline SiconosMatrix getD()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(DNode);
  }

  /** \fn void setH(SiconosMatrix *matrix)
   *   \brief Change the H matrix value (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for H matrix
   */
  void setH(const SiconosMatrix&);

  /** \fn void setB(SiconosVector *vector)
   *   \brief Change the b vector value (in xml file or external data file switch his origin position)
   *   \param SiconosVector vector : the new value for b vector
   */
  void setB(const SiconosVector&);

  /** \fn void setD(SiconosMatrix *matrix)
   *   \brief Change the D matrix value (in xml file or external data file switch his origin position)
   *   \param SiconosMatrix matrix : the new value for D matrix
   */
  void setD(const SiconosMatrix&);


  /** \fn bool hasB() const
   *   \brief return true if b is given in xmlfile
   */
  inline bool hasB() const
  {
    return (!(bNode == NULL));
  }

  /** \fn bool hasD() const
   *   \brief return true if D is given in xmlfile
   */
  inline bool hasD() const
  {
    return (!(DNode == NULL));
  }
};

#endif

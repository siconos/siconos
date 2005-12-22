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
/** \class FrictionContact2DXML
 *   \brief This class manages Lagrangian FrictionContact2D data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/18/2004
 *
 *
 *
 * FrictionContact2DXML allows to manage data of a FrictionContact2D DOM tree.
 */

#ifndef __FrictionContact2DXML__
#define __FrictionContact2DXML__

#include "OneStepNSProblemXML.h"

class OneStepNSProblem;

const std::string  FrictionContact2D_M = "M";
const std::string  FrictionContact2D_Q = "q";


class FrictionContact2DXML : public OneStepNSProblemXML
{
public:
  FrictionContact2DXML();

  /** \fn FrictionContact2DXML(xmlNode * FrictionContact2DNode)
   *   \brief Build a FrictionContact2DXML object from a DOM tree describing a FrictionContact2D
   *   \param FrictionContact2DNode : the FrictionContact2D DOM tree
   *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
   *   \exception XMLException : if a property of the FrictionContact2D lacks in the DOM tree
   */
  FrictionContact2DXML(xmlNode * FrictionContact2DNode, std::vector<int> definedInteractionNumbers);

  ~FrictionContact2DXML();

  /** \fn SiconosMatrix getM()
   *   \brief Return M
   *   \return The M SiconosMatrix of the FrictionContact2D
   */
  inline SiconosMatrix getM()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** \fn SimpleVector getQ()
   *   \brief Return vector q
   *   \return SimpleVector : q vector of the FrictionContact2D
   */
  inline SimpleVector getQ()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(qNode);
  }

  /** \fn void setM(const SiconosMatrix &m)
   *   \brief save M
   *   \param The M SiconosMatrix to save
   */
  inline void setM(const SiconosMatrix &m)
  {
    if (hasM() == false)
    {
      MNode = SiconosDOMTreeTools::createMatrixNode(problemTypeNode, FrictionContact2D_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, m);
  }

  /** \fn void setQ(const SiconosVector &q)
   *   \brief save q
   *   \param The q SiconosVector to save
   */
  inline void setQ(const SiconosVector& q)
  {
    if (hasQ() == false)
    {
      qNode = SiconosDOMTreeTools::createVectorNode(problemTypeNode, FrictionContact2D_Q, q);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, q);
  }

  /** \fn bool hasM()
   *  \brief returns true if MNode is defined
   *  \return true if MNode is defined
   */
  inline bool hasM()
  {
    return (MNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if qNode is defined
   *  \return true if qNode is defined
   */
  inline bool hasQ()
  {
    return (qNode != NULL);
  }

  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
   *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
   *   \param xmlNode* : the root node of the OneStepNSProblemXML
   *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
   */
  void updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb);


private:

  //Nodes
  xmlNode * MNode;
  xmlNode * qNode;
};


#endif

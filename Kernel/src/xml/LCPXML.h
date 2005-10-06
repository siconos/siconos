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
/** \class LCPXML
 *   \brief This class manages Lagrangian LCP data
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/18/2004
 *
 *
 *
 * LCPXML allows to manage data of a LCP DOM tree.
 */

#ifndef __LCPXML__
#define __LCPXML__

#include "OneStepNSProblemXML.h"

class OneStepNSProblem;

const std::string  LCP_M = "M";
const std::string  LCP_Q = "q";


class LCPXML : public OneStepNSProblemXML
{
public:
  LCPXML();

  /** \fn LCPXML(xmlNode * LCPNode)
   *   \brief Build a LCPXML object from a DOM tree describing a LCP
   *   \param LCPNode : the LCP DOM tree
   *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
   *   \exception XMLException : if a property of the LCP lacks in the DOM tree
   */
  LCPXML(xmlNode * LCPNode, std::vector<int> definedInteractionNumbers);

  ~LCPXML();

  /** \fn SiconosMatrix getM()
   *   \brief Return M
   *   \return The M SiconosMatrix of the LCP
   */
  inline SiconosMatrix getM()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** \fn SimpleVector getQ()
   *   \brief Return vector q
   *   \return SimpleVector : q vector of the LCP
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
      MNode = SiconosDOMTreeTools::createMatrixNode(rootNSProblemXMLNode, LCP_M, m);
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
      qNode = SiconosDOMTreeTools::createVectorNode(rootNSProblemXMLNode, LCP_Q, q);
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

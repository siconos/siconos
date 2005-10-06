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

/** \class QPXML
*   \brief This class manages Lagrangian QP data
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date 05/18/2004
*
*
*
* QPXML allows to manage data of a QP DOM tree.
*/


#ifndef __QPXMLDEF__
#define __QPXMLDEF__

#include "OneStepNSProblemXML.h"

class OneStepNSProblem;

const std::string  QP_Q = "Q";
const std::string  QP_P = "p";


class QPXML : public OneStepNSProblemXML
{
public:
  QPXML();

  /** \fn QPXML(xmlNode * QPNode)
  *   \brief Build a QPXML object from a DOM tree describing a QP
  *   \param QPNode : the QP DOM tree
  *   \param vector<int> definedInteractionNumbers : the Interaction numbers effectivly defined in the model
  *   \exception XMLException : if a property of the QP lacks in the DOM tree
  */
  QPXML(xmlNode * QPNode, std::vector<int> definedInteractionNumbers);

  ~QPXML();

  /** \fn SiconosMatrix getQ()
  *   \brief Return Q
  *   \return The Q SiconosMatrix of the QP
  */
  inline SiconosMatrix getQ()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->QNode);
  }

  /** \fn SimpleVector getP()
  *   \brief Return p
  *   \return SimpleVector :  vector p of the QP
  */
  inline /*SiconosVector*/SimpleVector getP()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->pNode);
  }

  /** \fn void setQ(const SiconosMatrix& m)
  *   \brief allows  to save Q
  *   \param The Q SiconosMatrix to save
  */
  inline void setQ(const SiconosMatrix& m)
  {
    if (this->hasQ() == false)
    {
      this->QNode = SiconosDOMTreeTools::createMatrixNode(this->rootNSProblemXMLNode, QP_Q, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->QNode, m);
  }

  /** \fn void setP(const SiconosVector&v)
  *   \brief allows to save p
  *   \param SimpleVector* : vector p to save
  */
  inline void setP(const SiconosVector&v)
  {
    if (this->hasP() == false)
    {
      this->pNode = SiconosDOMTreeTools::createVectorNode(this->rootNSProblemXMLNode, QP_P, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->pNode, v);
  }

  /** \fn bool hasP()
   *  \brief returns true if pNode is defined
   *  \return true if pNode is defined
   */
  inline bool hasP()
  {
    return (this->pNode != NULL);
  }

  /** \fn bool hasQ()
   *  \brief returns true if QNode is defined
   *  \return true if QNode is defined
   */
  inline bool hasQ()
  {
    return (this->QNode != NULL);
  }

  /** \fn void updateOneStepNSProblemXML( xmlNode* node, OneStepNSProblemXML* str )
  *   \brief makes the operations to create a OneStepNSProblemXML to the StrategyXML
  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  */
  void updateOneStepNSProblemXML(xmlNode* node, OneStepNSProblem* osnspb);


private:
  //Nodes
  xmlNode * QNode;
  xmlNode * pNode;

};

#endif

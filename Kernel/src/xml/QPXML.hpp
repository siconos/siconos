/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
/*! \file
 */



#ifndef __QPXMLDEF__
#define __QPXMLDEF__

#include "OneStepNSProblemXML.hpp"
#include "SimpleMatrix.hpp"
#include "SimpleVector.hpp"

class OneStepNSProblem;

const std::string  QP_Q = "Q";
const std::string  QP_P = "p";

/** XML management for QP
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/18/2004
 */
class QPXML : public OneStepNSProblemXML
{
public:
  QPXML();

  /** Build a QPXML object from a DOM tree describing a QP
   *   \param QPNode : the QP DOM tree
   *   \exception XMLException : if a property of the QP lacks in the DOM tree
   */
  QPXML(xmlNode * QPNode);

  ~QPXML();

  /** Return Q
   *   \return The Q SimpleMatrix of the QP
   */
  inline SimpleMatrix getQ()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(this->QNode);
  }

  /** Return p
   *   \return SimpleVector :  vector p of the QP
   */
  inline /*SiconosVector*/SimpleVector getP()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(this->pNode);
  }

  /** allows  to save Q
   *   \param The Q SiconosMatrix to save
   */
  inline void setQ(const SiconosMatrix& m)
  {
    if (this->hasQ() == false)
    {
      this->QNode = SiconosDOMTreeTools::createMatrixNode(this->rootNode, QP_Q, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(this->QNode, m);
  }

  /** allows to save p
   *   \param SimpleVector* : vector p to save
   */
  inline void setP(const SiconosVector&v)
  {
    if (this->hasP() == false)
    {
      this->pNode = SiconosDOMTreeTools::createVectorNode(this->rootNode, QP_P, v);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(this->pNode, v);
  }

  /** returns true if pNode is defined
   *  \return true if pNode is defined
   */
  inline bool hasP()
  {
    return (this->pNode);
  }

  /** returns true if QNode is defined
   *  \return true if QNode is defined
   */
  inline bool hasQ()
  {
    return (this->QNode);
  }

  /** makes the operations to create a OneStepNSProblemXML to the SimulationXML
   *   \param xmlNode* : the root node of the OneStepNSProblemXML
   *   \param SP::OneStepNSProblem : the OneStepNSProblem of this OneStepNSProblemXML
   */
  void updateOneStepNSProblemXML(xmlNode* node, SP::OneStepNSProblem osnspb);


private:
  //Nodes
  xmlNode * QNode;
  xmlNode * pNode;

};
DEFINE_SPTR(QPXML);
#endif

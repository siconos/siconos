/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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

#ifndef __FrictionContactXML__
#define __FrictionContactXML__

#include "OneStepNSProblemXML.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"

class OneStepNSProblem;

const std::string  FrictionContact_M = "M";
const std::string  FrictionContact_Q = "q";


/** XML management for FrictionContact
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.0.
 *   \date 05/18/2004
 *
 */
class FrictionContactXML : public OneStepNSProblemXML
{
public:
  FrictionContactXML();

  /** Build a FrictionContactXML object from a DOM tree describing a FrictionContact
  *   \param FrictionContactNode : the FrictionContact DOM tree
  *   \exception XMLException : if a property of the FrictionContact lacks in the DOM tree
  */
  FrictionContactXML(xmlNode * FrictionContactNode);

  ~FrictionContactXML();

  /** Return M
  *   \return The M SimpleMatrix of the FrictionContact
  */
  inline SimpleMatrix getM()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(MNode);
  }

  /** Return vector q
  *   \return SimpleVector : q vector of the FrictionContact
  */
  inline SimpleVector getQ()
  {
    return SiconosDOMTreeTools::getSiconosVectorValue(qNode);
  }

  /** save M
  *   \param The M SiconosMatrix to save
  */
  inline void setM(const SiconosMatrix &m)
  {
    if (hasM() == false)
    {
      MNode = SiconosDOMTreeTools::createMatrixNode(problemTypeNode, FrictionContact_M, m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(MNode, m);
  }

  /** save q
  *   \param The q SiconosVector to save
  */
  inline void setQ(const SiconosVector& q)
  {
    if (hasQ() == false)
    {
      qNode = SiconosDOMTreeTools::createVectorNode(problemTypeNode, FrictionContact_Q, q);
    }
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(qNode, q);
  }

  /** returns true if MNode is defined
  *  \return true if MNode is defined
  */
  inline bool hasM()
  {
    return (MNode != NULL);
  }

  /** returns true if qNode is defined
  *  \return true if qNode is defined
  */
  inline bool hasQ()
  {
    return (qNode != NULL);
  }

  /** makes the operations to create a OneStepNSProblemXML to the SimulationXML
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

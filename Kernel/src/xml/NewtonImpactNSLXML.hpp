/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#ifndef __NewtonImpactNSLXML__
#define __NewtonImpactNSLXML__

#include "NonSmoothLawXML.hpp"

const std::string NEWTON_E = "e";

/** XML management for NewtonImpactNSL
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) June 29, 2004
 *
 * NewtonImpactNSLXML manage input/output xml data for Newton Impact non-smooth law.
 *
 */

class NewtonImpactNSLXML : public NonSmoothLawXML
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonImpactNSLXML);

  /** node that handle restitution coefficient for the law */
  xmlNode * ENode;

  /** default constructor (private)
  */
  NewtonImpactNSLXML();

public:
  /** Build a NewtonImpactNSLXML object from a DOM tree describing a Law with Relay type
  *   \param NewtonImpactNSLLawNode : theNewtonImpactNSLLaw DOM tree
  *   \exception XMLException : if a property of the NewtonImpactNSL  lacks in the DOM tree
  */
  NewtonImpactNSLXML(xmlNode *);

  /** Return the e of the NSLaw
  *   \return a double
  */
  inline double getE() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(ENode);
  }

  /** set e value
  *   \param a double
  */
  inline void setE(const double& e)
  {
    if (!hasE())
      ENode = SiconosDOMTreeTools::createDoubleNode(rootNode, NEWTON_E, e);
    else SiconosDOMTreeTools::setDoubleContentValue(ENode, e);
  }

  /** returns true if ENode is defined
  *  \return true if ENode is defined
  */
  inline bool hasE() const
  {
    return (ENode);
  }

};

#endif

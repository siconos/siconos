/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
/** \class NewtonImpactLawNSLXML
 *  \brief  This class manages NewtonImpactLawNSL data part
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Creation) June 29, 2004
 *
 *
 *
 * NewtonImpactLawNSLXML allows to manage data of a NewtonImpactLawNSL DOM tree.
 * \bug
 *  \warning
 */

#ifndef __NewtonImpactLawNSLXML__
#define __NewtonImpactLawNSLXML__

#include "NonSmoothLawXML.h"

const std::string NEWTON_E = "e";


class NewtonImpactLawNSLXML : public NonSmoothLawXML
{
public:
  NewtonImpactLawNSLXML();

  /** \fn NewtonImpactLawNSLXML(xmlNode *NewtonImpactLawNSLLawNode)
  *   \brief Build a NewtonImpactLawNSLXML object from a DOM tree describing a Law with Relay type
  *   \param NewtonImpactLawNSLLawNode : theNewtonImpactLawNSLLaw DOM tree
  *   \exception XMLException : if a property of the NewtonImpactLawNSL  lacks in the DOM tree
  */
  NewtonImpactLawNSLXML(xmlNode * NewtonImpactLawNSLNode);

  /** \fn double getE()
  *   \brief Return the E of the NSLaw
  *   \return The E double value of the coefficient of restitution
  */
  inline double getE()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->ENode);
  }

  /** \fn void setE(double e)
  *   \brief Return the E of NSLaw
  *   \return The E double value of the coefficient of restitution
  */
  inline void setE(double e)
  {
    if (this->hasE() == false)
    {
      this->ENode = SiconosDOMTreeTools::createDoubleNode(this->rootNSLawXMLNode, NEWTON_E, e);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->ENode, e);
  }

  /** \fn bool hasE()
   *  \brief returns true if ENode is defined
   *  \return true if ENode is defined
   */
  inline bool hasE()
  {
    return (this->ENode != NULL);
  }

private:
  xmlNode * ENode;

};

#endif

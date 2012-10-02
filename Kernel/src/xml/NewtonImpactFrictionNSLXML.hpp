/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

#ifndef __NewtonImpactFrictionNSLXML__
#define __NewtonImpactFrictionNSLXML__

#include "SiconosPointers.hpp"
#include "NonSmoothLawXML.hpp"

/** XML management for NewtonImpactFrictionNSL
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) March 22, 2005
 *
 *
 *
 * NewtonImpactFrictionNSLXMLdata xml management for NewtonImpactFrictionNSL
 */
class NewtonImpactFrictionNSLXML : public NonSmoothLawXML
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonImpactFrictionNSLXML);

  xmlNode * enNode;
  xmlNode * etNode;
  xmlNode * muNode;

public:
  /** default constructor
  */
  NewtonImpactFrictionNSLXML();

  /** Build a NewtonImpactFrictionNSLXML object using DOM tree data loading
  *   \param xmlNode* : the NewtonImpactFrictionNSL node in the DOM tree
  */
  NewtonImpactFrictionNSLXML(xmlNodePtr);

  /** return the En of the NSLaw
  *   \return a double
  */
  inline double getEn() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(enNode);
  }

  /** set en value
  *   \param a double
  */
  inline void setEn(const double& en)
  {
    if (!hasEn())
      enNode = SiconosDOMTreeTools::createDoubleNode(rootNode, "en", en);
    else SiconosDOMTreeTools::setDoubleContentValue(enNode, en);
  }

  /** returns true if enNode is defined
  *  \return a bool
  */
  inline bool hasEn() const
  {
    return (enNode);
  }

  /** return the Et of the NSLaw
  *   \return a double
  */
  inline double getEt() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(etNode);
  }

  /** set et value
  *   \param a double
  */
  inline void setEt(const double& et)
  {
    if (!hasEt())
      etNode = SiconosDOMTreeTools::createDoubleNode(rootNode, "et", et);
    else SiconosDOMTreeTools::setDoubleContentValue(etNode, et);
  }

  /** returns true if etNode is defined
  *  \return a bool
  */
  inline bool hasEt() const
  {
    return (etNode);
  }

  /** return mu value
  *   \return a double
  */
  inline double getMu() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(muNode);
  }

  /** set mu value
  *   \param a double
  */
  inline void setMu(const double& mu)
  {
    if (!hasMu())
      muNode = SiconosDOMTreeTools::createDoubleNode(rootNode, "mu", mu);
    else SiconosDOMTreeTools::setDoubleContentValue(muNode, mu);
  }

  /** return true if muNode is defined
  *  \return a bool
  */
  inline bool hasMu() const
  {
    return (muNode);
  }

};

TYPEDEF_SPTR(NewtonImpactFrictionNSLXML)

#endif

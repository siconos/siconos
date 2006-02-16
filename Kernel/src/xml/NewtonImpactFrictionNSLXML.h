/* Siconos-Kernel version 1.1.1, Copyright INRIA 2005-2006.
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
/** \class NewtonImpactFrictionNSLXML
 *  \brief  This class manages NewtonImpactFrictionNSL data part
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.1.
 *  \date (Creation) March 22, 2005
 *
 *
 *
 * NewtonImpactFrictionNSLXMLdata xml management for NewtonImpactFrictionNSL
 * \bug
 */

#ifndef __NewtonImpactFrictionNSLXML__
#define __NewtonImpactFrictionNSLXML__

#include "NonSmoothLawXML.h"

class NewtonImpactFrictionNSLXML : public NonSmoothLawXML
{
private:
  xmlNode * enNode;
  xmlNode * etNode;
  xmlNode * muNode;

public:
  /** \fn NewtonImpactFrictionNSLXML()
   *   \brief default constructor
   */
  NewtonImpactFrictionNSLXML();

  /** \fn NewtonImpactFrictionNSLXML(xmlNodePtr)
   *   \brief Build a NewtonImpactFrictionNSLXML object using DOM tree data loading
   *   \param xmlNode* : the NewtonImpactFrictionNSL node in the DOM tree
   */
  NewtonImpactFrictionNSLXML(xmlNodePtr);

  /** \fn const double getEn() const
   *   \brief return the En of the NSLaw
   *   \return a double
   */
  inline const double getEn() const
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(enNode);
  }

  /** \fn void setEn(const double& en)
   *   \brief set en value
   *   \param a double
   */
  inline void setEn(const double& en)
  {
    if (!hasEn())
      enNode = SiconosDOMTreeTools::createDoubleNode(rootNSLawXMLNode, "en", en);
    else SiconosDOMTreeTools::setDoubleContentValue(enNode, en);
  }

  /** \fn bool hasEn() const
   *  \brief returns true if enNode is defined
   *  \return a bool
   */
  inline bool hasEn() const
  {
    return (enNode != NULL);
  }

  /** \fn const double getEt() const
   *   \brief return the Et of the NSLaw
   *   \return a double
   */
  inline const double getEt() const
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(etNode);
  }

  /** \fn void setEt(const double& et)
   *   \brief set et value
   *   \param a double
   */
  inline void setEt(const double& et)
  {
    if (!hasEt())
      etNode = SiconosDOMTreeTools::createDoubleNode(rootNSLawXMLNode, "et", et);
    else SiconosDOMTreeTools::setDoubleContentValue(etNode, et);
  }

  /** \fn bool hasEt() const
   *  \brief returns true if etNode is defined
   *  \return a bool
   */
  inline bool hasEt() const
  {
    return (etNode != NULL);
  }

  /** \fn const double getMu() const
   *   \brief return mu value
   *   \return a double
   */
  inline const double getMu() const
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(muNode);
  }

  /** \fn void setMu(const double& mu)
   *   \brief set mu value
   *   \param a double
   */
  inline void setMu(const double& mu)
  {
    if (!hasMu())
      muNode = SiconosDOMTreeTools::createDoubleNode(rootNSLawXMLNode, "mu", mu);
    else SiconosDOMTreeTools::setDoubleContentValue(muNode, mu);
  }

  /** \fn bool hasMu() const
   *  \brief return true if muNode is defined
   *  \return a bool
   */
  inline bool hasMu() const
  {
    return (muNode != NULL);
  }

};

#endif

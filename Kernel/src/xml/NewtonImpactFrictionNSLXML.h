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
/** \class NewtonImpactFrictionNSLXML
 *  \brief  This class manages NewtonImpactFrictionNSL data part
*  \author SICONOS Development Team - copyright INRIA
 *  \version 0.1
 *  \date (Creation) March 22, 2005
 *
 *
 *
 * NewtonImpactFrictionNSLXML allows to manage data of a NewtonImpactFrictionNSL DOM tree.
 * \bug
 */

#ifndef __NewtonImpactFrictionNSLXML__
#define __NewtonImpactFrictionNSLXML__

#include "NonSmoothLawXML.h"

const std::string NEWTON_EN = "en";
const std::string NEWTON_ET = "et";
const std::string NEWTON_MU = "mu";

class NewtonImpactFrictionNSLXML : public NonSmoothLawXML
{
public:
  NewtonImpactFrictionNSLXML();

  /** \fn NewtonImpactFrictionNSLXML(xmlNode *)
  *   \brief Build a NewtonImpactFrictionNSLXML object from a DOM tree describing a Law with Relay type
  *   \param xmlNode* : the NewtonImpactFrictionNSLXML node in the DOM tree
  *   \exception XMLException : if a property of the NewtonImpactFrictionNSL  lacks in the DOM tree
  */
  NewtonImpactFrictionNSLXML(xmlNode *);

  /** \fn double getEn()
  *   \brief Return the En of the NSLaw
  *   \return The En double value of the normal coefficient of restitution
  */
  inline double getEn()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->enNode);
  }

  /** \fn void setEn(double en)
  *   \brief Return the En of NSLaw
  *   \return The En double value of the normal coefficient of restitution
  */
  inline void setEn(double en)
  {
    if (this->hasEn() == false)
    {
      this->enNode = SiconosDOMTreeTools::createDoubleNode(this->rootNSLawXMLNode, NEWTON_EN, en);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->enNode, en);
  }

  /** \fn bool hasEn()
   *  \brief returns true if enNode is defined
   *  \return true if enNode is defined
   */
  inline bool hasEn()
  {
    return (this->enNode != NULL);
  }

  /** \fn double getEt()
  *   \brief Return the Et of the NSLaw
  *   \return The Et double value of the tangential coefficient of restitution
  */
  inline double getEt()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->etNode);
  }

  /** \fn void setEt(double et)
  *   \brief Return the Et of NSLaw
  *   \return The Et double value of the tangential coefficient of restitution
  */
  inline void setEt(double et)
  {
    if (this->hasEt() == false)
    {
      this->etNode = SiconosDOMTreeTools::createDoubleNode(this->rootNSLawXMLNode, NEWTON_ET, et);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->etNode, et);
  }

  /** \fn bool hasEt()
   *  \brief returns true if etNode is defined
   *  \return true if etNode is defined
   */
  inline bool hasEt()
  {
    return (this->etNode != NULL);
  }

  /** \fn double getMu()
  *   \brief Return the Mu of the NSLaw
  *   \return The Mu double value of the friction coefficient
  */
  inline double getMu()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->muNode);
  }

  /** \fn void setMu(double mu)
  *   \brief Return the Mu of NSLaw
  *   \return The Mu double value of the friction coefficient
  */
  inline void setMu(double mu)
  {
    if (this->hasMu() == false)
    {
      this->muNode = SiconosDOMTreeTools::createDoubleNode(this->rootNSLawXMLNode, NEWTON_MU, mu);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->muNode, mu);
  }

  /** \fn bool hasMu()
   *  \brief returns true if muNode is defined
   *  \return true if muNode is defined
   */
  inline bool hasMu()
  {
    return (this->muNode != NULL);
  }

private:
  xmlNode * enNode;
  xmlNode * etNode;
  xmlNode * muNode;

};

#endif

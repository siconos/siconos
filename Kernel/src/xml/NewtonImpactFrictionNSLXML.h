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


#include <libxml/tree.h>
#include <string>

#include "NonSmoothLawXML.h"
#include "SiconosDOMTreeTools.h"

using namespace std;

const string NEWTON_EN = "en";
const string NEWTON_ET = "et";
const string NEWTON_MU = "mu";


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

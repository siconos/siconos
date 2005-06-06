
/** \classRelayNSLawXML
 *   \brief This class manages RelayNSLaw data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/14/2004
 *
 *
 *
 * RelayNSLawXML allows to manage data of a RelayNSLaw DOM tree.
 */


#ifndef __RelayNSLawXML__
#define __RelayNSLawXML__

#include "NonSmoothLawXML.h"

const std::string RNSL_C = "c";
const std::string RNSL_D = "d";

class RelayNSLXML : public NonSmoothLawXML
{
public:
  RelayNSLXML();

  /** \fn RelayNSLXML(xmlNode * relayNSLawNode)
   *   \brief Build a RelayNSLXML object from a DOM tree describing a Law with Relay type
   *   \param relayNSLawNode : the relayNSLaw DOM tree
   *   \exception XMLException : if a property of the Relay NS Law lacks in the DOM tree
   */
  RelayNSLXML(xmlNode * relayNSLawNode);


  /** \fn double getC()
   *   \brief Return the C of a relayNSLaw
   *   \return The C double value of the relayNSLaw
   */
  inline double getC()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->CNode);
  }

  /** \fn double getD()
   *   \brief Return the D of a relayNSLaw
   *   \return The D double value of the relayNSLaw
   */
  inline double getD()
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(this->DNode);
  }

  /** \fn void setC(double c)
   *   \brief set the C of a relayNSLaw
   *   \param double : The C double value of the relayNSLaw
   */
  inline void setC(double c)
  {
    if (this->CNode == NULL)
    {
      this->CNode = SiconosDOMTreeTools::createDoubleNode(this->rootNSLawXMLNode, RNSL_C, c);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->CNode, c);
  }

  /** \fn void setD(double d)
   *   \brief set the D of a relayNSLaw
   *   \param double : The D double value of the relayNSLaw
   */
  inline void setD(double d)
  {
    if (this->DNode == NULL)
    {
      this->DNode = SiconosDOMTreeTools::createDoubleNode(this->rootNSLawXMLNode, RNSL_D, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(this->DNode, d);
  }


private:
  xmlNode * CNode;
  xmlNode * DNode;

};

#endif

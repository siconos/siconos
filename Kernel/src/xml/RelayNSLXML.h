
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


#include <libxml/tree.h>
#include <string>

#include "NonSmoothLawXML.h"
#include "SiconosDOMTreeTools.h"

using namespace std;

const string RNSL_C = "c";
const string RNSL_D = "d";


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
//$Log: RelayNSLXML.h,v $
//Revision 1.6  2004/09/14 13:49:59  jbarbier
//- files added in sample/ to run run the main_siconos test program
//
//- all the platform can now be saved in an XML file when it is created manually
//
//Revision 1.5  2004/07/29 14:25:45  jbarbier

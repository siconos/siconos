/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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

/** \classRelayNSLawXML
 *   \brief This class manages RelayNSLaw data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.2.0.
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
    return  SiconosDOMTreeTools::getContentValue<double>(this->CNode);
  }

  /** \fn double getD()
   *   \brief Return the D of a relayNSLaw
   *   \return The D double value of the relayNSLaw
   */
  inline double getD()
  {
    return  SiconosDOMTreeTools::getContentValue<double>(this->DNode);
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

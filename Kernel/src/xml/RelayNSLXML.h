/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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

#ifndef __RelayNSLawXML__
#define __RelayNSLawXML__

#include "NonSmoothLawXML.h"

const std::string RNSL_C = "c";
const std::string RNSL_D = "d";

/** XML management for RelayNSL
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.1.
 *   \date 05/14/2004
 *
 *
 *
 * RelayNSLawXML allows to manage data of a RelayNSLaw DOM tree.
 */
class RelayNSLXML : public NonSmoothLawXML
{
public:
  RelayNSLXML();

  /** Build a RelayNSLXML object from a DOM tree describing a Law with Relay type
  *   \param relayNSLawNode : the relayNSLaw DOM tree
  *   \exception XMLException : if a property of the Relay NS Law lacks in the DOM tree
  */
  RelayNSLXML(xmlNode * relayNSLawNode);


  /** Return the C of a relayNSLaw
  *   \return The C double value of the relayNSLaw
  */
  inline double getC()
  {
    return  SiconosDOMTreeTools::getContentValue<double>(CNode);
  }

  /** Return the D of a relayNSLaw
  *   \return The D double value of the relayNSLaw
  */
  inline double getD()
  {
    return  SiconosDOMTreeTools::getContentValue<double>(DNode);
  }

  /** set the C of a relayNSLaw
  *   \param double : The C double value of the relayNSLaw
  */
  inline void setC(double c)
  {
    if (CNode == NULL)
    {
      CNode = SiconosDOMTreeTools::createDoubleNode(rootNode, RNSL_C, c);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(CNode, c);
  }

  /** set the D of a relayNSLaw
  *   \param double : The D double value of the relayNSLaw
  */
  inline void setD(double d)
  {
    if (DNode == NULL)
    {
      DNode = SiconosDOMTreeTools::createDoubleNode(rootNode, RNSL_D, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(DNode, d);
  }


private:
  xmlNode * CNode;
  xmlNode * DNode;

};

#endif

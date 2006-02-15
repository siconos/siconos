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
/** \class MoreauXML
 *   \brief This class manages Moreau data part
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 05/17/2004
 *
 *
 * MoreauXML allows to manage data of a Moreau DOM tree.
 */

#ifndef __MOREAUXML__
#define __MOREAUXML__

#include "OneStepIntegratorXML.h"

const std::string MOREAU_R = "r";
const std::string MOREAU_W = "W";
const std::string MOREAU_THETA = "Theta";

class MoreauXML : public OneStepIntegratorXML
{
public:
  MoreauXML();

  /** \fn MoreauXML(xmlNode * MoreauNode)
   *   \brief Build a MoreauXML object from a DOM tree describing Moreau OneStepIntegrator
   *   \param xmlNode * MoreauNode : the Moreau DOM tree
   *   \param map<int, bool> definedDSNumbers : to know if DS numbers are not used by another OneStepIntegrator
   *   \exception XMLException : if the W property of the Moreau lacks in the DOM tree
   */
  MoreauXML(xmlNode * MoreauNode, std::map<int, bool> definedDSNumbers);

  // Destructor
  ~MoreauXML();


  /** \fn bool hasW()
   *  \brief return true if wNode is defined
   *  \return true if wNode is defined
   */
  inline bool hasW()
  {
    return (WNode != NULL);
  }

  /** \fn SiconosMatrix getW()
   *   \brief Return the w of the OneStepIntegratorXML
   *   \return SiconosMatrix : the w of the OneStepIntegratorXML
   */
  inline SiconosMatrix getW()
  {
    return  SiconosDOMTreeTools::getSiconosMatrixValue(WNode);
  }

  /** \fn void setW(SiconosMatrix *m)
   *   \brief allows to save the w of the OneStepIntegratorXML
   *   \param SiconosMatrix* : the w to save
   */
  inline void setW(SiconosMatrix *m)
  {
    if (hasW() == false)
    {
      WNode = SiconosDOMTreeTools::createMatrixNode(rootIntegratorXMLNode, MOREAU_W, *m);
    }
    else SiconosDOMTreeTools::setSiconosMatrixNodeValue(WNode, *m);
  }

  /** \fn bool hasTheta()
   *  \brief return true if ThetaNode is defined
   *  \return true if ThetaNode is defined
   */
  inline bool hasTheta()
  {
    return (ThetaNode != NULL);
  }

  /** \fn SiconosMatrix getTheta()
   *   \brief Return the theta of the OneStepIntegratorXML
   *   \return SiconosMatrix : the theta of the OneStepIntegratorXML
   */
  inline const double getTheta() const
  {
    return  SiconosDOMTreeTools::getDoubleContentValue(ThetaNode);
  }

  /** \fn void setTheta(double t)
   *   \brief allows to save  Theta of the OneStepIntegratorXML
   *   \param double t : the Theta to save
   */
  inline void setTheta(const double& t)
  {
    if (hasTheta() == false)
    {
      ThetaNode = SiconosDOMTreeTools::createDoubleNode(rootIntegratorXMLNode, MOREAU_THETA, t);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(ThetaNode, t);
  }
private:

  //Nodes
  xmlNode * WNode;
  xmlNode * ThetaNode;


};


#endif

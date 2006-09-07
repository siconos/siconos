/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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
 *   \brief read Moreau-related nodes in DOM tree
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.3.0.
 *   \date 05/17/2004
 *
 */

#ifndef __MOREAUXML__
#define __MOREAUXML__

#include "OneStepIntegratorXML.h"

class MoreauXML: public OneStepIntegratorXML
{
private:

  /** theta list node */
  xmlNode * thetaNode;

  /** W list node */
  xmlNode * WNode;

public:
  MoreauXML();

  /** \fn MoreauXML(xmlNodePtr MoreauNode)
   *   \brief Build a MoreauXML object from a DOM tree describing Moreau OneStepIntegrator
   *   \param xmlNode * MoreauNode : the Moreau DOM tree
   *   \exception XMLException : if the W property of the Moreau lacks in the DOM tree
   */
  MoreauXML(xmlNode * MoreauNode);

  /** \fn ~MoreauXML()
   * Destructor
   */
  ~MoreauXML();

  /** \fn bool hasWList()
   *  \brief return true if wNode is defined
   *  \return true if wNode is defined
   */
  inline bool hasWList()
  {
    return (WNode != NULL);
  }

  /** \fn SimpleMatrix getW()
   *   \brief Return the w of the OneStepIntegratorXML
   *   \return SimpleMatrix : the w of the OneStepIntegratorXML
   */
  // inline SimpleMatrix getW()
  // {
  //  return  SiconosDOMTreeTools::getSiconosMatrixValue(WNode);
  //}

  /** \fn void setW(SiconosMatrix *m)
   *   \brief allows to save the w of the OneStepIntegratorXML
   *   \param SiconosMatrix* : the w to save
   */
  //  inline void setW(SiconosMatrix *m)
  //{
  //  if( hasW() == false )
  //{
  //  WNode = SiconosDOMTreeTools::createMatrixNode(rootIntegratorXMLNode, MOREAU_W, *m);
  //}
  //  else SiconosDOMTreeTools::setSiconosMatrixNodeValue(WNode, *m);
  //}

  /** \fn bool hasThetaList()
   *  \brief return true if ThetaNode is defined
   *  \return true if ThetaNode is defined
   */
  inline bool hasThetaList() const
  {
    return (thetaNode != NULL);
  }

  /** \fn void getTheta(vector<double>&)
   *   \brief fill a vector<double> with given theta values
   *   \param: in-out vector<double>
   */
  void getTheta(std::vector<double>&) const;

  /** \fn void setTheta(const vector<double&>)
   *   \brief save  theta values in xml ouput file or DOMtree
   *   \param vector<double>
   */
  void setTheta(const std::vector<double>&);

  /** \fn const double getSingleTheta() const
   *   \brief get value of attribute all in theta node -> ie if one single value for all theta is given
   *   \return a double
   */
  const double getSingleTheta() const;

  /** \fn bool hasAllTheta()
   *  \brief attribute of the theta tag - all = val if all theta have the same value, whatever the ds is.
   *  \return a bool
   */
  inline bool hasAllTheta() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(thetaNode, ALL_ATTRIBUTE);
  }

};
#endif

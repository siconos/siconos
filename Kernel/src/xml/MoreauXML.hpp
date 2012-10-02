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

#ifndef __MOREAUXML__
#define __MOREAUXML__

#include "SiconosPointers.hpp"
#include "OneStepIntegratorXML.hpp"

/** XML management for Moreau
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/17/2004
 *
 */
class MoreauXML: public OneStepIntegratorXML
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MoreauXML);


  /** theta list node */
  xmlNode * thetaNode;

  /** W list node */
  xmlNode * WNode;

public:
  MoreauXML();

  /** Build a MoreauXML object from a DOM tree describing Moreau OneStepIntegrator
  *   \param xmlNode * MoreauNode : the Moreau DOM tree
  *   \exception XMLException : if the W property of the Moreau lacks in the DOM tree
  */
  MoreauXML(xmlNode * MoreauNode);

  /** Destructor */
  ~MoreauXML();

  /** return true if wNode is defined
   *  \return true if wNode is defined
  */
  inline bool hasWList()
  {
    return (WNode);
  };

  /** Return the w of the OneStepIntegratorXML
  *   \return SimpleMatrix : the w of the OneStepIntegratorXML
  */
  // inline SimpleMatrix getW()
  // {
  //  return  SiconosDOMTreeTools::getSiconosMatrixValue(WNode);
  //}

  /** allows to save the w of the OneStepIntegratorXML
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

  /** return true if ThetaNode is defined
  *  \return true if ThetaNode is defined
  */
  inline bool hasThetaList() const
  {
    return (thetaNode);
  }

  /** fill a vector<double> with given theta values
  *   \param: in-out vector<double>
  */
  void getTheta(std::vector<double>&) const;

  /** save  theta values in xml ouput file or DOMtree
  *   \param vector<double>
  */
  void setTheta(const std::vector<double>&);

  /** get value of attribute all in theta node -> ie if one single value for all theta is given
  *   \return a double
  */
  double getSingleTheta() const;

  /** attribute of the theta tag - all = val if all theta have the same value, whatever the ds is.
  *  \return a bool
  */
  inline bool hasAllTheta() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(thetaNode, ALL_ATTRIBUTE);
  }

};

DEFINE_SPTR(MoreauXML)
#endif

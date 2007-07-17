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

/*! \file FirstOrderRXML.h
  \brief XML data reading for first order non linear relations.
*/

#ifndef FirstOrderRXML_H
#define FirstOrderRXML_H

#include "RelationXML.h"
#include "SimpleMatrix.h"

class FirstOrderR;

/** XML management for FirstOrderR
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.1.
 *   \date 02/21/2007
 *
 */
class FirstOrderRXML : public RelationXML
{

protected:

  /** Plug-in name for input function g*/
  xmlNodePtr gNode;

  /** Plug-in name for output function h. */
  xmlNodePtr hNode;

  /** plug-in names for gradients of g.*/
  std::vector<xmlNodePtr> jacobianGNode;

  /** plug-in names for gradients of h.*/
  std::vector<xmlNodePtr> jacobianHNode;

  /**Default construcor. */
  FirstOrderRXML();

public:

  /** Basic data constructor
      \param the node that points to the FirstOrderR data.
  */
  FirstOrderRXML(xmlNodePtr);

  /** Destructor */
  virtual ~FirstOrderRXML();

  // ================== g ==================

  /** To set the name of the plug-in function used to compute g.
   *   \param a string.
   */
  void setGPlugin(const std::string&);

  /** To get the name of the plug-in function used to compute g.
   *   \return a string.
   */
  std::string getGPlugin() const ;

  /** Return true if g is calculated from a plugin
   * \return a bool. */
  inline bool isGPlugin() const
  {
    return xmlHasProp((xmlNodePtr)gNode, (xmlChar *)"plugin");
  }

  /** Tests if gNode is defined.
   *  \return a bool.
   */
  inline bool hasG() const
  {
    return (gNode != NULL);
  }
  // ================== h ==================

  /** To set the name of the plug-in function used to compute h.
   *   \param a string.
   */
  void setHPlugin(const std::string&);

  /** To get the name of the plug-in function used to compute h.
   *   \return a string.
   */
  std::string getHPlugin() const ;

  /** Return true if h is calculated from a plugin
   * \return a bool. */
  inline bool isHPlugin() const
  {
    return xmlHasProp((xmlNodePtr)hNode, (xmlChar *)"plugin");
  }

  /** Tests if hNode is defined.
   *  \return a bool.
   */
  inline bool hasH() const
  {
    return (hNode != NULL);
  }

  // ================== jacobianG ==================

  /** To set the name of the plug-in function used to compute jacobianG.
   *   \param a string.
   *   \param index of jacobian, default = 0
   */
  void setJacobianGPlugin(const std::string&, unsigned int = 0);

  /** To get the name of the plug-in function used to compute jacobianG.
   * \param index of jacobian, default = 0
   *   \return a string.
   */
  std::string getJacobianGPlugin(unsigned int = 0) const ;

  /** Return true if jacobianG is calculated from a plugin
      \param index of jacobian, default = 0
      * \return a bool. */
  inline bool isJacobianGPlugin(unsigned int i = 0) const
  {
    return xmlHasProp((xmlNodePtr)(jacobianGNode[i]), (xmlChar *)"matrixPlugin");
  }

  /** Tests if gNode is defined.
   *  \return a bool.
   */
  inline bool hasJacobianG() const
  {
    return (!jacobianGNode.empty());
  }

  /** Return the jacobianG[i] matrix.
   * \param an int: the required index for jacobianG.
   * \return a SimpleMatrix
   */
  SimpleMatrix getJacobianGMatrix(unsigned int = 0) const;

  // ================== jacobianH ==================

  /** To set the name of the plug-in function used to compute jacobianH.
   *   \param a string.
   * \param index of jacobian
   */
  void setJacobianHPlugin(const std::string&, unsigned int);

  /** To get the name of the plug-in function used to compute jacobianH.
      \param index of jacobian
      *   \return a string.
      */
  std::string getJacobianHPlugin(unsigned int) const ;

  /** Return true if jacobianH[i] is calculated from a plugin
      \param index of jacobian
      \return a bool. */
  inline bool isJacobianHPlugin(unsigned int i) const
  {
    return xmlHasProp((xmlNodePtr)(jacobianHNode[i]), (xmlChar *)"matrixPlugin");
  }

  /** Tests if hNode is defined.
   *  \return a bool.
   */
  inline bool hasJacobianH() const
  {
    return (!jacobianHNode.empty());
  }

  /** Return the jacobianH[i] matrix.
   * \param an int: the required index for jacobianH.
   * \return a SimpleMatrix
   */
  SimpleMatrix getJacobianHMatrix(unsigned int) const;

};

#endif

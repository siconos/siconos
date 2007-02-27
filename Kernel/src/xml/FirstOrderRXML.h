/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

class FirstOrderR;

/** XML management for FirstOrderR
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.1.
 *   \date 02/21/2007
 *
 */
class FirstOrderRXML : public RelationXML
{

protected:

  /** Plug-in name for input */
  xmlNodePtr computeInputNode;

  /** Plug-in name for output. */
  xmlNodePtr computeOutputNode;

  /**Default construcor. */
  FirstOrderRXML();

public:

  /** Basic data constructor
   \param the node that points to the FirstOrderR data.
  */
  FirstOrderRXML(xmlNodePtr);

  /** Destructor */
  virtual ~FirstOrderRXML();

  /** To get the name of the plug-in function used to compute input.
   *   \return a string.
   */
  virtual std::string getComputeInputPlugin() const ;

  /** To get the name of the plug-in function used to compute output.
   *   \return a string.
   */
  virtual std::string getComputeOutputPlugin() const ;

  /** To set the name of the plug-in function used to compute input.
   *   \param a string.
   */
  void setComputeInputPlugin(const std::string&);

  /** To set the name of the plug-in function used to compute outut.
   *   \param a string.
   */
  void setComputeOutputPlugin(const std::string&);

  /** Return true if computeInput is calculated from a plugin
   * \return a bool. */
  inline bool isComputeInputPlugin() const
  {
    return xmlHasProp((xmlNodePtr)computeInputNode, (xmlChar *)"plugin");
  }

  /** Return true if computeOutput is calculated from a plugin
   * \return a bool. */
  inline bool isComputeOutputPlugin() const
  {
    return xmlHasProp((xmlNodePtr)computeOutputNode, (xmlChar *)"plugin");
  }

  /** Tests if computeInputNode is defined.
   *  \return a bool.
   */
  inline bool hasComputeInput() const
  {
    return (computeInputNode != NULL);
  }

  /** Tests if computeOutputNode is defined.
  *  \return a bool.
  */
  inline bool hasComputeOutput() const
  {
    return (computeOutputNode != NULL);
  }

};

#endif

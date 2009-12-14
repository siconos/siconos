/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
/*! \file TimeDiscretisationXML.h

  \brief: xml management for TimeDiscretisation
*/

#ifndef __TIMEDISCRETISATIONXML__
#define __TIMEDISCRETISATIONXML__


#include "SiconosDOMTreeTools.hpp"

class TimeDiscretisation;

const std::string  TD_H = "h";
const std::string  TD_N = "N";
const std::string  TD_TK = "tk";

/** Data type used to save time instants of the discretisation */
typedef std::vector<double> TkVector;

/** XML management for Time Discretisation
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/17/2004
 *
 */
class TimeDiscretisationXML
{
private:

  //Nodes
  xmlNodePtr rootNode;
  xmlNodePtr hNode;
  xmlNodePtr NNode;
  xmlNodePtr tkNode;

  /** Default constructor */
  TimeDiscretisationXML();

public:
  /** Build a TimeDiscretisationXML object from a DOM tree describing a TimeDiscretisation
   *   \param timeDiscretisationNode : the TimeDiscretisation DOM tree
   */
  TimeDiscretisationXML(xmlNodePtr);


  /** Gets the rootNode of the TimeDiscretisationXML
   *   \return xmlNode* : the rootNode
   */
  inline xmlNodePtr getRootNode()
  {
    return rootNode;
  }

  /** sets the rootNode of the TimeDiscretisationXML
   *   \param a xmlNodePtr
   */
  inline void setRootNode(xmlNodePtr node)
  {
    rootNode = node;
  }

  /** Returns the chosen value in xml file for the time step
   *   \return a double
   */
  inline double geth() const
  {
    return  SiconosDOMTreeTools::getContentValue<double>(hNode);
  }

  /** Set time step value in xml output
   *  \param a double
   */
  inline void setH(double d)
  {
    if (!hasH())
    {
      hNode = SiconosDOMTreeTools::createDoubleNode(rootNode, TD_H, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(hNode, d);
  }

  /** returns the xml value of the chosen number of time steps
      \return an int
  */
  inline int getN() const
  {
    return SiconosDOMTreeTools::getContentValue<int>(NNode);
  }

  /** sets in xml file the value for the number of time steps
   * \param an int
   */
  inline void setN(int i)
  {
    if (!hasN())
    {
      NNode = SiconosDOMTreeTools::createIntegerNode(rootNode, TD_N, i);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(NNode, i);
  }

  /** Returns the vector tk given in the xml input file
   * \param[in,out] the TkVector to be set
   */
  inline void getTk(TkVector tk) const
  {
    return SiconosDOMTreeTools::getVector<double>(tkNode, tk);
  }

  /** set the tk Node value of the TimeDiscretisation
   *   \param The tk SiconosVector of the TimeDiscretisation
   */
  inline void setTkNode(const TkVector& v)
  {
    /*     if( !hasTk()) */
    /*       { */
    /*  tkNode = SiconosDOMTreeTools::createVectorNode(rootNode, TD_TK, v); */
    /*       } */
    //    else SiconosDOMTreeTools::setSiconosVectorNodeValue(tkNode, v);
    XMLException::selfThrow("TimeDiscretisationXML::setTkNode: set of vector<double> not yet implemented in SiconosDOMTreeTools.hpp");
  }

  /** returns true if hNode is defined
   *  \return true if hNode is defined
   */
  inline bool hasH() const
  {
    return (hNode);
  }

  /** returns true if NNode is defined
   *  \return true if NNode is defined
   */
  inline bool hasN() const
  {
    return (NNode);
  }

  /** returns true if tkNode is defined
   *  \return true if tkNode is defined
   */
  inline bool hasTk() const
  {
    return (tkNode);
  }

  /** makes the operations to create the TimeDiscretisation of the SimulationXML
   *   \param xmlNode* : the root node for the TimeDiscretisationXML
   *   \param NSDS* : the TimeDiscretisation of this TimeDiscretisationXML
   */
  void updateTimeDiscretisationXML(xmlNode*, TimeDiscretisation*);

};

#endif

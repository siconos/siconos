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
/*! \file TimeDiscretisationXML.hpp

  \brief: xml management for TimeDiscretisation
*/

#ifndef __TIMEDISCRETISATIONXML__
#define __TIMEDISCRETISATIONXML__


#include "SiconosDOMTreeTools.hpp"
#include "SiconosPointers.hpp"

//class TimeDiscretisation;

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
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(TimeDiscretisationXML);

  //Nodes
  xmlNodePtr _rootNode;
  xmlNodePtr _hNode;
  xmlNodePtr _NNode;
  xmlNodePtr _tkNode;

  /** Default constructor */
  TimeDiscretisationXML() {};

public:
  /** Build a TimeDiscretisationXML object from a DOM tree describing a TimeDiscretisation
   *   \param timeDiscretisationNode the TimeDiscretisation DOM tree
   */
  TimeDiscretisationXML(xmlNodePtr timeDiscretisationNode);

  /** Copy constructor
   * \param tdXML the TimeDiscretisationXML to copy
   */
  TimeDiscretisationXML(const TimeDiscretisationXML & tdXML);

  /** Destructor */
  virtual ~TimeDiscretisationXML() {};

  /** Gets the rootNode of the TimeDiscretisationXML
   *   \return a xmlNodePtr, the rootNode
   */
  inline xmlNodePtr getRootNode() const
  {
    return _rootNode;
  }

  /** sets the rootNode of the TimeDiscretisationXML
   *   \param node a xmlNodePtr
   */
  inline void setRootNode(xmlNodePtr node)
  {
    _rootNode = node;
  }

  /** Returns the chosen value in xml file for the time step
   *   \return a double
   */
  inline double geth() const
  {
    return SiconosDOMTreeTools::getContentValue<double>(_hNode);
  }

  /** Set time step value in xml output
   *  \param d a double
   */
  inline void setH(double d)
  {
    if (!hasH())
    {
      _hNode = SiconosDOMTreeTools::createDoubleNode(_rootNode, TD_H, d);
    }
    else SiconosDOMTreeTools::setDoubleContentValue(_hNode, d);
  }

  /** returns the xml value of the chosen number of time steps
      \return an int
  */
  inline int getN() const
  {
    return SiconosDOMTreeTools::getContentValue<int>(_NNode);
  }

  /** sets in xml file the value for the number of time steps
   * \param i an int
   */
  inline void setN(int i)
  {
    if (!hasN())
    {
      _NNode = SiconosDOMTreeTools::createIntegerNode(_rootNode, TD_N, i);
    }
    else SiconosDOMTreeTools::setIntegerContentValue(_NNode, i);
  }

  /** Returns the vector tk given in the xml input file
   * \param[in,out] tk the TkVector that will contain _tkNode
   */
  inline void getTk(TkVector &tk) const
  {
    SiconosDOMTreeTools::getVector<double>(_tkNode, tk);
  }

  /** set the tk Node value of the TimeDiscretisation
   *   \param v The tk SiconosVector of the TimeDiscretisation
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
    return (_hNode);
  };

  /** returns true if NNode is defined
   *  \return true if NNode is defined
   */
  inline bool hasN() const
  {
    return (_NNode);
  };

  /** returns true if tkNode is defined
   *  \return true if tkNode is defined
   */
  inline bool hasTk() const
  {
    return (_tkNode);
  };

  /** makes the operations to create the TimeDiscretisation of the SimulationXML
   *   \param node the root node for the TimeDiscretisationXML
   *   \param td the TimeDiscretisation of this TimeDiscretisationXML
   */
  void updateTimeDiscretisationXML(xmlNodePtr node, SP::TimeDiscretisation td);

};

#endif

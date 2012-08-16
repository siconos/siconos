/* Siconos-Kernel, Copyright INRIA 2005-2011.
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

/*! \file DynamicalSystemXML.hpp

 */

#ifndef __DynamicalSystemXML__
#define __DynamicalSystemXML__

#include "SiconosDOMTreeTools.hpp"
#include "DynamicalSystemTypes.hpp"

#include "SiconosVisitor.hpp"

/** XML management for DynamicalSystem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 04/04/2004
 *
 * Reading/Writing of xml data for the DynamicalSystem class and its derived classes.
 *
 */
class DynamicalSystemXML
{

protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(DynamicalSystemXML);


  xmlNodePtr rootNode;/**< Dynamical System Node  */
  xmlNodePtr stepsInMemoryNode; /**< size of memory */
  xmlNodePtr zNode; /**< z in \f$ M \dot x = f(x,t,z) \f$ */

  /** Default constructor */
  DynamicalSystemXML(): rootNode(NULL), stepsInMemoryNode(NULL), zNode(NULL) {};

public:

  /** Build a DynamicalSystemXML object from the DynamicalSystem node of the xml DOMtree
   *   \param an xmlNodePtr DynamicalSystemNode
   *   \param bool isBVP : if true, NonSmoothDynamicalSystem is a Boundary value problem
   */
  DynamicalSystemXML(xmlNodePtr DynamicalSystemNode, bool);

  /** Destructor */
  virtual inline ~DynamicalSystemXML() {};

  /** Return the number of the DynamicalSystem
      \return an integer
  */
  inline int number() const
  {
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, NUMBER_ATTRIBUTE);
  }

  /** Return the type of the DynamicalSystem
   *   \return a string
   */
  Type::Siconos getType() const;

  /** Returns the z vector, discret state of the DynamicalSystem
   *  \return SiconosVector
   */
  inline const SiconosVector getz() const
  {
    return  SiconosDOMTreeTools::getSiconosVectorValue(zNode);
  }

  /** save z of the DynamicalSystem
   *   \param a SiconosVector
   */
  inline void setz(const SiconosVector& v)
  {
    if (!hasz())
      zNode = SiconosDOMTreeTools::createVectorNode(rootNode, "z", v);
    else SiconosDOMTreeTools::setSiconosVectorNodeValue(zNode, v);
  }

  /** Returns the steps in memory for the DynamicalSystemXML
   *   \return The integer number of steps in memory for the DynamicalSystemXML
   */
  inline unsigned int getStepsInMemory() const
  {
    return  SiconosDOMTreeTools::getContentValue<unsigned int>(stepsInMemoryNode);
  }

  /** to save the steps in memory for the DynamicalSystemXML
   *   \param an unsigned int
   */
  void setStepsInMemory(const unsigned int&);

  /** returns true if stepsInMemoryNode is defined
   *  \return a bool
   */
  inline bool hasStepsInMemory() const
  {
    return (stepsInMemoryNode);
  }

  /** returns true if zNode is defined
   *  \return a bool
   */
  inline bool hasz() const
  {
    return (zNode);
  }

};

#endif

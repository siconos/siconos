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

/*! \file RelationXML.hpp
  \brief XML data reading for Relation.
*/

#ifndef RelationXML_H
#define RelationXML_H

#include "SiconosDOMTreeTools.hpp"
#include "RelationNamespace.hpp"
#include "SimpleVector.hpp"

using RELATION::TYPES;
using RELATION::SUBTYPES;

/** XML management for Relation
   \author SICONOS Development Team - copyright INRIA
    \version 3.0.0.
    \date 04/13/2004

    The tag for relation must belong to Interaction tag.

    It is specified with a type (FirstOrder or Lagrangian), a sub-type and the different operators.

    Example:

    \code
    <LagrangianRelation type="Scleronomous">
    ...
    </LagrangianRelation>
    \endcode

    The full content of each tag depends on the sub-type of the relation.


 */
class RelationXML
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(RelationXML);


private:
  RelationXML(const RelationXML&);

protected:

  /** Relation node.*/
  xmlNodePtr rootNode;

  /** Node related to output function h. */
  xmlNodePtr hNode;

  /** Node related to input function g*/
  xmlNodePtr gNode;

  /** Node related to hDot.*/
  xmlNodePtr hDotNode;

  /** Nodes related to jacobian functions of h.*/
  std::vector<xmlNodePtr> jacobianHNode;

  /** Nodes related to jacobian functions of g.*/
  std::vector<xmlNodePtr> jacobianGNode;

  /**Default constructor. */
  RelationXML(): rootNode(NULL), hNode(NULL), gNode(NULL), hDotNode(NULL) {};

public:

  /**Basic construcor.
     \param a pointer to the relation xml node.
  */
  RelationXML(xmlNodePtr);

  /** Destructor */
  virtual ~RelationXML() {};

  /** Returns the type of the RelationXML
   */
  RELATION::TYPES getType() const ;

  /** Returns the sub-type of the Relation
   *  \return a string.
   */
  RELATION::SUBTYPES getSubType() const;

  /** Returns the node of the RelationXML
   *   \return an xmlNodePtr.
   */
  inline xmlNodePtr getRootNode()const
  {
    return rootNode;
  }

  /** Gets the number of jacobian of h (index 0) or g (index 1) defined in XML file
      \param index value to choose among h and g
      \return an unsigned int.
  */
  unsigned int getNumberOfJacobians(unsigned int) const;

  // ================== h ==================

  /** To set the name of the plug-in function used to compute h.
   *   \param a string.
   */
  void setHPlugin(const std::string&);

  /** To get the name of the plug-in function used to compute h.
   *   \return a string.
   */
  std::string gethPlugin() const ;

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
    return (hNode);
  }

  /** checks if hDot is defined in the DOM tree
   *  \return a bool
   */
  inline bool hasHDot() const
  {
    return (hDotNode);
  };

  /** Return true if hDot is calculated with a plugin.
   *  \return a bool
   */
  inline bool isHDotPlugin() const
  {
    return xmlHasProp(hDotNode, (xmlChar*)"plugin");
  };

  /** Return hDot plug-in name, if it exists
   *   \return a string
   */
  inline std::string gethDotPlugin() const
  {
    if (!isHDotPlugin())
      XMLException::selfThrow("RelationRXML - gethDotPlugin : hDot is not plugged.");
    return  SiconosDOMTreeTools::getStringAttributeValue(hDotNode, "plugin");
  }

  /** Return hDot vector.
   *   \return a SimpleVector.
   */
  inline const SimpleVector gethDotVector() const
  {
    if (isHDotPlugin())
      XMLException::selfThrow("RelationRXML - gethDotVector : hDot is plugged, not set with a vector.");
    return  SiconosDOMTreeTools::getSiconosVectorValue(hDotNode);
  }

  /** to save the hDot plug-in name
   *   \param a string
   */
  void setHDotPlugin(const std::string&);

  /** to save the hDot vector.
   *   \param a SiconosVector
   */
  void setHDotVector(const SiconosVector&);

  // ================== g ==================

  /** To set the name of the plug-in function used to compute g.
   *   \param a string.
   */
  void setGPlugin(const std::string&);

  /** To get the name of the plug-in function used to compute g.
   *   \return a string.
   */
  std::string getgPlugin() const ;

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
    return (gNode);
  }

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

  /** To create the Relation of the Interaction. Useless?
   *   \param xmlNodePtr : the root node of the RelationXML
   *   \param Relation* : the Relation of this RelationXML
   */
  inline void updateRelationXML(xmlNodePtr node, Relation*)
  {
    rootNode = node;
  }

  /** xml utility: used to read a jacobian matrix or plugged matrix in xml file
      \param T object (must be a smart pointer): the jacobian to be filled in
      \param RelationXML, the xml pointer
      \param i, index of required jacobian
  */
  template <class T, class SPT> void readJacobianXML(SPT& jacobian, SP::RelationXML LRxml, unsigned int i)
  {
    if (isJacobianHPlugin(i))
      jacobian.reset(new T(getJacobianHPlugin(i)));
    else
      jacobian.reset(new T(getJacobianHMatrix(i)));
  }

};

#endif

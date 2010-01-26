/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
/*! \file OneStepNSProblemXML.h
  XML management for One Step Non Smooth Problem interface
*/

#ifndef __OneStepNSProblemXML__
#define __OneStepNSProblemXML__

#include "SiconosPointers.hpp"
#include "SiconosDOMTreeTools.hpp"


class OneStepNSProblem;
class NonSmoothSolverXML;
/** XML management for OneStepNSProblem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 3.0.0.
 *   \date 05/14/2004
 *
 * XML management for classes derived from OneStepNSProblem
 *
 *
 * XML tag for OneStepNSProblem (example for a LCP problem)

 \code
 <OneStepNSProblems_List>
 <LCP StorageType="0" Id="myLCP">
 <size>52</size>
 <Interactions_Concerned vectorSize='3'>1 4 5</Interactions_Concerned>
 <NonSmoothSolver Name="Lemke">
 <iparam  vectorSize='3'>1 3 4</iparam>
 <dparam  vectorSize='3'>1.3 1e-17 1e-12</dparam>
 </NonSmoothSolver>
 </LCP>
 </OneStepNSProblems_List>
 \endcode
 Interactions_Concerned is optional (default: all). \n
 NonSmoothSolver is optional (default: read XXX.opt file in Numerics, XXX being the formulation type (LCP ...) \n
 size optional (computed during preCompute() in OneStepNSProblem, according to simulation info).

 => minimum input in XML:
 \code
 <OneStepNSProblems_List>
 <LCP></LCP>
 </OneStepNSProblems_List>
 \endcode

 Remark: it is possible to define several OSNS for one Simulation like this:
 \code
 <OneStepNSProblems_List>
 <LCP>...</LCP>
 <FrictionContact>...</FrictionContact>
 </OneStepNSProblems_List>
 \endcode

*/
class OneStepNSProblemXML
{
protected:

  /** root node, name of the OSNS (LCP, FrictionContact ...) */
  xmlNodePtr rootNode;
  /** dimension of the problem (tag: size) */
  xmlNodePtr dimNode;
  /** node used to list the interactions involved in the OSNSS (tag Interactions_Concerned) */
  xmlNodePtr interactionsConcernedNode;
  /** node used to define solver (tag: NumericsSolverName) - Child of rootNode */
  xmlNodePtr numericsSolverNameNode;
  /** solverXML object */


public:

  /** Default constructor */
  inline OneStepNSProblemXML(): rootNode(NULL), dimNode(NULL), interactionsConcernedNode(NULL),
    numericsSolverNameNode(NULL) {};

  /** Build a OneStepNSProblemXML object from a DOM tree describing a OneStepNSProblem
   *   \param OneStepNSProblemNode : the OneStepNSProblem DOM tree
   *   \exception XMLException : if a property of the OneStepNSProblemXML lacks in the DOM tree
   */
  OneStepNSProblemXML(xmlNodePtr);

  /** Destructor
   */
  virtual ~OneStepNSProblemXML() {};

  /** Return the type of the OneStepNSProblem
   *   \return a string
   */
  inline std::string  getNSProblemType() const
  {
    //std::string  type((char*)problemTypeNode->name);
    return (char*)rootNode->name;
  }

  /** Return the dimension of the OneStepNSProblem
   *   \return an integer
   */
  inline int getDimNSProblem() const
  {
    return SiconosDOMTreeTools::getContentValue<int>(dimNode);
  }

  /** set dimension of the OneStepNSProblem
   *   \param an integer
   */
  void setDimNSProblem(int);

  /** returns true if dimNode is defined
   *  \return a bool
   */
  inline bool hasDim() const
  {
    return (dimNode);
  }

  /** All is an attribute of the Interactions_Concerned tag
   *  \return bool : true if attribute all is defined
   */
  bool hasAllInteractions() const;

  /** to set the attribute "all" of the Interaction_concerned tag
   *   \param bool : the value to assign to the attribute
   */
  void setAllInteractions(const bool&);

  /** Checks if attribute "storageType" is given in formalization tag
   *  \return a bool
   */
  inline bool hasStorageType() const
  {
    return (SiconosDOMTreeTools::hasAttributeValue(rootNode, "StorageType"));
  }

  /** Returns the value of attribute "storageType" in formalization tag
   *  \return an integer
   */
  inline int getStorageType() const
  {
    if (!hasStorageType())
      XMLException::selfThrow("OneStepNSProblemXML::getStorageType - Attribute named storageType does not exists in tag formalization.");
    return SiconosDOMTreeTools::getAttributeValue<int>(rootNode, "StorageType");
  }


  /** makes the operations to create a OneStepNSProblemXML to the SimulationXML
   *   \param xmlNodePtr : the root node of the OneStepNSProblemXML
   *   \param SP::OneStepNSProblem : the OneStepNSProblem of this OneStepNSProblemXML
   */
  void updateOneStepNSProblemXML(xmlNodePtr , SP::OneStepNSProblem);

  /** Return the id of the OneStepNSProblem (attribute of the root node)
   *   \return a string
   */
  inline std::string getId() const
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, "Id");
  }

  /** Return true if the id attribute of the rootNode is present.
     *   \return a bool   */
  inline bool hasId() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(rootNode, "Id");
  }

  /** Return the  NumericsSolverName of the OneStepNSProblem (attribute of the root node)
   *   \return a string
   */
  inline std::string getNumericsSolverName() const
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, "NumericsSolverName");
  }

  /** Return true if the NumericsSolverName attribute of the rootNode is present.
  *   \return a bool   */
  inline bool hasNumericsSolverName() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(rootNode, "NumericsSolverName");
  }

  /** return a vector<int> of the number of the interactions related to the OSNS
   *  \param in-out vector<int>
   */
  void getInteractionsNumbers(std::vector<int>&);
};

#endif

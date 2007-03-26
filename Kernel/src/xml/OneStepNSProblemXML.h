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
/*! \file
*/

#ifndef __OneStepNSProblemXML__
#define __OneStepNSProblemXML__

#include "SiconosDOMTreeTools.h"

class OneStepNSProblem;
class SolverXML;
/** XML management for OneStepNSProblem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.0.1.
 *   \date 05/14/2004
 *
 * OneStepNSProblemXML class id dedicated data management for OneStepNSProblem DOM.
 *
 * There are two steps in OneStepNSProblem xml loading:
 *  - Problem formalisation part ( LCP, FrictionContact ...) => done in derived class constructor
 *    This is the way the problem is written in Siconos Kernel. Two each formalisation corresponds
 *    a OneStepNSProblem derived class.
 *
 *  - Solver data loading => done in OneStepNS top class constructor.
 *    First data is the solving formalisation, ie how to problem is written in order to be solved:
 *    LcpSolving, PrimalRelaySolving etc ...
 *    This formalisation can be different from the one used in Siconos Kernel, ie from the OneStepNSProblem
 *    type.
 *    Second data is the solver type (NLGS ...) and its parameters.
 *
 */
class OneStepNSProblemXML
{
protected:

  /** root node named "OneStepNSProblem" */
  xmlNode* rootNode;

  /** root node for formalisation type of the problem in Siconos/Kernel (LCP, QP ...)
   * This corresponds to the OneStepNSProblem type - Child of rootNode */
  xmlNode* problemTypeNode;
  /** dimension of the problem - Child of problemTypeNode */
  xmlNode* dimNode;
  /** node named interaction_Concerned  - Child of problemTypeNode */
  xmlNode* interactionConcernedNode;
  /** Interactions list (those concerned by the problem) - Child of interactionConcernedNode */
  xmlNode* interactionListNode;

  /** node named Solver - Child of rootNode */
  xmlNode* solverNode;
  /** solverXML object */
  SolverXML* solverXML;
  /** bool to check whether solverXML has been allocated inside the class or not*/
  bool isSolverXMLAllocatedIn;

public:

  /** Default constructor
  */
  OneStepNSProblemXML();

  /** Build a OneStepNSProblemXML object from a DOM tree describing a OneStepNSProblem
  *   \param OneStepNSProblemNode : the OneStepNSProblem DOM tree
  *   \exception XMLException : if a property of the OneStepNSProblemXML lacks in the DOM tree
  */
  OneStepNSProblemXML(xmlNode *);

  /** Destructor
  */
  virtual ~OneStepNSProblemXML();

  /** Return the type of the OneStepNSProblem
  *   \return a string
  */
  inline std::string  getNSProblemType() const
  {
    //std::string  type((char*)problemTypeNode->name);
    return (char*)problemTypeNode->name;
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
  void setDimNSProblem(const int&);

  /** returns true if dimNode is defined
  *  \return a bool
  */
  inline bool hasDim() const
  {
    return (dimNode != NULL);
  }

  /** All is an attribute of the DS_Concerned tag
  *  \return bool : true if attribute all is defined
  */
  bool hasAll() const;

  /** to set the attribute "all" of the Interaction_concerned tag
  *   \param bool : the value to assign to the attribute
  */
  void setAll(const bool&);

  /** to xml object that handles solver
  *   \return a pointer to a SolverXML
  */
  inline SolverXML* getSolverXMLPtr() const
  {
    return solverXML;
  }

  /** set xml object that handles solver
  *   \param a pointer to a SolverXML
  */
  void setSolverXMLPtr(SolverXML *);

  //============================================================
  //      Solver tags and attributes of the OneStepNSProblem
  //               ( according to SICONOS/Numerics )
  //============================================================
  /** checks if tag Solver exists
  *  \return bool : true if tag Solver exists
  */
  inline bool hasSolver() const
  {
    return (solverNode != NULL);
  }

  /** set the solver for the OneStepNSProblem
  *   \param string : the type of solver
  *   \param string : the norm type used by the solver
  *   \param double : the tolerance parameter used by the solver
  *   \param unsigned int : the maximum iteration parameter used by the solver
  *   \param double : the search direction parameter used by the solver
  *   \param double : the regularization parameter
  */
  void setSolver(const std::string&, const std::string&, double, unsigned int, double, double);

  /** makes the operations to create a OneStepNSProblemXML to the SimulationXML
  *   \param xmlNode* : the root node of the OneStepNSProblemXML
  *   \param OneStepNSProblem* : the OneStepNSProblem of this OneStepNSProblemXML
  */
  void updateOneStepNSProblemXML(xmlNode* , OneStepNSProblem*);

  /** Return the id of the OneStepNSProblem (attribute of the root node)
  *   \return a string
  */
  inline std::string getId() const
  {
    return SiconosDOMTreeTools::getStringAttributeValue(rootNode, "id");
  }

  /** Return true if the id attribute of the rootNode is present.
  *   \return a bool   */
  inline bool hasId() const
  {
    return SiconosDOMTreeTools::hasAttributeValue(rootNode, "id");
  }

};


#endif

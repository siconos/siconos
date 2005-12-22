/* Siconos version 1.0, Copyright INRIA 2005.
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
/** \class SolverXML
 *   \brief XML data management for class Solver
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 20/12/2005
 *
 */

#ifndef __SOLVERXML__
#define __SOLVERXML__

#include "SiconosDOMTreeTools.h"
#include "Solver.h"

class Solver;

class SolverXML
{
protected:

  /** root node named "Solver" */
  xmlNodePtr rootNode;

  /** root node for SOLVER formalisation type of the problem (LcpSolving, ...) - Child of solverNode */
  xmlNodePtr solvingFormalisationNode;

  /** Name of algorithm used to solved the problem - Child of solvingFormalisationNode */
  xmlNodePtr solverAlgorithmNode;

  /** \fn SolverXML()
   *  \brief default constructor
   */
  SolverXML();

public:

  /** \fn SolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr)
   *  \brief constructor using xmlNode input
   *  \param : rootNode
   */
  SolverXML(xmlNodePtr);

  /** \fn SolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr)
   *  \brief constructor using xmlNode input for all xmlNode
   *  This constructor should be removed thanks to factories
   *  \param : rootNode
   *  \param : solvingFormalisationNode
   *  \param : solverAlgorithmNode
   */
  SolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr);

  /** \fn ~SolverXML()
   *  \brief destructor
   */
  virtual ~SolverXML();

  /** \fn xmlNode* getRootNode()
   *   \brief Return the node corresponding to the Solver
   *   \return an xmlNodePtr
   */
  inline xmlNodePtr getRootNode() const
  {
    return rootNode;
  }

  /** \fn string getSolvingFormalisation() const
   *   \brief Return contain of node solving formalisation
   *   \return a string
   */
  inline std::string getSolvingFormalisation() const
  {
    std::string type((char*) solvingFormalisationNode->name);
    return type;
  }

  /** \fn bool hasSolvingFormalisation()
   *  \brief checks if tag Solver contains a solving formalisation type
   *  \return a bool
   */
  inline bool hasSolvingFormalisation() const
  {
    return (solvingFormalisationNode != NULL);
  }

  /** \fn string getSolverAlgorithmName() const
   *   \brief return name of solve algorithm
   *   \return a string
   */
  inline std::string getSolverAlgorithmName() const
  {
    std::string type((char*) solverAlgorithmNode->name);
    return type;
  }

  /** \fn bool hasSolverAlgorithm()
   *  \brief checks if tag Solver contains a solving method
   *  \return a bool
   */
  inline bool hasSolverAlgorithm() const
  {
    return (solverAlgorithmNode != NULL);
  }

};


#endif

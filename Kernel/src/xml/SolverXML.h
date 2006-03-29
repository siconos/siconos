/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
 *   \version 1.1.4.
 *   \date 20/12/2005
 *
 *  The only member of this class is the adress of node corresponding to tag "Solver",
 *  child of node "OneStepNSProblem.
 *  Member functions provide reading of attributes of this node:
 *   - Solver type (among those delivered by Numerics)
 *   - Solver parameters
 *   The tag solver is optional and, if omitted, solver type and parameters are set to default values (see Solver.h)
 */

#ifndef __SOLVERXML__
#define __SOLVERXML__

#include "SiconosDOMTreeTools.h"
#include "Solver.h"

class Solver;

class SolverXML
{
protected:

  /** root node named "Solver" - Child of OneStepNSProblem node */
  xmlNodePtr rootNode;

public:

  /** \fn SolverXML(xmlNodePtr = NULL)
   *  \brief default constructor
   *  \param : rootNode
   */
  SolverXML(xmlNodePtr = NULL);

  /** \fn ~SolverXML()
   *  \brief destructor
   */
  ~SolverXML();

  /** \fn xmlNode* getRootNode()
   *  \brief Return the node corresponding to the Solver
   *  \return an xmlNodePtr
   */
  inline xmlNodePtr getRootNode() const
  {
    return rootNode;
  }

  // Functions to get attributes of node "Solver"

  /** \fn string getType() const
   *  \brief return the name of the solver algorithm
   *  \return a string
   */
  std::string getType() const;

  /** \fn unsigned int getMaxIter() const
   *  \brief Return the maximum number of iteration the algorithm can do
   *  \return unsigned int
   */
  unsigned int getMaxIter() const;

  /** \fn double getTolerance() const
   *  \brief Return the tolerance value for the algorithm
   *  \return double
   */
  double getTolerance() const;

  /** \fn string getNormType() const
   *  \brief Return norm type of the algorithm
   *  \return a string
   */
  std::string getNormType() const;

  /** \fn double getSearchDirection() const
   *  \brief Return the searchDirection value
   *  \return double
   */
  double getSearchDirection() const;
};


#endif

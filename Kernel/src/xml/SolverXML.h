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
/*! \file
*/

#ifndef __SOLVERXML__
#define __SOLVERXML__

#include "SiconosDOMTreeTools.h"

class Solver;

/** XML management for Solver
 *
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 2.1.1.
 *   \date 20/12/2005
 *
 *  The only member of this class is the adress of node corresponding to tag "Solver",
 *  child of node "OneStepNSProblem.
 *  Member functions provide reading of attributes of this node:
 *   - Solver type (among those delivered by Numerics)
 *   - Solver parameters
 *   The tag solver is optional and, if omitted, solver type and parameters are set to default values (see Solver.h)
 */
class SolverXML
{
protected:

  /** root node named "Solver" - Child of OneStepNSProblem node */
  xmlNodePtr rootNode;

public:

  /** default constructor
  *  \param : rootNode
  */
  SolverXML(xmlNodePtr = NULL);

  /** destructor
  */
  ~SolverXML();

  /** Return the node corresponding to the Solver
  *  \return an xmlNodePtr
  */
  inline xmlNodePtr getRootNode() const
  {
    return rootNode;
  }

  // Functions to get attributes of node "Solver"

  /** return the name of the solver algorithm
  *  \return a string
  */
  std::string getType() const;

  /** Return the maximum number of iteration the algorithm can do
  *  \return unsigned int
  */
  unsigned int getMaxIter() const;

  /** Return the tolerance value for the algorithm
  *  \return double
  */
  double getTolerance() const;

  /** Return the verbose mode variable
  *  \return unsigned int
  */
  unsigned int getVerbose() const;

  /** Return norm type of the algorithm
  *  \return a string
  */
  std::string getNormType() const;

  /** Return the searchDirection value
  *  \return double
  */
  double getSearchDirection() const;

  /** Return the rho value for the algorithm
  *  \return double
  */
  double getRho() const;

};


#endif

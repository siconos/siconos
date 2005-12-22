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
/** \class NSQPSolverXML
 *   \brief XML data management for class NSQPSolver
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 20/12/2005
 *
 */


#ifndef __NSQPSOLVERXML__
#define __NSQPSOLVERXML__

#include "SolverXML.h"
#include "NSQPSolver.h"

class NSQPSolver;

class NSQPSolverXML : public SolverXML
{
private:

  /** \fn NSQPSolverXML(xmlNodePtr)
   *  \brief default constructor
   */
  NSQPSolverXML::NSQPSolverXML();

public:

  /** \fn NSQPSolverXML(xmlNodePtr)
   *  \brief constructor using xmlNode input
   *  \param a pointer to xmlNode
   */
  NSQPSolverXML(xmlNodePtr);

  /** \fn NSQPSolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr)
   *  \brief constructor using xmlNode input for all xmlNode
   *  This constructor should be removed thanks to factories
   *  \param : rootNode
   *  \param : solvingFormalisationNode
   *  \param : solverAlgorithmNode
   */
  NSQPSolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr);

  /** \fn ~NSQPSolverXML()
   *  \brief destructor
   */
  ~NSQPSolverXML();

  /** \fn double getTolerance() const
   *   \brief Return the tolerance value for the algorithm
   *   \return double
   */
  double getTolerance() const;
};


#endif

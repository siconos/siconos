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
/** \class CPGSolverXML
 *   \brief XML data management for class CPGSolver
 *  \author SICONOS Development Team - copyright INRIA
 *   \version 1.0
 *   \date 20/12/2005
 *
 */


#ifndef __CPGSOLVERXML__
#define __CPGSOLVERXML__

#include "SolverXML.h"
#include "CPGSolver.h"

class CPGSolver;

class CPGSolverXML : public SolverXML
{
private:

  /** \fn CPGSolverXML(xmlNodePtr)
   *  \brief default constructor
   */
  CPGSolverXML::CPGSolverXML();

public:

  /** \fn CPGSolverXML(xmlNodePtr)
   *  \brief constructor using xmlNode input
   *  \param a pointer to xmlNode
   */
  CPGSolverXML(xmlNodePtr);

  /** \fn CPGSolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr)
   *  \brief constructor using xmlNode input for all xmlNode
   *  This constructor should be removed thanks to factories
   *  \param : rootNode
   *  \param : solvingFormalisationNode
   *  \param : solverAlgorithmNode
   */
  CPGSolverXML(xmlNodePtr, xmlNodePtr, xmlNodePtr);

  /** \fn ~CPGSolverXML()
   *  \brief destructor
   */
  ~CPGSolverXML();

  /** \fn unsigned int getMaxIter() const
   *   \brief Return the maximum number of iteration the algorithm can do
   *   \return unsigned int
   */
  unsigned int getMaxIter() const;

  /** \fn double getTolerance() const
   *   \brief Return the tolerance value for the algorithm
   *   \return double
   */
  double getTolerance() const;
};


#endif

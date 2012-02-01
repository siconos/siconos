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
/*! \file MLCPProjectOnConstraints.hpp
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef MLCPProjectOnConstraints_H
#define MLCPProjectOnConstraints_H

#include "MLCP.hpp"
#include "OSNSMatrixProjectOnConstraints.hpp"
/** Formalization and Resolution of a Mixed Linear Complementarity Problem (MLCP)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.3.0.
 *  \date (Creation) 01, 2011
 *
 * This class is devoted to the formalization and the resolution of the
 * Mixed Linear Complementarity Problem (MLCP) for the specific problem
 * of the projection onto the constraints in Mechanics
 *
 */
class MLCPProjectOnConstraints : public MLCP
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(MLCPProjectOnConstraints);

public:

  /** compute the number of inequality and equality for a given tuple of URs
   * update the global number of equality(_n) and inequality (_m)
   * set up _numerics_problem parameters (blocksLine and blocksIsComp )
   */
  virtual void computeOptions(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2);

  /** constructor from data
  *  \param Solver* pointer to object that contains solver algorithm and formulation \n
  *  (optional, default = NULL => read .opt file in Numerics)
  *  \param String: id of the problem (default = "unamed")
  */
  MLCPProjectOnConstraints(const int newNewNumericsSolverId = SICONOS_MLCP_ENUM);

  /** destructor
  */
  ~MLCPProjectOnConstraints() {};

  /** print the data to the screen
  */
  void display() const;
  void saveInMemory()
  {
    ;
  }
  virtual void initOSNSMatrix();
  virtual void updateUnitaryBlocks();
  virtual void computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor&);
  virtual void computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor&);
  virtual void computeqBlock(SP::UnitaryRelation, unsigned int);
  virtual void postCompute();
  virtual void postComputeLagrangianR(SP::UnitaryRelation, unsigned int);
  virtual void postComputeNewtonEulerR(SP::UnitaryRelation, unsigned int);



};

TYPEDEF_SPTR(MLCPProjectOnConstraints);
#endif // MLCPProjectOnConstraints_H

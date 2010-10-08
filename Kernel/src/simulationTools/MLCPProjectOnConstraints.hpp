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
/*! \file MLCPProjectOnConstraints.h
\brief Linear Complementarity Problem formulation and solving
*/

#ifndef MLCPProjectOnConstraints_H
#define MLCPProjectOnConstraints_H

#include "MLCP.hpp"


class MLCPProjectOnConstraints : public MLCP
{

public:


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
  virtual void computeDiagonalUnitaryBlock(const UnitaryRelationsGraph::VDescriptor&);
  virtual void computeUnitaryBlock(const UnitaryRelationsGraph::EDescriptor&);
  virtual void computeqBlock(SP::UnitaryRelation, unsigned int);
  virtual void postCompute();

};

TYPEDEF_SPTR(MLCPProjectOnConstraints);
#endif // MLCPProjectOnConstraints_H

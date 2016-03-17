/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

/*! \file SiconosBroadphase.hpp
  \brief Definition of an abstract broadphase algorithm.
*/

#ifndef SiconosBroadphase_h
#define SiconosBroadphase_h

#include <SiconosVisitor.hpp>
#include <MechanicsFwd.hpp>

class SiconosBroadphase : public SiconosVisitor
{
protected:
  virtual void visit(SP::SiconosPlane plane) = 0;
  virtual void visit(SP::SiconosSphere sphere) = 0;
  virtual void visit(SP::SiconosBox box) = 0;
  virtual void visit(SP::Contactor contactor) = 0;

public:
  virtual void buildGraph(SP::Contactor contactor) = 0;
  virtual void updateGraph() = 0;
  virtual void performBroadphase() = 0;
};

#endif /* SiconosBroadphase_h */

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
#include <Model.hpp>

class SiconosBroadphase : public SiconosVisitor
{
protected:
  virtual void visit(SP::SiconosPlane plane) = 0;
  virtual void visit(SP::SiconosSphere sphere) = 0;
  virtual void visit(SP::SiconosBox box) = 0;
  //  virtual void visit(SP::BodyDS body) = 0;
  virtual void visit(const BodyDS &body) = 0;

  SP::Model _model;
  
  void link(SP::Interaction inter,
            SP::DynamicalSystem ds1,
            SP::DynamicalSystem ds2 = SP::DynamicalSystem());

  void unlink(SP::Interaction inter);

public:
  virtual void buildGraph(SP::Model model) = 0;
  virtual void buildGraph(std::vector<SP::BodyDS> bodies) = 0;
  virtual void updateGraph() = 0;
  virtual void performBroadphase() = 0;

  SP::Model model() { return _model; }
};

#endif /* SiconosBroadphase_h */

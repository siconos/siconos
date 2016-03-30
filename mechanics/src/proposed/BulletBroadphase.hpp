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

/*! \file BulletBroadphase.hpp
  \brief Definition of a Bullet-based broadphase algorithm.
*/

#ifndef BulletBroadphase_h
#define BulletBroadphase_h

#include <MechanicsFwd.hpp>

#include <SiconosBroadphase.hpp>
#include <SiconosShape.hpp>
#include <SiconosContactor.hpp>

#include <map>

DEFINE_SPTR(BulletBroadphase_impl);

class BulletBroadphase : public SiconosBroadphase, public std11::enable_shared_from_this<BulletBroadphase>
{
protected:
  SP::BulletBroadphase_impl impl;

  // callback for contact point removal, and a global for context
  static bool bulletContactClear(void* userPersistentData);
  static BulletBroadphase *gBulletBroadphase;

public:
  BulletBroadphase();
  ~BulletBroadphase();

protected:
  void visit(SP::SiconosPlane plane);
  void visit(SP::SiconosSphere sphere);
  void visit(SP::SiconosBox box);
  void visit(const BodyDS &body);

  void update(SP::SiconosPlane plane);
  void update(SP::SiconosSphere sphere);
  void update(SP::SiconosBox box);

  template<typename ST, typename BT>
  void visit_helper(ST& shape, BT& btshape,
                    std::map<ST,BT>& shapemap);
  
public:
  // TODO: default implementations of buildGraph to SiconosBroadphase?
  //       encountered problems with shared_from_this() when doing so.
  void buildGraph(SP::Model model);
  void buildGraph(std::vector<SP::BodyDS> bodies);

  void updateGraph();
  void performBroadphase();
};

#endif /* BulletBroadphase.hpp */

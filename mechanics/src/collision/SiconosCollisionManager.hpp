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

/*! \file SiconosCollisionManager.hpp
\brief A mechanics world is a Siconos InteractionManager that supports
  static contactors and dynamic contactors attached to special
  Dynamical Systems (BodyDS, derived from NewtonEulerDS) found in the
  NonSmoothDynamicalSystem.
*/

#ifndef SiconosCollisionManager_h
#define SiconosCollisionManager_h

#include <InteractionManager.hpp>
#include <SiconosContactor.hpp>

class SiconosCollisionManager : public InteractionManager
{
public:
  SiconosCollisionManager() : InteractionManager() {}
  virtual ~SiconosCollisionManager() {}

  /** An opaque handle can be used to refer to a specific static
   * contactor set previously added to the collision manager. */
  typedef void* StaticContactorSetID;

public:
  virtual StaticContactorSetID insertStaticContactorSet(
    SP::SiconosContactorSet cs, SP::SiconosVector position = SP::SiconosVector()) = 0;

  virtual bool removeStaticContactorSet(StaticContactorSetID id) = 0;
};

#endif /* SiconosCollisionManager.hpp */

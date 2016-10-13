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

/*! \file SiconosMechanicsWorld.hpp
\brief A mechanics world is a Siconos InteractionManager that supports
  static contactors and dynamic contactors attached to special
  Dynamical Systems (BodyDS, derived from NewtonEulerDS) found in the
  NonSmoothDynamicalSystem.
*/

#ifndef SiconosMechanicsWorld_h
#define SiconosMechanicsWorld_h

#include <InteractionManager.hpp>

class SiconosMechanicsWorld : public InteractionManager
{
public:
  SiconosMechanicsWorld() {}
  virtual ~SiconosMechanicsWorld() {}

public:
  virtual void insertStaticContactor(SP::SiconosContactor contactor) = 0;

  virtual void insertNonSmoothLaw(SP::NonSmoothLaw, int group1, int group2) = 0;
  virtual SP::NonSmoothLaw nonSmoothLaw(int group1, int group2) = 0;
};

#endif /* SiconosMechanicsWorld.hpp */

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

#ifndef BulletTimeSteppingProjectOnConstraints_hpp
#define BulletTimeSteppingProjectOnConstraints_hpp

#include "TimeSteppingProjectOnConstraints.hpp"
#include "BulletSpaceFilter.hpp"

class BulletTimeSteppingProjectOnConstraints : public TimeSteppingProjectOnConstraints
{

  SP::BulletSpaceFilter _spaceFilter;

public:
  BulletTimeSteppingProjectOnConstraints(SP::BulletSpaceFilter sf,
                                         SP::TimeDiscretisation t,
                                         SP::OneStepIntegrator osi,
                                         SP::OneStepNSProblem osnspb_velo,
                                         SP::OneStepNSProblem osnspb_pos,
                                         unsigned int level = 1) :
    TimeSteppingProjectOnConstraints(t, osi, osnspb_velo, osnspb_pos),
    _spaceFilter(sf) {};


  void updateWorldFromDS();
};

TYPEDEF_SPTR(BulletTimeSteppingProjectOnConstraints);
#endif



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

/*! \file SiconosContactor.hpp
  \brief Definition of an abstract contactor
*/


#ifndef SiconosContactor_h
#define SiconosContactor_h

#include <vector>
#include <utility>

#include "MechanicsFwd.hpp"

#include <SiconosSerialization.hpp>

#include "SiconosShape.hpp"

/** Class to hold the shape assigned to a body, and to associate each
 *  shape with an offset and collision group. */

struct SiconosContactor
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosContactor);

public:
  SiconosContactor(SP::SiconosShape _shape,
                   SP::SiconosVector _offset = SP::SiconosVector(),
                   int _collision_group = 0);

  SP::SiconosShape shape;
  SP::SiconosVector offset;
  int collision_group;
};

class SiconosContactorSet : public std::vector< SP::SiconosContactor >
{
public:
  typedef std::vector< SP::SiconosContactor >::iterator iterator;

  void append(SP::SiconosContactor b) { push_back(b); }
  void append(std::vector<SP::SiconosContactor> b) { insert(end(), b.begin(), b.end()); }
  void append(const SiconosContactorSet& b) { insert(end(), b.begin(), b.end()); }
  void append(const SP::SiconosContactorSet& b) { insert(end(), b->begin(), b->end()); }
};

#endif /* SiconosContactor_h */

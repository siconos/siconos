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

#include "MechanicsFwd.hpp"

#include <SiconosSerialization.hpp>
#include <SiconosVisitor.hpp>

#include "SiconosShape.hpp"

// NEW APPROACH: No inheritance on SiconosContactor/Sphere, etc.  Just create SiconosContactor
// and Shape descriptors, then in buildInteractions, "compile" this down to a
// Bullet-specific graph.

// Can we share SiconosContactors between BodyDS instances?  It would be best if the
// SiconosContactor did *not* have pointers back to the BodyDS.  And yet, a
// btCollisionShape must be associated with each shape...

// Support groups, NSLs per surface, ..

class SiconosContactor : public std11::enable_shared_from_this<SiconosContactor>
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(SiconosContactor);

public:
  virtual ~SiconosContactor() {}

  struct ShapeOffset {
    ShapeOffset(SP::SiconosShape _shape, SP::SiconosVector _offset)
      : shape(_shape), offset(_offset) {}
    SP::SiconosShape shape;
    SP::SiconosVector offset;
  };

  const std::vector<ShapeOffset> &shapes() const
    { return _shapes; }

  virtual void addShape(SP::SiconosShape shape,
                        SP::SiconosVector position);

  void setPosition(const SP::SiconosVector position);

  /** visitors hook
   */
  ACCEPT_VISITORS();

protected:
  std::vector<ShapeOffset> _shapes;
};

#endif /* SiconosContactor_h */

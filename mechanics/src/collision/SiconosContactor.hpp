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

  void addShape(SP::SiconosPlane shape, SP::SiconosVector offset);
  void addShape(SP::SiconosSphere shape, SP::SiconosVector offset);
  void addShape(SP::SiconosBox shape, SP::SiconosVector offset);
  void addShape(SP::SiconosConvexHull shape, SP::SiconosVector offset);

  const std::vector< std::pair<SP::SiconosPlane,
                               SP::SiconosVector> >& planes() const
    { return _planes; }

  const std::vector< std::pair<SP::SiconosSphere,
                               SP::SiconosVector> >& spheres() const
    { return _spheres; }

  const std::vector< std::pair<SP::SiconosBox,
                               SP::SiconosVector> >& boxes() const
    { return _boxes; }

  const std::vector< std::pair<SP::SiconosConvexHull,
                               SP::SiconosVector> >& convexhulls() const
    { return _chs; }

  /** visitors hook
   */
  ACCEPT_VISITORS();

protected:
  std::vector< std::pair<SP::SiconosPlane,  SP::SiconosVector> > _planes;
  std::vector< std::pair<SP::SiconosSphere, SP::SiconosVector> > _spheres;
  std::vector< std::pair<SP::SiconosBox, SP::SiconosVector> > _boxes;
  std::vector< std::pair<SP::SiconosConvexHull, SP::SiconosVector> > _chs;
};

#endif /* SiconosContactor_h */

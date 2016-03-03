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

/*! \file SiconosShape.hpp
  \brief Definition of an abstract rigid shape
*/


#ifndef SiconosShape_h
#define SiconosShape_h

#include "MechanicsFwd.hpp"
#include <SiconosSerialization.hpp>
#include <SiconosVisitor.hpp>

class SiconosShapeHandler
{
public:
  virtual void onChanged(SP::SiconosSphere) = 0;
  virtual void onChanged(SP::SiconosBox) = 0;
};

class SiconosShape
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosShape);

  SP::SiconosShapeHandler _handler;
  
  float _x, _y, _z;

  SiconosShape(float x, float y, float z)
    : _x(x), _y(y), _z(z) {}
  
  // TODO orientation

public:
  void setHandler(SP::SiconosShapeHandler handler)
    { _handler = handler; }

  /** visitors hook
   */
  // Get around the fact that ACCEPT_VISITORS() requires
  // enabled_shared_from_this, blocking derivative classes from using it.  Which
  // is fine since this base class' visitor should never be called anyway.
  virtual void acceptSP(SP::SiconosVisitor tourist) = 0;
};

struct SiconosSphere : public SiconosShape, public std11::enable_shared_from_this<SiconosSphere>
{
protected:
  float _radius;

  void setDirty()
    { if (_handler) _handler->onChanged(shared_from_this()); }

public:
  SiconosSphere(float x, float y, float z, float radius)
    : SiconosShape(x,y,z), _radius(radius) {}
  float radius() const { return _radius; }
  void setRadius(float r) { _radius = r; setDirty(); }

  /** visitors hook
   */
  ACCEPT_BASE_VISITORS(SiconosShape);
};

struct SiconosBox : public SiconosShape, public std11::enable_shared_from_this<SiconosBox>
{
protected:

  float _width;
  float _height;
  float _depth;

  void setDirty()
    { /*if (_handler) _handler->onChanged(shared_from_this());*/ }

public:
  SiconosBox(float x, float y, float z,
             float width, float height, float depth)
    : SiconosShape(x,y,z),
      _width(width), _height(height), _depth(depth) {}
};

#endif /* SiconosShape_h */

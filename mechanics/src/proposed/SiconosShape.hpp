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
#include <SiconosVector.hpp>

class SiconosShapeHandler
{
public:
  virtual void onChanged(SP::SiconosSphere) = 0;
  virtual void onChanged(SP::SiconosBox) = 0;
  virtual void onChanged(SP::SiconosPlane) = 0;
};

class SiconosShape
{
protected:
  /** serialization hooks
   */
  ACCEPT_SERIALIZATION(SiconosShape);

  /* We use a weak_ptr to avoid cycles, since _handler may point back to the
   * structures which contain references to shapes. */
  std11::weak_ptr<SiconosShapeHandler> _handler;

  virtual void onChanged() = 0;

  SP::SiconosVector _position;

  SiconosShape(float x, float y, float z)
    : _position(new SiconosVector(3))
  {
    (*_position)(0) = x;
    (*_position)(1) = y;
    (*_position)(2) = z;
  }

public:
  SP::SiconosVector position() const { return _position; }

  void setPosition(SP::SiconosVector pos)
  {
    (*_position)(0) = (*pos)(0);
    (*_position)(1) = (*pos)(1);
    (*_position)(2) = (*pos)(2);
    onChanged();
  }

  void setPosition(float x, float y, float z)
  {
    (*_position)(0) = x;
    (*_position)(1) = y;
    (*_position)(2) = z;
    onChanged();
  }

  // TODO orientation

  void setHandler(SP::SiconosShapeHandler handler)
    { _handler = handler; }

  /** visitors hook
   */
  // Get around the fact that ACCEPT_VISITORS() requires
  // enabled_shared_from_this, blocking derivative classes from using it.  Which
  // is fine since this base class' visitor should never be called anyway.
  virtual void acceptSP(SP::SiconosVisitor tourist) = 0;
};

struct SiconosPlane : public SiconosShape, public std11::enable_shared_from_this<SiconosPlane>
{
protected:
  float _distance;

  virtual void onChanged()
    { SP::SiconosShapeHandler h(_handler.lock());
      if (h) h->onChanged(shared_from_this()); }

public:
  SiconosPlane(float x, float y, float z, float distance)
    : SiconosShape(x,y,z), _distance(distance) {}
  float distance() const { return _distance; }
  void setDistance(float d) { _distance = d; onChanged(); }
  SP::SiconosVector normal() { return position(); }
  void setNormal(float x, float y, float z) { setPosition(x,y,z); }

  /** visitors hook
   */
  ACCEPT_BASE_VISITORS(SiconosShape);
};

struct SiconosSphere : public SiconosShape, public std11::enable_shared_from_this<SiconosSphere>
{
protected:
  float _radius;

  virtual void onChanged()
    { SP::SiconosShapeHandler h(_handler.lock());
      if (h) h->onChanged(shared_from_this()); }

public:
  SiconosSphere(float x, float y, float z, float radius)
    : SiconosShape(x,y,z), _radius(radius) {}
  float radius() const { return _radius; }
  void setRadius(float r) { _radius = r; onChanged(); }

  /** visitors hook
   */
  ACCEPT_BASE_VISITORS(SiconosShape);
};

struct SiconosBox : public SiconosShape, public std11::enable_shared_from_this<SiconosBox>
{
protected:

  SP::SiconosVector _dimensions;

  virtual void onChanged()
    { SP::SiconosShapeHandler h(_handler.lock());
      if (h) h->onChanged(shared_from_this()); }

public:
  SiconosBox(float x, float y, float z,
             float width, float height, float depth)
    : SiconosShape(x,y,z), _dimensions(new SiconosVector(3))
  {
    (*_dimensions)(0) = width;
    (*_dimensions)(1) = height;
    (*_dimensions)(2) = depth;
  }

  SP::SiconosVector dimensions() const { return _dimensions; }

  void setDimensions(SP::SiconosVector pos)
  {
    (*_dimensions)(0) = (*pos)(0);
    (*_dimensions)(1) = (*pos)(1);
    (*_dimensions)(2) = (*pos)(2);
    onChanged();
  }
};

#endif /* SiconosShape_h */

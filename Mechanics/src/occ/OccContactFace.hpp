#ifndef OccContactFace_hpp
#define OccContactFace_hpp

#include "OccContactShape.hpp"

struct OccContactFace : public OccContactShape
{
  OccContactFace() : OccContactShape() {};

  OccContactFace(const OccContactShape& reference_shape, unsigned int index);

  virtual const TopoDS_Shape& contact() const;

  virtual void computeUVBounds();

  unsigned int _index;
  SP::TopoDS_Face _face;

};
#endif

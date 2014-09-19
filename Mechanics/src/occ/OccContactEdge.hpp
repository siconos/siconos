#ifndef OccContactEdge_hpp
#define OccContactEdge_hpp

#include "OccContactShape.hpp"

struct OccContactEdge : public OccContactShape
{
  OccContactEdge() : OccContactShape() {};

  OccContactEdge(const OccContactShape& reference_shape, unsigned int index);

  virtual const TopoDS_Shape& contact() const;

  virtual void computeUVBounds();

  unsigned int _index;
  SP::TopoDS_Shape _edge;

};
#endif

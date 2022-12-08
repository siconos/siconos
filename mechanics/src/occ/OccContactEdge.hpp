#ifndef OccContactEdge_hpp
#define OccContactEdge_hpp

#include "OccContactShape.hpp"

struct OccContactEdge : public OccContactShape
{
  OccContactEdge() : OccContactShape() {};

  OccContactEdge(const OccContactShape& shape, unsigned int index);

  virtual const SPC::TopoDS_Edge contact() const;

  virtual void computeUVBounds();

  unsigned int _index;
  SPC::TopoDS_Edge _edge;

  ACCEPT_STD_VISITORS();

};
#endif

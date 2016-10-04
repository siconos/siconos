#ifndef OccContactFace_hpp
#define OccContactFace_hpp

#include "OccContactShape.hpp"

struct OccContactFace : public OccContactShape
{
  OccContactFace() : OccContactShape() {};

  OccContactFace(const OccContactShape& shape, unsigned int index);

  virtual SPC::TopoDS_Face contact() const;

  virtual void computeUVBounds();

  unsigned int _index;
  SPC::TopoDS_Face _face;

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();
};
#endif

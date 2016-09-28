#ifndef OccContactFace_hpp
#define OccContactFace_hpp

#include "OccContactShape.hpp"
#include "SiconosVisitor.hpp"

struct OccContactFace : public OccContactShape
{
  OccContactFace() : OccContactShape() {};

  OccContactFace(const OccContactShape& shape, unsigned int index);

  virtual SPC::TopoDS_Face contact() const;

  virtual void computeUVBounds();

  using OccContactShape::distance;

  /** Distance to a contact face.
      \param sh2 : the contact face.
      \param normalFromFace1 : normal on first contact shape, default on second.
      \return the distance, contact points and normal in ContactShapeDistance
   */
  virtual SP::ContactShapeDistance distance(
    const OccContactFace& sh2, bool normalFromFace1=false) const;

  /** Distance to a contact edge.
      \param sh2 : the contact edge.
      \param normalFromFace1 : normal on first contact shape, default on second.
      \return the distance, contact points and normal in ContactShapeDistance
   */
  virtual SP::ContactShapeDistance distance(
    const OccContactEdge& sh2, bool normalFromFace1=false) const;

  unsigned int _index;
  SPC::TopoDS_Face _face;

  /** visitors hook
   */
  ACCEPT_STD_VISITORS();
};
#endif

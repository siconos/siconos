#include "OccContactFace.hpp"
#include "ContactShapeDistance.hpp"
#include <cadmbtb.hpp>

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

#include <limits>

OccContactFace::OccContactFace(const OccContactShape& reference_shape,
                               unsigned int index) :
  OccContactShape(reference_shape),
  _index(index),
  _face(reference_shape.face(index))
{
};


SPC::TopoDS_Face OccContactFace::contact() const
{
  return this->face(this->_index);
}

void OccContactFace::computeUVBounds()
{
  BRepTools::UVBounds(*this->_face,
                      this->binf1[0],
                      this->bsup1[0],
                      this->binf1[1],
                      this->bsup1[1]);
}

#include "OccBody.hpp"

SP::ContactShapeDistance OccContactFace::distance(
  const OccContactFace& sh2, bool normalFromFace1) const
{

  SP::OccBody body(new OccBody());

  SP::ContactShapeDistance pdist;
  pdist.reset(new ContactShapeDistance());
  ContactShapeDistance& dist = *pdist;

  dist.value = std::numeric_limits<double>::infinity();

  cadmbtb_distanceFaceFace(*this, sh2,
                           dist.x1, dist.y1, dist.z1,
                           dist.x2, dist.y2, dist.z2,
                           dist.nx, dist.ny, dist.nz,
                           normalFromFace1,
                           dist.value);

  return pdist;
}

SP::ContactShapeDistance OccContactFace::distance(
  const OccContactEdge& sh2, bool normalFromFace1) const
{

  SP::ContactShapeDistance pdist(new ContactShapeDistance());
  ContactShapeDistance& dist = *pdist;

  dist.value = std::numeric_limits<double>::infinity();

  cadmbtb_distanceFaceEdge(*this, sh2,
                           dist.x1, dist.y1, dist.z1,
                           dist.x2, dist.y2, dist.z2,
                           dist.nx, dist.ny, dist.nz,
                           normalFromFace1,
                           dist.value);

  return pdist;

}

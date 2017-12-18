#include "OccContactFace.hpp"
#include "OccUtils.hpp"
#include "ContactShapeDistance.hpp"
#include "cadmbtb.hpp"

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>


OccContactFace::OccContactFace(const OccContactShape& reference_shape,
                               unsigned int index) :
  OccContactShape(reference_shape),
  _index(index),
  _face(reference_shape.face(index))
{
  this->computeUVBounds();
};


SPC::TopoDS_Face OccContactFace::contact() const
{
  return this->face(this->_index);
}

void OccContactFace::computeUVBounds()
{

  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_FACE);
  for (unsigned int i=0; i<_index; ++i, exp.Next());
  if (exp.More())
  {
    const TopoDS_Face& face = TopoDS::Face(exp.Current());
    BRepTools::UVBounds(face,
                        this->binf1[0],
                        this->bsup1[0],
                        this->binf1[1],
                        this->bsup1[1]);
  }
}


#include "OccContactFace.hpp"

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>


OccContactFace::OccContactFace(const OccContactShape& reference_shape,
                               unsigned int index) :
  OccContactShape(reference_shape),
  _index(index), _face(new TopoDS_Face())
{
  *this->_face = reference_shape.face(this->_index);
};


const TopoDS_Shape& OccContactFace::contact() const
{
  return *this->_face;
}

void OccContactFace::computeUVBounds()
{
  BRepTools::UVBounds(*this->_face,
                      this->binf1[0],
                      this->bsup1[0],
                      this->binf1[1],
                      this->bsup1[1]);
}


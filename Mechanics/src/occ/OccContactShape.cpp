#include "OccContactShape.hpp"
#include "ContactShapeDistance.hpp"
#include "cadmbtb.hpp"

#include <RuntimeException.hpp>

#include <TopoDS.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Quaternion.hxx>

#include <limits>

//#define DEBUG_MESSAGES 1
#include <debug.h>

OccContactShape::OccContactShape() : _shape(new TopoDS_Shape())
{
}


Geometer::Geometer(const OccContactShape& base) :
  base(base),
  answer(new ContactShapeDistance())
{
  answer->value = std::numeric_limits<double>::infinity();
};


OccContactShape::ContactTypeValue OccContactShape::contactType() const
{
  switch (this->_shape->ShapeType())
  {
  case TopAbs_EDGE: { return OccContactShape::Edge; }
  case TopAbs_FACE: { return OccContactShape::Face; }
  default:
    return OccContactShape::Unknown;
  };


};

void OccContactShape::computeUVBounds()
{
  RuntimeException::selfThrow(
    "OccContactShape::computeUVBounds() : cannot compute UV bounds for this contact shape"
    );
}


SP::ContactShapeDistance OccContactShape::distance(
  const OccContactShape& sh2, bool normalFromFace1) const
{
  Geometer geometer(*this);
  sh2.accept(geometer);

  return geometer.answer;

}

SP::ContactShapeDistance OccContactShape::distance(
  const OccContactFace& sh2, bool normalFromFace1) const
{
  RuntimeException::selfThrow(
    "OccContactShape::distance() : cannot compute distance to this contact face"
    );
  return SP::ContactShapeDistance();
}

SP::ContactShapeDistance OccContactShape::distance(
  const OccContactEdge& sh2, bool normalFromFace1) const
{
  RuntimeException::selfThrow(
    "OccContactShape::distance() : cannot compute distance to this contact edge"
    );

  return SP::ContactShapeDistance();
}


std::string OccContactShape::exportBRepToString() const
{
  std::stringstream out;

  BRepTools::Write(this->data(), out);

  return out.str();
}

void OccContactShape::importBRepFromString(const std::string& brepstr)
{
  std::stringstream in;
  BRep_Builder brep_builder;

  in << brepstr;

  BRepTools::Read(this->data(), in, brep_builder);

}

#include <SiconosVector.hpp>

SPC::TopoDS_Face OccContactShape::face(unsigned int index) const
{
  SP::TopoDS_Face return_value(new TopoDS_Face());

  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_FACE);
  for (unsigned int i=0; i<index; ++i, exp.Next());
  if (exp.More())
  {
    // taking a ref fail!
    *return_value = TopoDS::Face(exp.Current());
  }
  else
  {
    RuntimeException::selfThrow("OccContactShape::face failed");
  }

  return return_value;
}

SPC::TopoDS_Edge OccContactShape::edge(unsigned int index) const
{
  SP::TopoDS_Edge return_value(new TopoDS_Edge());

  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_EDGE);
  for (unsigned int i=0; i<index; ++i, exp.Next());
  if (exp.More())
  {
    // taking a ref fail!
    *return_value = TopoDS::Edge(exp.Current());
  }
  else
  {
    RuntimeException::selfThrow("OccContactShape::edge failed");
  }

  return return_value;
}


void OccContactShape::move(const SiconosVector& q)
{

  const gp_Vec translat = gp_Vec(q(0), q(1), q(2));

  const gp_Quaternion rota = gp_Quaternion(q(4), q(5), q(6), q(3));

  gp_Trsf transfo;

  transfo.SetRotation(rota);
  transfo.SetTranslationPart(translat);

  // this->data().Move(transfo);
  this->data().Location(TopLoc_Location(transfo));

  // cf code from Olivier
  // reset Location to avoid accumulation of TopLoc_Datum3D
//  const TopLoc_Location& aLoc = this->data().Location();
//  const gp_Trsf& T = aLoc.Transformation();
//  TopLoc_Location aLocWithoutList(T);
//  this->data().Location(aLocWithoutList);

}



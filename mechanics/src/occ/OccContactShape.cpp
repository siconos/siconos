#include "OccContactShape.hpp"

#include <RuntimeException.hpp>

#include <TopoDS.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Vec.hxx>
#include <gp_Lin.hxx>
#include <gp_Ax3.hxx>
#include <limits>

//#define DEBUG_MESSAGES 1
#include <debug.h>

#include <boost/math/quaternion.hpp>

OccContactShape::OccContactShape() : _shape(new TopoDS_Shape())
{
}

OccContactShape::ContactTypeValue OccContactShape::contactType() const
{
  switch(this->_shape->ShapeType())
  {
  case TopAbs_EDGE:
  {
    return OccContactShape::Edge;
  }
  case TopAbs_FACE:
  {
    return OccContactShape::Face;
  }
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

  this->computeUVBounds();
}

#include <SiconosVector.hpp>

SPC::TopoDS_Face OccContactShape::face(unsigned int index) const
{
  SP::TopoDS_Face return_value(new TopoDS_Face());

  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_FACE);
  for(unsigned int i=0; i<index; ++i, exp.Next());
  if(exp.More())
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
  for(unsigned int i=0; i<index; ++i, exp.Next());
  if(exp.More())
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



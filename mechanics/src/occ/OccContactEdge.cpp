#include "OccContactEdge.hpp"
#include "ContactShapeDistance.hpp"
#include "cadmbtb.hpp"

#include <RuntimeException.hpp>

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

#include <limits>

OccContactEdge::OccContactEdge(const OccContactShape& shape,
                               unsigned int index) :
  OccContactShape(shape),
  _index(index),
  _edge(shape.edge(index))
{
  this->computeUVBounds();
};

const SPC::TopoDS_Edge OccContactEdge::contact() const
{
  return this->edge(this->_index);
}

void OccContactEdge::computeUVBounds()
{
  TopExp_Explorer exp;
  exp.Init(this->data(),TopAbs_EDGE);
  for (unsigned int i=0; i<_index; ++i, exp.Next());
  if (exp.More())
  {
    const TopoDS_Edge& edge = TopoDS::Edge(exp.Current());
    BRepAdaptor_Curve SC(edge);
    this->binf1[0]=SC.FirstParameter();
    this->bsup1[0]=SC.LastParameter();
    this->binf1[1]=0.;
    this->bsup1[1]=0.;
  }
}



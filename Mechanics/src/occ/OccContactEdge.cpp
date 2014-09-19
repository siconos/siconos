#include "OccContactEdge.hpp"

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>

OccContactEdge::OccContactEdge(const OccContactShape& reference_shape,
                               unsigned int index) :
  OccContactShape(reference_shape),
  _index(index),
  _edge(new TopoDS_Edge())
{
  *_edge = reference_shape.edge(this->_index);
};

const TopoDS_Shape& OccContactEdge::contact() const
{
  return *this->_edge;
}

void OccContactEdge::computeUVBounds()
{
  TopExp_Explorer Ex1;
  Ex1.Init(this->contact(),TopAbs_EDGE);
  const TopoDS_Edge& edge = TopoDS::Edge(Ex1.Current());
  BRepAdaptor_Curve SC(edge);
  this->binf1[0]=SC.FirstParameter();
  this->bsup1[0]=SC.LastParameter();
  this->binf1[1]=0.;
  this->bsup1[1]=0.;
}


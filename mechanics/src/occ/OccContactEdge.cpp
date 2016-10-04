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
};

const SPC::TopoDS_Edge OccContactEdge::contact() const
{
  return this->edge(this->_index);
}

void OccContactEdge::computeUVBounds()
{
  TopExp_Explorer Ex1;
  Ex1.Init(*this->contact(),TopAbs_EDGE);
  const TopoDS_Edge& edge = TopoDS::Edge(Ex1.Current());
  BRepAdaptor_Curve SC(edge);
  this->binf1[0]=SC.FirstParameter();
  this->bsup1[0]=SC.LastParameter();
  this->binf1[1]=0.;
  this->bsup1[1]=0.;
}



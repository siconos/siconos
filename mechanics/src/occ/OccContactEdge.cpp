/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "OccContactEdge.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <SiconosException.hpp>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <limits>

#include "ContactShapeDistance.hpp"
#include "cadmbtb.hpp"

OccContactEdge::OccContactEdge(const OccContactShape& shape, unsigned int index)
    : OccContactShape(shape), _index(index), _edge(shape.edge(index)) {
  computeUVBounds();
};

const SPC::TopoDS_Edge OccContactEdge::contact() const { return edge(_index); }

void OccContactEdge::computeUVBounds() {
  TopExp_Explorer exp;
  exp.Init(data(), TopAbs_EDGE);
  for (unsigned int i = 0; i < _index; ++i, exp.Next());
  if (exp.More()) {
    const TopoDS_Edge& edge = TopoDS::Edge(exp.Current());
    BRepAdaptor_Curve SC(edge);
    binf1[0] = SC.FirstParameter();
    bsup1[0] = SC.LastParameter();
    binf1[1] = 0.;
    bsup1[1] = 0.;
  }
}

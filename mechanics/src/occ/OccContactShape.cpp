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
#include "OccContactShape.hpp"

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Ax3.hxx>
#include <gp_Lin.hxx>
#include <gp_Vec.hxx>
#include <limits>

#include "SiconosException.hpp"

// #define DEBUG_MESSAGES 1
#include <boost/math/quaternion.hpp>

#include "siconos_debug.h"

OccContactShape::OccContactShape() : _shape(new TopoDS_Shape()) {}

OccContactShape::ContactTypeValue OccContactShape::contactType() const {
  switch (this->_shape->ShapeType()) {
    case TopAbs_EDGE: {
      return OccContactShape::Edge;
    }
    case TopAbs_FACE: {
      return OccContactShape::Face;
    }
    default:
      return OccContactShape::Unknown;
  };
};

void OccContactShape::computeUVBounds() {
  THROW_EXCEPTION(
      "OccContactShape::computeUVBounds() : cannot compute UV bounds for this contact shape");
}

std::string OccContactShape::exportBRepToString() const {
  std::stringstream out;

  BRepTools::Write(this->data(), out);

  return out.str();
}

void OccContactShape::importBRepFromString(const std::string& brepstr) {
  std::stringstream in;
  BRep_Builder brep_builder;

  in << brepstr;

  BRepTools::Read(this->data(), in, brep_builder);

  this->computeUVBounds();
}

#include <SiconosVector.hpp>

SPC::TopoDS_Face OccContactShape::face(unsigned int index) const {
  SP::TopoDS_Face return_value(new TopoDS_Face());

  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_FACE);
  for (unsigned int i = 0; i < index; ++i, exp.Next())
    ;
  if (exp.More()) {
    // taking a ref fail!
    *return_value = TopoDS::Face(exp.Current());
  } else {
    THROW_EXCEPTION("OccContactShape::face failed");
  }

  return return_value;
}

SPC::TopoDS_Edge OccContactShape::edge(unsigned int index) const {
  SP::TopoDS_Edge return_value(new TopoDS_Edge());

  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_EDGE);
  for (unsigned int i = 0; i < index; ++i, exp.Next())
    ;
  if (exp.More()) {
    // taking a ref fail!
    *return_value = TopoDS::Edge(exp.Current());
  } else {
    THROW_EXCEPTION("OccContactShape::edge failed");
  }

  return return_value;
}

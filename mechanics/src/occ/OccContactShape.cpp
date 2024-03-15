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
#include <boost/math/quaternion.hpp>

#include "SiconosException.hpp"
#include "SiconosPointers.hpp"  // For createSPtr
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

siconos::mechanics::occ::OccContactShape::OccContactShape(TopoDS_Shape& shape)
    : _shape{siconos::pointers::createSPtr(shape)} {};

siconos::mechanics::occ::OccContactShape::OccContactShape(const TopoDS_Shape& shape)
    : _shape{siconos::pointers::createSPtr(const_cast<TopoDS_Shape&>(shape))} {};

void siconos::mechanics::occ::OccContactShape::setShape(std::shared_ptr<TopoDS_Shape> shape) {
  _shape = shape;
  computeUVBounds();
}

void siconos::mechanics::occ::OccContactShape::setData(TopoDS_Shape& data) {
  _shape = siconos::pointers::createSPtr(data);
  computeUVBounds();
}

siconos::mechanics::occ::OccContactShape::ContactTypeValue
siconos::mechanics::occ::OccContactShape::contactType() const {
  switch (_shape->ShapeType()) {
    case TopAbs_EDGE: {
      return ContactTypeValue::Edge;
    }
    case TopAbs_FACE: {
      return ContactTypeValue::Face;
    }
    default:
      return ContactTypeValue::Unknown;
  };
};

void siconos::mechanics::occ::OccContactShape::computeUVBounds() {
  THROW_EXCEPTION(
      "siconos::mechanics::occ::OccContactShape::computeUVBounds() : cannot compute UV bounds "
      "for this "
      "contact shape");
}

std::string siconos::mechanics::occ::OccContactShape::exportBRepToString() const {
  std::stringstream out;

  BRepTools::Write(data(), out);

  return out.str();
}

void siconos::mechanics::occ::OccContactShape::importBRepFromString(
    const std::string& brepstr) {
  std::stringstream in;
  BRep_Builder brep_builder;

  in << brepstr;

  BRepTools::Read(data(), in, brep_builder);

  computeUVBounds();
}

std::shared_ptr<TopoDS_Face> siconos::mechanics::occ::OccContactShape::face(
    unsigned int index) const {
  auto return_value = std::make_shared<TopoDS_Face>();

  TopExp_Explorer exp{data(), TopAbs_FACE};
  for (unsigned int i = 0; i < index; ++i, exp.Next())
    ;
  if (exp.More()) {
    // taking a ref fail!
    *return_value = TopoDS::Face(exp.Current());
  } else {
    THROW_EXCEPTION("siconos::mechanics::occ::OccContactShape::face failed");
  }

  return return_value;
}

std::shared_ptr<TopoDS_Edge> siconos::mechanics::occ::OccContactShape::edge(
    unsigned int index) const {
  auto return_value = std::make_shared<TopoDS_Edge>();

  TopExp_Explorer exp;
  exp.Init(data(), TopAbs_EDGE);
  for (unsigned int i = 0; i < index; ++i, exp.Next())
    ;
  if (exp.More()) {
    // taking a ref fail!
    *return_value = TopoDS::Edge(exp.Current());
  } else {
    THROW_EXCEPTION("siconos::mechanics::occ::OccContactShape::edge failed");
  }

  return return_value;
}

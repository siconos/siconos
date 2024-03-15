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
#include "OccContactFace.hpp"

#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

#include "ContactShapeDistance.hpp"
#include "OccUtils.hpp"
#include "cadmbtb.hpp"

OccContactFace::OccContactFace(const OccContactShape& reference_shape, unsigned int index)
    : OccContactShape(reference_shape), _index(index), _face(reference_shape.face(index)) {
  // Note FP: what's the point of copying the input shape rather than just "pointer-link" it?
  this->computeUVBounds();
};

SPC::TopoDS_Face OccContactFace::contact() const { return this->face(this->_index); }

void OccContactFace::computeUVBounds() {
  TopExp_Explorer exp;
  exp.Init(this->data(), TopAbs_FACE);
  for (unsigned int i = 0; i < _index; ++i, exp.Next())
    ;
  if (exp.More()) {
    const TopoDS_Face& face = TopoDS::Face(exp.Current());
    BRepTools::UVBounds(face, this->binf1[0], this->bsup1[0], this->binf1[1], this->bsup1[1]);
  }
}

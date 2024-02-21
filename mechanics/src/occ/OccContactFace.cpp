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

#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
// #include "OccUtils.hpp"
// #include "ContactShapeDistance.hpp"
// #include "cadmbtb.hpp"

// #include <BRepAdaptor_Surface.hxx>
// #include <BRep_Tool.hxx>

siconos::mechanics::occ::OccContactFace::OccContactFace(const OccContactShape& shape, unsigned int index)
    : OccContactShape(shape), _index(index), _face(shape.face(index)) {

      std::cout << "FACE O?DEX CONSTRUC \n";


  computeUVBounds();
};

std::shared_ptr<const TopoDS_Face> siconos::mechanics::occ::OccContactFace::contact() const {
  return face(_index);
}

void siconos::mechanics::occ::OccContactFace::computeUVBounds() {
  TopExp_Explorer exp;
  exp.Init(data(), TopAbs_FACE);
  for (unsigned int i = 0; i < _index; ++i, exp.Next())
    ;
  if (exp.More()) {
    const TopoDS_Face& face = TopoDS::Face(exp.Current());
    BRepTools::UVBounds(face, binf1[0], bsup1[0], binf1[1], bsup1[1]);
  }
}

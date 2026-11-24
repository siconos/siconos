/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */
#include "OccContactFace.hpp"

#include <BRepTools.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>

siconos::mechanics::occ::OccContactFace::OccContactFace(const OccContactShape& shape,
                                                        int index)
    : OccContactShape(shape), _index(index), _face(shape.face(index)) {
  // Note FP: what's the point of copying the input shape rather than just "pointer-link" it?

  computeUVBounds();
};

std::shared_ptr<TopoDS_Face> siconos::mechanics::occ::OccContactFace::contact() const {
  return face(_index);
}

void siconos::mechanics::occ::OccContactFace::computeUVBounds() {
  TopExp_Explorer exp{data(), TopAbs_FACE};
  for (auto i = 0; i < _index; ++i, exp.Next());
  if (exp.More()) {
    const TopoDS_Face& face = TopoDS::Face(exp.Current());
    BRepTools::UVBounds(face, binf1[0], bsup1[0], binf1[1], bsup1[1]);
  }
}

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
#include "OccUtils.hpp"

#include <BRepExtrema_DistShapeShape.hxx>
#include <cadmbtb.hpp>
#include <gp_Quaternion.hxx>

#include "ContactShapeDistance.hpp"
#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "SiconosException.hpp"
#include "SiconosVector.hpp"

void siconos::mechanics::occ::occ_move(TopoDS_Shape& shape,
                                       const siconos::algebra::SiconosVector& q) {
  const gp_Vec translat{q(0), q(1), q(2)};
  const gp_Quaternion rota{q(4), q(5), q(6), q(3)};

  gp_Trsf transfo;
  transfo.SetRotation(rota);
  transfo.SetTranslationPart(translat);

  shape.Move(transfo);
  shape.Location(TopLoc_Location(transfo));
}

auto siconos::mechanics::occ::occ_distanceFaceFace(const OccContactFace& csh1,
                                                   const OccContactFace& csh2) -> ContactShapeDistance {
  // need the 2 sp pointers to keep memory
  const TopoDS_Face& face1 = *csh1.contact();
  const TopoDS_Face& face2 = *csh2.contact();

  BRepExtrema_DistShapeShape measure{face1, face2};
  auto isDone = measure.Perform();
  ContactShapeDistance dist{};

  if (isDone) {
    /* we look for the first solution on a face */
    auto nb_solutions =
        measure.NbSolution();  // the number of solutions satisfying the minimum distance
    for (Standard_Integer i = 1; i <= nb_solutions; ++i) {
      if (measure.SupportTypeShape2(i) == BRepExtrema_IsInFace) {
        // if the i-th solution on the second shape is inside a face

        dist.point1 = measure.PointOnShape1(i);
        dist.point2 = measure.PointOnShape2(i);

        Standard_Real u, v;

        measure.ParOnFaceS2(i, u, v);
        dist.normal = cadmbtb::tools::FaceNormal(face2, u, v);
        /**check orientation of normal from face 2**/
        dist.oriantates();
        dist.value = measure.Value();
        break;
      }
    }
  } else
    THROW_EXCEPTION("occ distance: BRepExtrema_DistShapeShape failed");

  return dist;  // RVO, no copy
}

auto siconos::mechanics::occ::occ_distanceFaceEdge(const OccContactFace& csh1,
                                                   const OccContactEdge& csh2) -> ContactShapeDistance {
  // need the 2 sp pointers to keep memory
  const TopoDS_Face& face = *csh1.contact();
  const TopoDS_Edge& edge = *csh2.contact();

  ContactShapeDistance dist{};
  BRepExtrema_DistShapeShape measure{face, edge};
  auto isDone = measure.Perform();

  if (isDone) {
    // the number of solutions satisfying the minimum distance
    auto nb_solutions = measure.NbSolution();
    for (Standard_Integer i = 1; i <= nb_solutions; ++i) {
      /* we look for the first solution on a face */
      if (measure.SupportTypeShape1(i) == BRepExtrema_IsInFace) {
        dist.point1 = measure.PointOnShape1(i);
        dist.point2 = measure.PointOnShape2(i);

        Standard_Real u, v;

        measure.ParOnFaceS1(1, u, v);
        dist.normal = cadmbtb::tools::FaceNormal(face, u, v);
        dist.oriantates();
        dist.value = measure.Value();
        break;
      }
    }
    // what to do now if MinDist is not changed ?
  } else
    THROW_EXCEPTION("occ distance: BRepExtrema_DistShapeShape failed");

  return dist;  // RVO, no copy
}

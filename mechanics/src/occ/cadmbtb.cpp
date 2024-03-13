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
#include "cadmbtb.hpp"

#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
// #include <BRep_Tool.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <gp_Dir.hxx>
#include <gp_Quaternion.hxx>
#include <gp_Vec.hxx>
#include <iostream>

#include "ContactShapeDistance.hpp"
#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "SiconosConfig.h"
#include "SiconosException.hpp"
#include "SiconosFortran.h"  // For n2qn1
// #define DEBUG_MESSAGES 1
#include "siconos_debug.h"

// adapted from _CADMBTB_getMinDistanceFace*_using_n2qn1  (Olivier Bonnefon)
auto siconos::mechanics::occ::cadmbtb::distanceFaceFace(const OccContactFace& csh1,
                                                        const OccContactFace& csh2)
    -> ContactShapeDistance {
  // need the 2 sp pointers to keep memory
  const TopoDS_Face& face1 = *csh1.contact();
  const TopoDS_Face& face2 = *csh2.contact();

  constexpr int N = 4;
  std::array<double, N> bsup{csh1.bsup1[0], csh1.bsup1[1], csh2.bsup1[0], csh2.bsup1[1]};
  std::array<double, N> binf{csh1.binf1[0], csh1.binf1[1], csh2.binf1[0], csh2.binf1[1]};

  std::array<double, N> dxim{};
  std::transform(bsup.begin(), bsup.end(), binf.begin(), dxim.begin(),
                 [](double e1, double e2) { return 1.e-6 * (e1 - e2); });
  std::array<double, N> x{};
  std::transform(bsup.begin(), bsup.end(), binf.begin(), x.begin(),
                 [](double e1, double e2) { return 0.5 * (e1 + e2); });

  auto [f, g] = tools::myf_FaceFace(x, face1, face2);

  auto df1 = f;

  int mode = 1;

  double epsabs = 0;

  // Note FP: work array sizes are not those given in the paper
  // https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m2qn1/m2qn1.pdf
  // Why ?
  std::array<double, N*(N + 9) / 2 + 1> rz;
  std::array<int, 2 * N + 1> iz;

  siconos::fortran::optim::n2qn1(N, x.data(), &f, g->data(), dxim.data(), &df1, &epsabs, &mode,
                                 binf.data(), bsup.data(), iz.data(), rz.data());
  while (mode > 7) {
    auto [f, g] = tools::myf_FaceFace(x, face1, face2);
    siconos::fortran::optim::n2qn1(N, x.data(), &f, g->data(), dxim.data(), &df1, &epsabs,
                                   &mode, binf.data(), bsup.data(), iz.data(), rz.data());
  }

  DEBUG_PRINTF("mode=%d and min value at u=%e,v=%e f=%e\n", mode, x[0], x[1], sqrt(f));
  DEBUG_PRINTF("siconos::mechanics::occ::cadmbtb::distanceFaceFace( dist = %e\n", sqrt(f));
  ContactShapeDistance dist{};
  dist.value = sqrt(f);

  dist.normal = tools::FaceNormal(face2, x[2], x[3]);
  dist.point1 = tools::FacePoint(face1, x[0], x[1]);
  dist.point2 = tools::FacePoint(face2, x[2], x[3]);
  /**check orientation of normal from face 2**/
  dist.oriantates();
  return dist;  // no copy, RVO
}

// template <std::size_t N>
// auto prepare_data_for_n2qn1(const std::array<double, N>& binf,
//                             const std::array<double, N>& bsup) {
//   constexpr std::array<double, N> dxim;
//   std::transform(bsup.begin(), bsup.end(), binf.begin(), dxim.begin(),
//                  [](double e1, double e2) { return 1.e-6 * (e1 - e2); });

//   constexpr std::array<double, N> x;
//   std::transform(bsup.begin(), bsup.end(), binf.begin(), x.begin(),
//                  [](double e1, double e2) { return 0.5 * (e1 + e2); });

//   return std::make_tuple{&dxim, &x};
// }

// template <std::size_t N>
// auto siconos::mechanics::occ::cadmbtb::tools::wrap_n2qn1(const std::array<double, N>& binf,
//                                                          const std::array<double, N>& bsup,
//                                                          const TopoDS_Face& face1,
//                                                          const TopoDS_Edge& edge2) {
//   std::array<double, N> dxim;
//   std::transform(bsup.begin(), bsup.end(), binf.begin(), dxim.begin(),
//                  [](double e1, double e2) { return 1.e-6 * (e1 - e2); });

//   std::array<double, N> x;
//   std::transform(bsup.begin(), bsup.end(), binf.begin(), x.begin(),
//                  [](double e1, double e2) { return 0.5 * (e1 + e2); });

//   double f = 0.;
//   std::array<double, N> g;
//   double epsabs = 0;  // ??

//   myf_FaceEdge(x, f, g, face1, edge2);

//   int mode = 1;
//   double df1 = f;

//   std::array<double, N*(N + 9) / 2 + 1> rz;
//   std::array<int, 2 * N + 1> iz;

//   siconos::fortran::optim::n2qn1(N, x.data(), &f, g.data(), dxim.data(), &df1, &epsabs,
//   &mode,
//                                  binf.data(), bsup.data(), iz.data(), rz.data());
//   while (mode > 7) {
//     // mode == 7 means: failure to factor the Hessian matrix
//     myf_FaceEdge(x, f, g, face1, edge2);

//     siconos::fortran::optim::n2qn1(N, x.data(), &f, g.data(), dxim.data(), &df1, &epsabs,
//                                    &mode, binf.data(), bsup.data(), iz.data(), rz.data());
//   }

//   auto minDist = sqrt(f);

//   return std::make_tuple(x, minDist, g);
// }

auto siconos::mechanics::occ::cadmbtb::distanceFaceEdge(const OccContactFace& csh1,
                                                        const OccContactEdge& csh2)
    -> ContactShapeDistance {
  // need the 2 sp pointers to keep memory
  const TopoDS_Face& face = *csh1.contact();
  const TopoDS_Edge& edge = *csh2.contact();

  constexpr int N = 3;
  std::array<double, N> bsup{csh1.bsup1[0], csh1.bsup1[1], csh2.bsup1[0]};
  std::array<double, N> binf{csh1.binf1[0], csh1.binf1[1], csh2.binf1[0]};

  std::array<double, N> dxim{};
  std::transform(bsup.begin(), bsup.end(), binf.begin(), dxim.begin(),
                 [](double e1, double e2) { return 1.e-6 * (e1 - e2); });
  std::array<double, N> x{};
  std::transform(bsup.begin(), bsup.end(), binf.begin(), x.begin(),
                 [](double e1, double e2) { return 0.5 * (e1 + e2); });

  auto [f, g] = tools::myf_FaceEdge(x, face, edge);

  auto df1 = f;

  int mode = 1;

  double epsabs = 0;

  // Note FP: work array sizes are not those given in the paper
  // https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m2qn1/m2qn1.pdf
  // Why ?
  std::array<double, 4 * (4 + 9) / 2 + 1> rz;
  std::array<int, 2 * 4 + 1> iz;

  siconos::fortran::optim::n2qn1(N, x.data(), &f, g->data(), dxim.data(), &df1, &epsabs, &mode,
                                 binf.data(), bsup.data(), iz.data(), rz.data());
  while (mode > 7) {
    auto [f, g] = tools::myf_FaceEdge(x, face, edge);
    siconos::fortran::optim::n2qn1(N, x.data(), &f, g->data(), dxim.data(), &df1, &epsabs,
                                   &mode, binf.data(), bsup.data(), iz.data(), rz.data());
  }

  ContactShapeDistance dist{};
  dist.value = sqrt(f);

  DEBUG_PRINTF("mode=%d and min value at u=%e,v=%e f=%e\n", mode, x[0], x[1], dist.value);
  DEBUG_PRINTF("cadmbtb_getMinDistanceFaceEdge_using_n2qn1 dist = %e\n", dist.value);

  /** V.A. Normal is always computed from the surface which is safer  */
  dist.normal = tools::FaceNormal(face, x[0], x[1]);
  dist.point1 = tools::FacePoint(face, x[0], x[1]);
  dist.point2 = tools::EdgePoint(edge, x[2]);
  /** check orientation of normal*/
  dist.oriantates();
  return dist;  // RVO, no copy
}

auto siconos::mechanics::occ::cadmbtb::tools::myf_FaceFace(std::span<const double> x,
                                                           const TopoDS_Face& face1,
                                                           const TopoDS_Face& face2)
    -> std::tuple<double, std::unique_ptr<std::array<double, 4>>> {
  assert(x.size() == 4);

  // Face of the BRep topology as a 3D surface
  BRepAdaptor_Surface SF1{face1};
  // Computes the point and the first derivatives on the surface, from parameters x[0] and
  // x[1]
  gp_Pnt point_on_face1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  SF1.D1(x[0], x[1], point_on_face1, aV1u, aV1v);

  // Face of the BRep topology as a 3D surface
  BRepAdaptor_Surface SF2{face2};
  // Computes the point and the first derivatives on the surface, from parameters x[2] and
  // x[3]
  gp_Pnt point_on_face2;
  gp_Vec aV2u;
  gp_Vec aV2v;
  SF2.D1(x[2], x[3], point_on_face2, aV2u, aV2v);

  gp_Vec aVP2P1{point_on_face1.Coord() - point_on_face2.Coord()};
  auto fx = aVP2P1.SquareMagnitude();

  auto gx = std::make_unique<std::array<double, 4>>();

  // Gradient computation
  (*gx)[0] = 2 * aV1u.Dot(aVP2P1);
  (*gx)[1] = 2 * aV1v.Dot(aVP2P1);
  (*gx)[2] = -2 * aV2u.Dot(aVP2P1);
  (*gx)[3] = -2 * aV2v.Dot(aVP2P1);

  return std::make_tuple(fx, std::move(gx));  // No copy, RVO
}

auto siconos::mechanics::occ::cadmbtb::tools::myf_FaceEdge(std::span<const double> x,
                                                           const TopoDS_Face& face,
                                                           const TopoDS_Edge& edge)
    -> std::tuple<double, std::unique_ptr<std::array<double, 3>>> {
  assert(x.size() == 3);

  // Face of the BRep topology as a 3D surface
  BRepAdaptor_Surface SF1{face};
  // Computes the point and the first derivatives on the surface, from parameters x[0] and
  // x[1]
  gp_Pnt point_on_face;
  gp_Vec aV1u;
  gp_Vec aV1v;
  SF1.D1(x[0], x[1], point_on_face, aV1u, aV1v);

  // Creates a Curve to access the geometry of edge
  BRepAdaptor_Curve SC(edge);
  // Computes the point of parameter x[2] on the curve with its first derivative.
  gp_Pnt point_on_curve;
  gp_Vec aV2u;
  SC.D1(x[2], point_on_curve, aV2u);

  gp_Vec aVP2P1{point_on_face.Coord() - point_on_curve.Coord()};
  auto fx = aVP2P1.SquareMagnitude();
  DEBUG_PRINTF("myf %e %e %e %e --> %e\n", x[0], x[1], x[2], x[3], fx);
  auto gx = std::make_unique<std::array<double, 3>>();
  // Gradient computation
  (*gx)[0] = 2 * aV1u.Dot(aVP2P1);
  (*gx)[1] = 2 * aV1v.Dot(aVP2P1);
  (*gx)[2] = -2 * aV2u.Dot(aVP2P1);
  return std::make_tuple(fx, std::move(gx));  // No copy, RVO
}

gp_Pnt siconos::mechanics::occ::cadmbtb::tools::FacePoint(const TopoDS_Face& face,
                                                          Standard_Real u, Standard_Real v) {
  // Face of the BRep topology as a 3D surface
  BRepAdaptor_Surface SF(face);

  // compute the point of parameters u,v on the surface
  gp_Pnt aPaux;
  SF.D0(u, v, aPaux);
  return aPaux;
}

gp_Pnt siconos::mechanics::occ::cadmbtb::tools::EdgePoint(const TopoDS_Edge& edge,
                                                          Standard_Real u) {
  // Creates a Curve to access the geometry of edge
  BRepAdaptor_Curve SC(edge);

  // Computes the point of parameter U.
  gp_Pnt aPaux;
  SC.D0(u, aPaux);

  return aPaux;
}

gp_Dir siconos::mechanics::occ::cadmbtb::tools::FaceNormal(const TopoDS_Face& face,
                                                           Standard_Real u, Standard_Real v) {
  // Get the geometric surface of the face
  Handle(Geom_Surface) surf = BRep_Tool::Surface(face);

  // Initializes the local properties of the surface
  GeomLProp_SLProps props{surf, u, v, 1, 0.01};  // ..., u, v, N, resolution)
  // Note FP (from OCCT doc): N indicates the maximum number of derivations to be done (0, 1,
  // 2 or 3). For example, to compute only the tangent, N should be equal to 1. Resolution is
  // the linear tolerance (it is used to test if a vector is null).

  // Get the normal direction.
  return props.Normal();
}

// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#ifndef ContactShapeDistance_hpp
#define ContactShapeDistance_hpp

#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
namespace siconos::mechanics::occ {

struct ContactShapeDistance {
  double value{0.};

  // double x1{0.};
  // double y1{0.};
  // double z1{0.};

  // double x2{0.};
  // double y2{0.};
  // double z2{0.};

  // double nx{0.};
  // double ny{0.};
  // double nz{0.};

  gp_Pnt point1;
  gp_Pnt point2;
  gp_Dir normal;

  bool orientates() {
    if (gp_Vec{point1.Coord() - point2.Coord()}.Dot(normal) < 0.) {
      normal.Reverse();
      return true;
    }
    return false;
  }
};
}  // namespace siconos::mechanics::occ

#endif

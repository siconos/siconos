// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MBTB_Body.hpp"

#include "StorageTools.hpp"

siconos::mechanisms::MBTB_Body::MBTB_Body(
    Eigen::Ref<siconos::algebra::SiconosVector> q0,
    Eigen::Ref<siconos::algebra::SiconosVector> v0, double mass,
    Eigen::Ref<siconos::algebra::SiconosMatrix> I,
    Eigen::Ref<siconos::algebra::SiconosVector> centerOfMass, const std::string& BodyName,
    const std::string& CADFile)
    : NewtonEulerDS{q0, v0, mass, I, siconos::algebra::alias_t},
      _mBodyName{BodyName},
      _cadFileName{CADFile} {
  centerOfMass_view_ =
      std::make_shared<siconos::algebra::MapVectorType>(centerOfMass.data(), ndof_);
}

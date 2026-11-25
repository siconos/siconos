// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file BodyShapeRecord.hpp
  \brief Definition of an abstract Body shape record
  The objective of this class is to keep associate ds and static body woth the
shape in a contactor.
*/
#include "BodyShapeRecord.hpp"

#include <iostream>

#include "SecondOrderDS.hpp"
#include "SiconosShape.hpp"
#include "SiconosVector.hpp"
#include "StaticBody.hpp"

siconos::collision::BodyShapeRecord::BodyShapeRecord(
    std::shared_ptr<siconos::algebra::SiconosVector> b,
    std::shared_ptr<siconos::modeling::SecondOrderDS> d, std::shared_ptr<SiconosShape> sh,
    std::shared_ptr<SiconosContactor> con, std::shared_ptr<StaticBody> staticCSR)
    : base(b),
      ds(d),
      sshape(sh),
      contactor(con),
      shape_version(sh->version()),
      staticBody(staticCSR) {}

void siconos::collision::BodyShapeRecord::display() const {
  std::cout << "BodyShapeRecord display\n";

  if (ds) {
    std::cout << "ds number: " << ds->number() << "\n";
  }
  if (staticBody) {
    std::cout << "static Body number :" << staticBody->number << "\n";
  }
};

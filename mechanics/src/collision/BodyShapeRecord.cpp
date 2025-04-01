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

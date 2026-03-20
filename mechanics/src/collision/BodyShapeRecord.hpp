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
  The objective of this class is to keep associate ds and static body with the
shape in a contactor.
*/

#ifndef BodyShapeRecord_h
#define BodyShapeRecord_h

#include <memory>

#include "SiconosVector.hpp"

namespace siconos::modeling {

class SecondOrderDS;
}

namespace siconos::collision {

class SiconosShape;
class SiconosContactor;
class StaticBody;

// We need to maintain a record to associate each body with a shape,
// a contactor, and a collision object for each shape type.  We also need
// to access generic shape stuff (group, margin) by a pointer from the
// collision callback, so we need a record base class.
class BodyShapeRecord {
 public:
  BodyShapeRecord(std::shared_ptr<siconos::algebra::SiconosVector> b,
                  std::shared_ptr<siconos::modeling::SecondOrderDS> d,
                  std::shared_ptr<SiconosShape> sh, std::shared_ptr<SiconosContactor> con,
                  std::shared_ptr<StaticBody> staticCSR);

  virtual ~BodyShapeRecord() noexcept = default;

  std::shared_ptr<siconos::algebra::SiconosVector> base{nullptr};
  std::shared_ptr<siconos::modeling::SecondOrderDS> ds{nullptr};
  std::shared_ptr<SiconosShape> sshape{nullptr};
  std::shared_ptr<SiconosContactor> contactor{nullptr};
  unsigned int shape_version{0};
  std::shared_ptr<StaticBody> staticBody{nullptr};

  void display() const;
};
}  // namespace siconos::collision
#endif /* BodyShapeRecord_h */

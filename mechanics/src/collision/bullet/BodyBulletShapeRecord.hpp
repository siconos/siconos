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

/*! \file BodyBulletShapeRecord.hpp
  \brief Body shape record
  The objective of this class is to keep associate ds and static body with the
shape in a contactor.
*/

#ifndef BodyBulletShapeRecord_h
#define BodyBulletShapeRecord_h

#include "BodyShapeRecord.hpp"
#include "BulletDeclarations.h"

namespace siconos::algebra {
class SiconosVector;
}

namespace siconos::collision {
class StaticBody;
class SiconosContactor;

}  // namespace siconos::collision

namespace siconos::collision::bullet {

class BodyBulletShapeRecord : public siconos::collision::BodyShapeRecord {
 private:
  BodyBulletShapeRecord() = delete;
  BodyBulletShapeRecord(const BodyBulletShapeRecord&) = delete;
  BodyBulletShapeRecord(BodyBulletShapeRecord&&) = delete;
  BodyBulletShapeRecord operator=(const BodyBulletShapeRecord&) = delete;
  BodyBulletShapeRecord operator=(BodyBulletShapeRecord&&) = delete;

 public:
  BodyBulletShapeRecord(std::shared_ptr<siconos::algebra::SiconosVector> b,
                        std::shared_ptr<siconos::modeling::SecondOrderDS> d,
                        std::shared_ptr<SiconosShape> sh,
                        std::shared_ptr<btCollisionObject> btobj,
                        std::shared_ptr<SiconosContactor> con,
                        std::shared_ptr<StaticBody> staticCSR)
      : BodyShapeRecord(b, d, sh, con, staticCSR), btobject(btobj) {}
  ~BodyBulletShapeRecord() noexcept = default;

  std::shared_ptr<btCollisionObject> btobject{nullptr};
};

template <class DS, class SHAPE, class BT_SHAPE>
class BodyRecord : public BodyBulletShapeRecord,
                   public std::enable_shared_from_this<BodyRecord<DS, SHAPE, BT_SHAPE>> {
 public:
  std::shared_ptr<SHAPE> shape{nullptr};
  std::shared_ptr<BT_SHAPE> btshape{nullptr};

  BodyRecord(std::shared_ptr<siconos::algebra::SiconosVector> base, std::shared_ptr<DS> ds,
             std::shared_ptr<SHAPE> sh, std::shared_ptr<BT_SHAPE> btsh,
             std::shared_ptr<btCollisionObject> btobj,
             std::shared_ptr<siconos::collision::SiconosContactor> con,
             std::shared_ptr<siconos::collision::StaticBody> staticCSR)
      : BodyBulletShapeRecord(base, ds, sh, btobj, con, staticCSR), shape(sh), btshape(btsh) {}

  ~BodyRecord() noexcept = default;
};

}  // namespace siconos::collision::bullet
#endif

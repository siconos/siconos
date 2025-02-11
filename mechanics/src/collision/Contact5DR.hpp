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

#ifndef Contact5DR_hpp
#define Contact5DR_hpp

#include "NewtonEuler5DR.hpp"

namespace siconos::collision {
class BodyShapeRecord;

class Contact5DR : public siconos::modeling::NewtonEuler5DR {
 private:
  ACCEPT_SERIALIZATION(Contact5DR);

 public:
  /* For users that may require extra information about contacts. */
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordA{nullptr};
  std::shared_ptr<siconos::collision::BodyShapeRecord> bodyShapeRecordB{nullptr};

  /**
       to compute the output y = h(q) of the Relation

       \param[in] q generalized coordinates vector of the dynamical systems (at most 2)
     involved in the relation \param[in,out] y the resulting vector
   */
  void computeh(const siconos::algebra::BlockVector& q,
                Eigen::Ref<siconos::algebra::SiconosVector> y) override;

  /** Update this contact point information.
   *
   *  \param pos1 Position on ds1 in ds1 frame.
   *  \param pos2 Position on ds2 in ds2 frame (or world frame if ds2=null).
   *  \param normal Normal in ds2 frame (or world frame if ds2=null).
   */
  virtual void updateContactPoints(const siconos::algebra::SiconosVector& pos1,
                                   const siconos::algebra::SiconosVector& pos2,
                                   const siconos::algebra::SiconosVector& normal);

  virtual void preDelete() {}
  virtual void accept(modeling::relations::Visitor &tourist) const override { tourist.visit(*this); }

};
}  // namespace siconos::collision
#endif

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

/*! \file SiconosContactor.hpp
  \brief Definition of an abstract contactor and a set of such objects.
*/

#ifndef SiconosContactor_h
#define SiconosContactor_h

#include <atomic>
#include <memory>
#include <vector>

#include "SiconosSerialization.hpp"
#include "SiconosVector.hpp"

namespace siconos::collision {

class SiconosShape;

/** @brief Class to hold the shape assigned to a body
 *
 * Associates each shape with an offset and a collision group
 *
 * Offset and collision group are immutable.
 */
class SiconosContactor {
 private:
  /** global counter for contactors */
  static inline std::atomic<uint32_t> counter_{0};

  /** id of the current contactor */
  const uint32_t id_{counter_++};

  /** The collision shape associated with this contactor. */
  std::shared_ptr<SiconosShape> shape_;

  /** Offset in the body frame */
  const siconos::algebra::SiconosVector offset_;

  /** index to identify the collision group */
  const int collision_group_{0};

  // Rule of five
  SiconosContactor() = delete;
  SiconosContactor(const SiconosContactor&) = delete;
  SiconosContactor(SiconosContactor&&) = delete;
  SiconosContactor& operator=(const SiconosContactor&) = delete;
  SiconosContactor& operator=(SiconosContactor&&) = delete;

 protected:
  ACCEPT_SERIALIZATION(SiconosContactor);

 public:
  /**
   * @brief Constructs a SiconosContactor.
   *
   * @param shape The collision shape. Must not be nullptr.
   * @param offset Position/orientation offset of the contactor relative to the body frame.
               Defaults to identity (no translation, no rotation) = [0, 0, 0, 1, 0, 0, 0].
   * @param collision_group Collision group index. Defaults to 0.
   */
  SiconosContactor(std::shared_ptr<SiconosShape> shape,
                   const siconos::algebra::SiconosVector& input_offset =
                       siconos::algebra::SiconosVector::Unit(7, 3),
                   int collision_gp = 0)
      : shape_{std::move(shape)}, offset_{input_offset}, collision_group_{collision_gp} {}

  /** @returns the unique identifier of this contactor. */
  uint32_t id() const { return id_; }

  /** @returns the collision shape (shared_ptr) */
  std::shared_ptr<SiconosShape> shape() { return shape_; }

  // NoteFP : should we really allow that shape can be changed outside of this class?

  virtual ~SiconosContactor() noexcept = default;
  /** @returns the offset of this contactor in the body frame */
  const siconos::algebra::SiconosVector& offset() const { return offset_; }

  /** @returns the collision group index */
  int collision_group() const { return collision_group_; }
};

using SiconosContactorSet = std::vector<std::shared_ptr<SiconosContactor>>;

}  // namespace siconos::collision
#endif /* SiconosContactor_h */

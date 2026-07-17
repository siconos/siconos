/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2026 INRIA.
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

/** @file namespaces.hpp description of all siconos mechanics namespaces
 * Only useful for doc.
 */

/** @namespace siconos::collision
 *  @brief Collision description and detection tools
 */
namespace siconos::collision {

/** @namespace siconos::collision::native
 * @brief native implementation of collision detection
 */
namespace native {

/** @namespace siconos::collision::native::bodies
 * @brief bodies (that may collide) definitions
 */
namespace bodies {}
}  // namespace native

/** @namespace siconos::collision::bullet
 * @brief bullet (optional) implementation of collision detection
 */
namespace bullet {}

/** @namespace siconos::collision::shape
 * @brief
 */
namespace shape {}

}  // namespace siconos::collision

/** @namespace siconos::fem
 *  @brief Finite-element definitions and tools
 *
 */
namespace siconos::fem {

/** @namespace siconos::fem::cable
 *  @brief Cable modeling and simulation with FEM
 */
namespace cable {}

}  // namespace siconos::fem

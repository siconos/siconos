/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
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

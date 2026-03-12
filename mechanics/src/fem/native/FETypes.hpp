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
#ifndef FETypesH
#define FETypesH
#include <cstddef>  // For size_t

namespace siconos::mechanics::fem {

/** Ids for element types in Siconos */
enum class FiniteElementType : int  // we (try to) follow the gmsh numbering convention.
{
  Unknown = 0,
  L2 = 1,  /** 2-node line */
  T3 = 2,  /** 3-node triangle */
  Q4 = 3,  /** 4-node quadrangle */
  TH4 = 4, /** 4-node tetrahedron */
  H8 = 5,  /** 8-node hexahedron */
  P6 = 6,  /** 6-node prism */
  PY5 = 7, /** 5-node pyramid */
  L3 = 8,  /** 3-node second order line (2 nodes associated with the vertices and 1 with the
              edge) */
  T6 = 9, /** 6-node second order triangle (3 nodes associated with the vertices and 3 with the
             edges) */
  B2 = 10, /** 2D beam (2 nodes) --> warning, not gmsh convention*/
  B3 = 11, /** 3D beam  (2 nodes) --> warning, not gmsh convention*/

};

/** @return the number of dof per node of a given element type
 *
 *  If 0 is returned, it means the element is unknown or not properly implemented in Siconos
 *
 * @param type element type
 */
inline constexpr std::size_t number_of_dof_per_node(FiniteElementType type) {
  switch (type) {
    case FiniteElementType::T3:
    case FiniteElementType::Q4:
      return 2;
    case FiniteElementType::TH4:
      return 3;
    case FiniteElementType::B2:
      return 3;
    case FiniteElementType::B3:
      return 6;
    default:  // i.e. elem not properly implemented in Siconos
      return 0;
  }
}

/** @return the order of an element, from its type
 *
 *  If 0 is returned, it means the element is unknown or not properly implemented in Siconos
 *
 * @param type element type
 */
inline constexpr int element_order(FiniteElementType type) {
  switch (type) {
    case FiniteElementType::T3:
    case FiniteElementType::TH4:
      return 1;
    default:  // i.e. elem not properly implemented in Siconos
      return 0;
  }
}

/** Ids for element families */
enum class FiniteElementFamily { isoparametric };

}  // namespace siconos::mechanics::fem

#endif

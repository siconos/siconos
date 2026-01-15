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
#include <unordered_set>

namespace siconos::mechanics::fem {

/** Ids for element types, following gmsh numbering convention */
enum class FiniteElementType : int  // we follow the gmsh numbering convention.
{
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
};

/** Ids for element families */
enum class FiniteElementFamily { isoparametric };

// Helper function to check if an int is a valid Color
inline bool is_valid_element(int value) {
  static const std::unordered_set<int> valid_values = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  return valid_values.count(value) > 0;
}

}  // namespace siconos::mechanics::fem

#endif

/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright INRIA.
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

#include "siconos/algebra/eigen.hpp"

namespace siconos::algebra {
static_assert(vector<int, 3>::SizeAtCompileTime != Eigen::Dynamic);
static_assert(unbounded_vector<int>::SizeAtCompileTime == Eigen::Dynamic);
static_assert(vector<int, 3>::ColsAtCompileTime == 1);
static_assert(match::fixed_size_matrix<matrix<int, 1, 1>>);
static_assert(!match::fixed_size_matrix<Eigen::Matrix<int, Eigen::Dynamic, 1>>);
static_assert(!match::fixed_size_matrix<unbounded_matrix<int>>);
static_assert(match::fixed_size_vector<vector<int, 3>>);
static_assert(!match::fixed_size_vector<unbounded_vector<int>>);
static_assert(match::variable_size_matrix<unbounded_matrix<int>>);
static_assert(match::variable_size_vector<unbounded_vector<int>>);
static_assert(match::any_matrix<matrix<int, 2, 2>>);
static_assert(match::any_matrix<unbounded_matrix<int>>);
static_assert(match::any_matrix<Eigen::DiagonalMatrix<int, 3>>);
static_assert(!match::matrix<int>);
static_assert(!match::vector<int>);

static_assert(std::is_same_v<matrix<int, 1, 2>, trans_t<matrix<int, 2, 1>>>);
static_assert(std::is_same_v<prod_t<matrix<int, 2, 2>, matrix<int, 2, 1>>, matrix<int, 2, 1>>);
static_assert(std::is_same_v<prod_t<matrix<int, 1, 2>, matrix<int, 2, 1>>, matrix<int, 1, 1>>);
static_assert(
    std::is_same_v<prod_t<matrix<int, 1, 2>, trans_t<matrix<int, 1, 2>>>, matrix<int, 1, 1>>);

}  // namespace siconos::algebra
int main() {};

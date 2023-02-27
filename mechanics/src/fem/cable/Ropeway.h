/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2023 INRIA.
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

/*! \file Ropeway.h

  A sequence of ropes used to described the whole cable between two pulleys
*/

#pragma once

#include <nlohmann/json.hpp>
#include <vector>

namespace siconos::fem::cable {

class Cable;
class Pile;
class Rope;
class Support;
class Point;

class Ropeway {
 private:
  std::vector<Rope> m_ropes = {};
  bool m_down{false};

  Ropeway(const Ropeway &) = delete;
  Ropeway &operator=(const Ropeway &) = delete;
  Ropeway(Ropeway &&) = delete;
  Ropeway &operator=(Ropeway &&) = delete;

  /** Create new support from the lower pile of a given rope

     \param[in] a_rope the considered Rope
     \param[in,out] a_supports the vector of all supports
     \param[in,out] a_pulleyIdx current number of added support (internal counter)
   */
  void addSupport(const Rope &a_rope, std::vector<std::shared_ptr<Support>> &a_supports,
                  int &a_pulleyIdx) const;

 public:
  using ojson = nlohmann::ordered_json;
  Ropeway() = default;
  virtual ~Ropeway() noexcept = default;

  /** Initialize/compute profile of rope segments and save them

     \param a_meca object which handles the geometric and material description of the cable
     \param a_piles the vector of pylons, must be ordered!
     \param nb_nodes number of nodes by segment
     \param a_tol tolerance used to compute initial profile
     \apram a_nmax max number of iterations used to compute initial profile
  */
  void compute(const Cable &a_meca, const std::vector<Pile> &a_piles, int nb_nodes,
               double a_tol = 1e-20, int a_nmax = 20);

  /**
     Create and prepare all supports, from the list of declared ropes
     (a support for each 'lower' pile of each rope).

     \param[in,out] a_supports the vector of all supports
     \param[in,out] a_pulleyIdx current number of added support (internal counter)
   */
  void prepareSupport(std::vector<std::shared_ptr<Support>> &a_supports,
                      int &a_pulleyIdx) const;

  int computeNbNodes(int nb_elem, double L);
  int computeMesh(std::vector<Point> &a_q, std::vector<Point> &a_R, std::vector<double> &a_TS,
                  int q_offset);

  const Pile &get_FirstPile();
  const Pile &get_LastPile();
  double get_T0();
  double get_LastT();
  double get_L();
  const Cable &get_meca0() const;

  int to_json(ojson &j);
  void set_Down(bool a_value);
};
}  // namespace siconos::fem::cable

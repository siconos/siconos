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
#include "Ropeway.h"

#include <algorithm>
#include <iostream>

#include "Cable.h"
#include "Pile.h"
#include "Pulley.h"
#include "Rope.h"

void siconos::fem::cable::Ropeway::compute(const Cable &a_meca,
                                           const std::vector<Pile> &a_piles, int nb_nodes,
                                           double a_tol, int a_nmax) {
  /*
    Returns the positions in the shape of a 3*nb_nodes vector [x,y,z,...,x,y,z],
  the tension in the shape of a nb_nodes vector [T,...T],
  the internal force vector in the shape of a 3*nb_nodes vector [H,V,B,...,H,V,B]
  the support reaction of piles in the shape of a 6*nb_piles vector [px,py,pz,H,V,B]
  Moreover, it provides with admissibilities for each subspan with a np.array which i-th row
  contains cable unknowns for each subspan np.array([L, etaY, etaZ])
  */
  size_t n = a_piles.size();
  m_ropes.clear();
  double T0 = a_meca.get_T0();
  Point R0;
  for (size_t k = 0; k < n; k++) {
    size_t k1 = k;
    if (k < n - 1) {
      k1++;
    }
    Rope r(a_piles[k], a_piles[k1], a_tol, a_nmax);
    r.compute(a_meca, nb_nodes, T0, R0);
    T0 = r.get_LastT();
    R0 = r.get_LastR();
    m_ropes.push_back(r);  //  Copy
  }
}

void siconos::fem::cable::Ropeway::prepareSupport(
    std::vector<std::shared_ptr<Support>> &a_supports, int &a_pulleyIdx) const {
  if (!m_down) {
    for (auto &r : m_ropes) {
      addSupport(r, a_supports, a_pulleyIdx);
    }
  } else {
    for (auto r = m_ropes.rbegin(); r != m_ropes.rend(); r++) {
      addSupport(*r, a_supports, a_pulleyIdx);
    }
  }
}

void siconos::fem::cable::Ropeway::addSupport(
    const Rope &a_rope, std::vector<std::shared_ptr<Support>> &a_supports,
    int &a_pulleyIdx) const {
  const Pile &rPile0 = a_rope.get_pile0();
  if (!rPile0.isStation()) {
    a_supports.push_back(std::make_shared<Support>(rPile0));
    a_supports[a_supports.size() - 1]->prepare(a_rope);
  } else if (a_pulleyIdx < 0) {  // non déjà ajouté
    a_pulleyIdx = a_supports.size();
    a_supports.push_back(std::make_shared<Pulley>(rPile0));
  }
}

int siconos::fem::cable::Ropeway::computeNbNodes(int nb_elem, double L) {
  int N = 0;
  for (auto &r : m_ropes) {
    N += r.computeNbNodes(nb_elem, L);
  }
  return N;
}

int siconos::fem::cable::Ropeway::computeMesh(std::vector<Point> &a_q, std::vector<Point> &a_R,
                                              std::vector<double> &a_TS, int q_offset) {
  int offset = q_offset;
  if (!m_down) {
    for (auto &r : m_ropes) {
      offset += r.computeMesh(a_q, a_R, a_TS, offset);
    }
  } else {
    for (auto r = m_ropes.rbegin(); r != m_ropes.rend(); r++) {
      offset += r->computeMesh(a_q, a_R, a_TS, offset, true);
    }
  }
  return offset;
}

const siconos::fem::cable::Pile &siconos::fem::cable::Ropeway::get_FirstPile() {
  return m_ropes.front().get_pile0();
}

const siconos::fem::cable::Pile &siconos::fem::cable::Ropeway::get_LastPile() {
  return m_ropes.back().get_pile0();
}

double siconos::fem::cable::Ropeway::get_T0() {
  if (m_ropes.size()) {
    return m_ropes.front().get_T0();
  } else
    return 0.0;
}

double siconos::fem::cable::Ropeway::get_LastT() {
  if (m_ropes.size()) {
    return m_ropes.back().get_LastT();
  } else
    return 0.0;
}

double siconos::fem::cable::Ropeway::get_L() {
  double l = 0.0;
  for (auto &r : m_ropes) {
    l += r.get_L();
  }
  return l;
}

const siconos::fem::cable::Cable &siconos::fem::cable::Ropeway::get_meca0() const {
  return m_ropes.front().get_meca();
}

int siconos::fem::cable::Ropeway::to_json(ojson &j) {
  j = {{"ropeway_inc", ojson::array()}, {"q", ojson::array()},
       {"TS", ojson::array()},          {"R", ojson::array()},
       {"SR", ojson::array()},          {"meca_global", ojson::array()}};
  for (auto &r : m_ropes) {
    r.to_json(j);
  }
  return 0;
}

void siconos::fem::cable::Ropeway::set_Down(bool a_value) { m_down = a_value; }

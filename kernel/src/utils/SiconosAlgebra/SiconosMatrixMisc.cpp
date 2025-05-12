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

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "io.hpp"

void siconos::algebra::normInfByColumn(const SiconosMatrix &m, SiconosVector &v) {
  if (v.size() != m.cols()) THROW_EXCEPTION("the given vector does not have the right length");
  for (auto i = 0; i < m.cols(); i++) {
    v(i) = m.col(i).norm();
  }
}

bool siconos::algebra::checkSymmetry(SiconosMatrix &m, double tol) {
  auto m_trans = m.transpose();
  auto tmp = m - m_trans;
  double err = siconos::algebra::normInf(tmp);
  if (siconos::algebra::normInf(m_trans) > 0.0) {
    err /= siconos::algebra::normInf(m_trans);
  }
  // std::cout << "err_rel  ="<< err <<"\n";
  return (err < tol);
}

siconos::algebra::SiconosMatrix siconos::algebra::readMatrixFromFile(
    const std::string &filename, bool ascii) {
  SiconosMatrix m;
  if (ascii) {
    io::read(filename, m, io::ASCII_IN);
  } else {
    io::read(filename, m, io::BINARY_IN);
  }
  return m;
}

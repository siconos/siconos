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

#include "Rope.h"

#include <cmath>  // for pow

#include "MechanicalProperties.h"
#include "Pylon.h"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::fem::cable::Rope::Rope(const Pylon &start_pylon, const Pylon &end_pylon,
                                MechanicalProperties meca_prop, double T0,
                                const siconos::algebra::SiconosVector3 &R0, int nb_nodes,
                                double tol, int max_iter)

    : mechanicalProp_{meca_prop},
      start_pylon_(start_pylon),
      end_pylon_(end_pylon),
      tol_(tol),
      max_iter_(max_iter) {
  m_last_ = (&start_pylon_ == &end_pylon_);

  // Ref: case 2 in C. Bertrand PHD, annexe A.
  mechanicalProp_.set_T(T0);
  if (m_last_) {
    // Last segment specific case (if start and end pylons are the same)
    // Note FP : fake part between last and ... last pylon?
    support_reaction_ = R0;
    nodes_coords_ = start_pylon_.coords();
    TS_.resize(1);
    TS_(0) = T0;
  } else {
    // compute length and slopes of the cable and ensures it satisfies the cable
    // (catenary) equation
    catenaryUnknowns_ = siconos::fem::cable::computeAdmissibilityConditions(
        mechanicalProp_, start_pylon_.coords(), end_pylon_.coords(), max_iter_, tol_);

    nodes_coords_.resize(3 * (nb_nodes - 1));
    // we do not save the coordinates of the last node
    // since it corresponds to the starting point of next rope.
    nodes_coords_.setZero();
    TS_.resize(nb_nodes);
    TS_.setZero();
    internal_forces_.resize(3 * nb_nodes);
    internal_forces_.setZero();
    // computes the initial profile of the
    // --> updates nodes_coords, internal_forces and TS
    computeCatenary(mechanicalProp_, end_pylon_.coords(), catenaryUnknowns_, nb_nodes,
                    nodes_coords_, internal_forces_, TS_);

    // Determines nature of the roller battery (compression or support) according
    // to the sign of the reaction (Ch. Bertrand phd, p169)
    auto drk = nodes_coords_.segment<3>(3) - nodes_coords_.head<3>();
    double dot = drk.dot(internal_forces_.head<3>());

    if (dot < 0) {
      internal_forces_ *= -1.;
    } else if (dot == 0) {
      internal_forces_.setZero();
    }

    // Compute support reaction (at first node of the rope)
    // SR = (R at last node of the previous rope) - (R at first node of this rope)
    support_reaction_ = R0 - internal_forces_.head<3>();
  }
}

int siconos::fem::cable::Rope::computeNumberOfElements(double element_length, double L) {
  if (!m_last_) {
    number_of_elements_ = (int)rint(catenaryUnknowns_[0] / element_length);
  }
  return number_of_elements_;
}

int siconos::fem::cable::Rope::initializeFEM(siconos::algebra::SiconosVector &a_q,
                                             siconos::algebra::SiconosVector &a_R,
                                             siconos::algebra::SiconosVector &a_TS,
                                             int q_offset, bool a_reverse) const {
  if (!m_last_) {
    computeCatenary(mechanicalProp_, end_pylon_.coords(), catenaryUnknowns_,
                    number_of_elements_ + 1, a_q, a_R, a_TS, q_offset, a_reverse);
  }
  return number_of_elements_;
}

double siconos::fem::cable::Rope::initialTension() {
  if (TS_.size() > 0)
    return TS_(0);
  else
    return 0.0;
}

double siconos::fem::cable::Rope::getTensionAtLastNode() {
  if (TS_.size() > 0)
    return TS_.tail(1)(0);
  else
    return 0.0;
}

Eigen::Ref<const siconos::algebra::SiconosVector3> siconos::fem::cable::Rope::getLastR()
    const {
  assert(internal_forces_.size() >= 3);
  return internal_forces_.tail<3>();
}

void siconos::fem::cable::Rope::display() const {
  std::cout << "--- Rope: \n";
  if (m_last_) std::cout << "(last one in the line)\n";
  mechanicalProp_.display();
  start_pylon_.display();
  end_pylon_.display();
  std::cout << "- Tolerance and max iter used to compute admissibility conditions: " << tol_
            << ", " << max_iter_;
  std::cout << "- Number of elements used to discretize the rope: " << number_of_elements_;
  std::cout << "\n- Suppport reaction at contact: ";
  siconos::algebra::print(support_reaction_.transpose());
}

siconos::algebra::SiconosVector3 siconos::fem::cable::guess(
    const siconos::algebra::SiconosVector3 &start,
    const siconos::algebra::SiconosVector3 &end) {
  siconos::algebra::SiconosVector3 delta = end - start;
  auto L = delta.norm();
  delta.normalize();
  delta(0) = L;
  return delta;
}

siconos::algebra::SiconosVector3 siconos::fem::cable::computeAdmissibilityConditions(
    const MechanicalProperties &a_meca, const siconos::algebra::SiconosVector3 &start,
    const siconos::algebra::SiconosVector3 &end, int max_iter, double tol) {
  // Initial guess for length and slopes
  auto cable_inc = siconos::fem::cable::guess(start, end);

  auto diff_coords = end - start;
  double alpha = a_meca.get_alpha();
  double beta = a_meca.get_beta();
  siconos::algebra::SiconosVector3 residu;
  siconos::algebra::SiconosMatrix33 jac;

  bool test = true;
  int current_iteration = 0;
  while (test) {
    current_iteration++;

    //  Compute the residual equation and the jacobian of a fixed-fixed cable with
    // imposed initial tension

    // L, etaY, etaZ = cable_inc
    double L = cable_inc(0);
    double etaY = cable_inc(1);
    double etaZ = cable_inc(2);

    double etaY2 = etaY * etaY;
    double etaZ2 = etaZ * etaZ;
    double sq0 = sqrt(1 + etaZ2 + etaY2);
    double sqL = sqrt(1 + etaY2 + (etaZ + sq0 * alpha * L) * (etaZ + sq0 * alpha * L));
    double logL = log((etaZ + sq0 * alpha * L + sqL) / (etaZ + sq0));

    double inv_alpha_sq0 = 1. / (alpha * sq0);
    double alpha_sq0_L = alpha * sq0 * L;
    double sq03 = sq0 * sq0 * sq0;

    residu(0) = (-((beta * L) / sq0) + diff_coords(0) - logL * inv_alpha_sq0);
    residu(1) = (-((beta * etaY * L) / sq0) + diff_coords(1) - (etaY * logL) * inv_alpha_sq0);
    residu(2) = (-((beta * etaZ * L) / sq0) - (alpha * beta * L * L) / 2.0 +
                 (sq0 - sqL) * inv_alpha_sq0 + diff_coords(2));

    // Elastic catenary jacobian
    // See Ch. Bertrand Phd, Appendix A
    jac(0, 0) = (-(beta / sq0) - 1 / sqL);
    jac(0, 1) = ((etaY * (beta * L -
                          (etaZ * sq0 *
                           (sq0 + alpha * etaZ * L - alpha_sq0_L + alpha * alpha_sq0_L * L -
                            sqL + alpha * L * sqL)) /
                              (alpha * (etaZ + sq0) * sqL * (etaZ + alpha_sq0_L + sqL)) +
                          logL / alpha)) /
                 sq03);
    jac(0, 2) = ((-1 - etaZ2 - etaY2 - alpha * etaZ * sq0 * L + sq0 * sqL +
                  alpha * beta * etaZ * L * sqL + etaZ * sqL * logL) /
                 (alpha * sq03 * sqL));

    jac(1, 0) = (-((beta * etaY) / sq0) - etaY / sqL);
    jac(1, 1) = ((beta * etaY2 * L - beta * (1 + etaZ2 + etaY2) * L -
                  (etaZ * etaY2 * sq0 *
                   (sq0 + alpha * etaZ * L - alpha_sq0_L + alpha * alpha * sq0 * L * L - sqL +
                    alpha * L * sqL)) /
                      (alpha * (etaZ + sq0) * sqL * (etaZ + alpha_sq0_L + sqL)) +
                  (etaY2 * logL) / alpha - ((1 + etaZ2 + etaY2) * logL) / alpha) /
                 sq03);
    jac(1, 2) = (-((etaY * (1 + etaZ2 + etaY2 + alpha * etaZ * sq0 * L - sq0 * sqL -
                            alpha * beta * etaZ * L * sqL - etaZ * sqL * logL)) /
                   (alpha * sq03 * sqL)));

    jac(2, 0) = -((beta * etaZ) / sq0) - alpha * beta * L - (etaZ + alpha_sq0_L) / sqL;
    jac(2, 1) = ((etaZ * etaY * L * (sq0 + beta * sqL)) / (sq03 * sqL));
    jac(2, 2) = (beta * ((etaZ2 * L) / sq03 - L / sq0) -
                 ((1 + etaY2) * L) / ((1 + etaZ2 + etaY2) * sqL));

    // Solve admissibility conditions for the catenary
    auto delta = jac.partialPivLu().solve(residu);

    // Update stop criteria
    test = (current_iteration < max_iter) && (residu.norm() > tol) && (delta.norm() > tol);

    // Update cable_inc
    cable_inc -= delta;
  }
  return cable_inc;  // RVO
}

void siconos::fem::cable::computeCatenary(const MechanicalProperties &a_meca,
                                          const siconos::algebra::SiconosVector3 &end,
                                          const siconos::algebra::SiconosVector3 &cable_inc,
                                          int nb_nodes,
                                          siconos::algebra::SiconosVector &positions,
                                          siconos::algebra::SiconosVector &internal_forces,
                                          siconos::algebra::SiconosVector &tension,
                                          int q_offset, bool a_reverse) {
  //  a_meca:
  //           T0  : initial tension
  //           EA  : rigidity of the cable
  //           rho : density of the cable

  // cable_inc [L, etaY, etaZ]
  //         L    : Unstretched length of the cable
  //         etaY : Initial y-slope
  //         etaZ : Initial z-slope

  double L = cable_inc(0);
  double etaY = cable_inc(1);
  double etaZ = cable_inc(2);
  double etaY2 = etaY * etaY;
  double etaZ2 = etaZ * etaZ;
  double sq0 = sqrt(1 + etaY2 + etaZ2);
  double alpha = a_meca.get_alpha();
  double beta = a_meca.get_beta();
  double sq0_alpha = sq0 * alpha;
  double sqL = sqrt(1 + etaY2 + (etaZ + sq0_alpha * L) * (etaZ + sq0_alpha * L));
  double coef = a_meca.initialTension() / sq0;
  double coeff2 = 1. / (alpha * sq0);

  for (int i = 0; i < nb_nodes; i++) {
    int j = (a_reverse) ? nb_nodes - 1 - i : i;
    double S = L / (nb_nodes - 1) * j;
    double sqS = sqrt(1 + etaY2 + (etaZ + sq0_alpha * S) * (etaZ + sq0_alpha * S));
    double logS = log((etaZ + sq0_alpha * L + sqL) / (etaZ + sq0_alpha * S + sqS));
    double coeff = beta / sq0 * (L - S);
    int pos = 3 * (i + q_offset);
    // Elastic catenary positions
    if (i < nb_nodes - 1) {
      // we do not save the coordinates of the last node
      // since it corresponds to the starting point of next rope.
      positions(pos) = end(0) - coeff - coeff2 * logS;
      positions(pos + 1) = end(1) - etaY * coeff - etaY * coeff2 * logS;
      positions(pos + 2) =
          end(2) - etaZ * coeff - 0.5 * alpha * beta * (L * L - S * S) - coeff2 * (sqL - sqS);
    }
    // Tension
    tension(i + q_offset) = coef * sqS;

    // internal forces
    internal_forces(pos) = coef;
    internal_forces(pos + 1) = coef * etaY;
    internal_forces(pos + 2) = coef * (etaZ + alpha * sq0 * S);
  }
}

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
#include <iostream>

#include "MechanicalProperties.h"
#include "Pylon.h"

siconos::fem::cable::Rope::Rope(const Pylon &a_pile0, const Pylon &a_pile1, double a_tol,
                                int a_nmax)
    : left_pylon_(a_pile0), right_pylon_(a_pile1), tol_(a_tol), max_iter_(a_nmax) {
  m_last_ = (&left_pylon_ == &right_pylon_);
}

void siconos::fem::cable::Rope::compute(const class MechanicalProperties &a_meca, int nb_nodes,
                                        double a_T0, const Point &a_R0) {
  // Ref: case 2 in C. Bertrand PHD, annexe A.

  meca = a_meca;
  meca.set_T(a_T0);
  if (m_last_) {
    // Last segment specific case (if left and right pylons are the same)
    // Note FP : fake part between last and ... last pylon?
    support_reaction_ = a_R0;
    nodes_coords_ = {left_pylon_};
    TS_ = {a_T0};
  } else {
    // ensures that ropeway_inc satisfies the cable equation
    ropeway_inc = get_adm_1C(meca, {left_pylon_, right_pylon_});

    nodes_coords_.resize(nb_nodes - 1);
    TS_.resize(nb_nodes);
    internal_forces_.resize(nb_nodes);

    // computes the initial profile of the
    // --> updates nodes_coords, internal_forces and TS
    get_profile_1C(meca, ropeway_inc, nb_nodes, nodes_coords_, internal_forces_, TS_);

    // updates support reaction
    Point drk;
    drk.diff(nodes_coords_[1], nodes_coords_[0]);
    double dot = drk.dot(internal_forces_[0]);
    if (dot < 0) {
      for (auto &r : internal_forces_) {
        r.opposite();
      }
    } else if (dot == 0) {
      for (auto &r : internal_forces_) {
        r.zero();
      }
    }
    support_reaction_.diff(a_R0, internal_forces_[0]);
  }
}

int siconos::fem::cable::Rope::computeNbNodes(int nb_elem, double L) {
  number_of_nodes_ = 0;
  if (!m_last_) {
    number_of_nodes_ = (int)rint(nb_elem * get_L() / L);
  }
  return number_of_nodes_;
}

int siconos::fem::cable::Rope::computeMesh(std::vector<Point> &a_q, std::vector<Point> &a_R,
                                           std::vector<double> &a_TS, int q_offset,
                                           bool a_reverse) {
  if (!m_last_) {
    get_profile_1C(meca, ropeway_inc, number_of_nodes_ + 1, a_q, a_R, a_TS, q_offset,
                   a_reverse);
  }
  return number_of_nodes_;
}

double siconos::fem::cable::Rope::get_T0() {
  if (TS_.size())
    return TS_.front();
  else
    return 0.0;
}

double siconos::fem::cable::Rope::get_LastT() {
  if (TS_.size())
    return TS_.back();
  else
    return 0.0;
}

siconos::fem::cable::Point siconos::fem::cable::Rope::get_LastR() {
  if (internal_forces_.size())
    return internal_forces_.back();
  else
    return Point();
}

double siconos::fem::cable::Rope::get_L() { return ropeway_inc.x; }

const siconos::fem::cable::Point &siconos::fem::cable::Rope::supportReaction() const {
  return support_reaction_;
}

const siconos::fem::cable::MechanicalProperties &siconos::fem::cable::Rope::get_meca() const {
  return meca;
}

const siconos::fem::cable::Pylon &siconos::fem::cable::Rope::left_pylon() const {
  return left_pylon_;
}

void siconos::fem::cable::Rope::to_json(ojson &j) {
  if (TS_.size()) {
    j["ropeway_inc"].push_back(ropeway_inc);

    for (auto &elem : nodes_coords_) {
      j["q"].push_back(elem);
    }
    for (auto &elem : TS_) {
      j["TS"].push_back(elem);
    }
    for (auto &elem : internal_forces_) {
      j["R"].push_back(elem);
    }

    j["SR"].push_back(support_reaction_);

    Point vMeca(meca.get_T0(), meca.get_EA(), meca.get_rho());
    j["meca_global"].push_back(vMeca);
  } else {
    for (auto &elem : nodes_coords_) {
      j["q"].push_back(elem);
    }
    j["SR"].push_back(support_reaction_);
  }
}

siconos::fem::cable::Point siconos::fem::cable::Rope::get_adm_1C(
    const MechanicalProperties &a_meca, const std::vector<Pylon> &bc) {
  // Initial guess for lenght and slopes
  Point cable_inc = guess(bc);

  bool test = true;
  int n = 0;
  Point residu;
  std::vector<std::vector<double>> jacobian;
  while (test) {
    n++;
    // compute initial residual eq and jacobian
    cable_eq(a_meca, bc, cable_inc, residu, jacobian);

    //    delta = np.linalg.solve(jacobian,r)
    Point delta;
    double det = jacobian[0][0] * jacobian[1][1] * jacobian[2][2] +
                 jacobian[1][0] * jacobian[2][1] * jacobian[0][2] +
                 jacobian[0][1] * jacobian[1][2] * jacobian[2][0] -
                 (jacobian[2][0] * jacobian[1][1] * jacobian[0][2] +
                  jacobian[1][0] * jacobian[0][1] * jacobian[2][2] +
                  jacobian[0][0] * jacobian[2][1] * jacobian[1][2]);
    if (det != 0) {
      std::vector<double> trans1{jacobian[0][0], jacobian[1][0], jacobian[2][0]};
      std::vector<double> trans2{jacobian[0][1], jacobian[1][1], jacobian[2][1]};
      std::vector<double> trans3{jacobian[0][2], jacobian[1][2], jacobian[2][2]};
      std::vector<std::vector<double>> trans{trans1, trans2, trans3};
      double m00 = trans[1][1] * trans[2][2] - trans[2][1] * trans[1][2];
      double m01 = trans[1][0] * trans[2][2] - trans[2][0] * trans[1][2];
      double m02 = trans[1][0] * trans[2][1] - trans[2][0] * trans[1][1];
      double m10 = trans[0][1] * trans[2][2] - trans[2][1] * trans[0][2];
      double m11 = trans[0][0] * trans[2][2] - trans[2][0] * trans[0][2];
      double m12 = trans[0][0] * trans[2][1] - trans[2][0] * trans[0][1];
      double m20 = trans[0][1] * trans[1][2] - trans[1][1] * trans[0][2];
      double m21 = trans[0][0] * trans[1][2] - trans[1][0] * trans[0][2];
      double m22 = trans[0][0] * trans[1][1] - trans[1][0] * trans[0][1];
      std::vector<double> inv1{m00 / det, -1 * m01 / det, m02 / det};
      std::vector<double> inv2{-1 * m10 / det, m11 / det, -1 * m12 / det};
      std::vector<double> inv3{m20 / det, -1 * m21 / det, m22 / det};
      std::vector<std::vector<double>> inverse{inv1, inv2, inv3};
      delta.x = inverse[0][0] * residu.x + inverse[0][1] * residu.y + inverse[0][2] * residu.z;
      delta.y = inverse[1][0] * residu.x + inverse[1][1] * residu.y + inverse[1][2] * residu.z;
      delta.z = inverse[2][0] * residu.x + inverse[2][1] * residu.y + inverse[2][2] * residu.z;
    }
    //    cable_inc = cable_inc - delta
    cable_inc.diff(cable_inc, delta);

    //    test = n<max_iter_ and np.linalg.norm(r)>tol_ and np.linalg.norm(delta)>tol_
    test = (n < max_iter_) && (residu.norm() > tol_) && (delta.norm() > tol_);
  }
  return cable_inc;
}

siconos::fem::cable::Point siconos::fem::cable::Rope::guess(const std::vector<Pylon> &bc) {
  Point res;
  Point delta;
  delta.diff(bc[1], bc[0]);
  // L = np.linalg.norm(pl-po)
  res.x = delta.norm();
  // v = (1/L)*(pl-po)
  res.y = (1 / res.x) * (delta.y);
  res.z = (1 / res.x) * (delta.z);
  return res;
}

void siconos::fem::cable::Rope::cable_eq(const MechanicalProperties &a_meca,
                                         const std::vector<Pylon> &bc, const Point &cable_inc,
                                         Point &residu,
                                         std::vector<std::vector<double>> &jacobian) {
  //  Returns the residual equation and the jacobian of a fixed-fixed cable with
  // imposed initial tension

  // po, pl = bc
  const Point &po = bc[0];
  const Point &pl = bc[1];
  // xo, yo, zo = po
  double xo = po.x;
  double yo = po.y;
  double zo = po.z;
  // xl, yl, zl = pl
  double xl = pl.x;
  double yl = pl.y;
  double zl = pl.z;
  // L, etaY, etaZ = cable_inc
  double L = cable_inc.x;
  double etaY = cable_inc.y;
  double etaZ = cable_inc.z;
  // sq0 = np.sqrt(1 + etaY**2 + etaZ**2)
  double sq0 = sqrt(1 + std::pow(etaZ, 2) + std::pow(etaY, 2));
  double alpha = a_meca.get_alpha();
  double beta = a_meca.get_beta();
  // sqL = np.sqrt(1 + etaY**2 + (etaZ + sq0*alpha*L)**2)
  double sqL = sqrt(1 + std::pow(etaY, 2) + std::pow((etaZ + sq0 * alpha * L), 2));
  // logL = np.log((etaZ + alpha*sq0*L + sqL)/(etaZ + sq0))
  double logL = log((etaZ + alpha * sq0 * L + sqL) / (etaZ + sq0));

  // Elastic catenary equations
  residu.x = (-((beta * L) / sq0) + xl - xo - logL / (alpha * sq0));
  residu.y = (-((beta * etaY * L) / sq0) + yl - yo - (etaY * logL) / (alpha * sq0));
  residu.z = (-((beta * etaZ * L) / sq0) - (alpha * beta * std::pow(L, 2)) / 2.0 +
              (sq0 - sqL) / (alpha * sq0) + zl - zo);

  // Elastic catenary jacobian
  double eq_x_L = (-(beta / sq0) - 1 / sqL);
  double eq_x_etaY =
      ((etaY * (beta * L -
                (etaZ * sq0 *
                 (sq0 + alpha * etaZ * L - alpha * sq0 * L +
                  std::pow(alpha, 2) * sq0 * std::pow(L, 2) - sqL + alpha * L * sqL)) /
                    (alpha * (etaZ + sq0) * sqL * (etaZ + alpha * sq0 * L + sqL)) +
                logL / alpha)) /
       std::pow(sq0, 3));
  double eq_x_etaZ = ((-1 - std::pow(etaZ, 2) - std::pow(etaY, 2) - alpha * etaZ * sq0 * L +
                       sq0 * sqL + alpha * beta * etaZ * L * sqL + etaZ * sqL * logL) /
                      (alpha * std::pow(sq0, 3) * sqL));

  double eq_y_L = (-((beta * etaY) / sq0) - etaY / sqL);
  double eq_y_etaY =
      ((beta * std::pow(etaY, 2) * L - beta * (1 + std::pow(etaZ, 2) + std::pow(etaY, 2)) * L -
        (etaZ * std::pow(etaY, 2) * sq0 *
         (sq0 + alpha * etaZ * L - alpha * sq0 * L +
          std::pow(alpha, 2) * sq0 * std::pow(L, 2) - sqL + alpha * L * sqL)) /
            (alpha * (etaZ + sq0) * sqL * (etaZ + alpha * sq0 * L + sqL)) +
        (std::pow(etaY, 2) * logL) / alpha -
        ((1 + std::pow(etaZ, 2) + std::pow(etaY, 2)) * logL) / alpha) /
       std::pow(sq0, 3));
  double eq_y_etaZ =
      (-((etaY * (1 + std::pow(etaZ, 2) + std::pow(etaY, 2) + alpha * etaZ * sq0 * L -
                  sq0 * sqL - alpha * beta * etaZ * L * sqL - etaZ * sqL * logL)) /
         (alpha * std::pow(sq0, 3) * sqL)));

  double eq_z_L = -((beta * etaZ) / sq0) - alpha * beta * L - (etaZ + alpha * sq0 * L) / sqL;
  double eq_z_etaY = ((etaZ * etaY * L * (sq0 + beta * sqL)) / (std::pow(sq0, 3) * sqL));
  double eq_z_etaZ =
      (beta * ((std::pow(etaZ, 2) * L) / std::pow(sq0, 3) - L / sq0) -
       ((1 + std::pow(etaY, 2)) * L) / ((1 + std::pow(etaZ, 2) + std::pow(etaY, 2)) * sqL));

  // jac = np.array([[eq_x_L,eq_x_etaY,eq_x_etaZ], [eq_y_L,eq_y_etaY,eq_y_etaZ],
  // [eq_z_L,eq_z_etaY,eq_z_etaZ]])
  std::vector<double> j1{eq_x_L, eq_x_etaY, eq_x_etaZ};
  std::vector<double> j2{eq_y_L, eq_y_etaY, eq_y_etaZ};
  std::vector<double> j3{eq_z_L, eq_z_etaY, eq_z_etaZ};
  jacobian = {j1, j2, j3};
}

void siconos::fem::cable::Rope::get_profile_1C(const MechanicalProperties &a_meca,
                                               const Point &cable_inc, int nb_nodes,
                                               std::vector<Point> &a_q,
                                               std::vector<Point> &a_R,
                                               std::vector<double> &a_TS, int q_offset,
                                               bool a_reverse) {
  /*
   @author: charl

  (cable_param,cable_inc,nb_nodes)

  where

  cable_param is a list containing [T0, EA, rho]
          T0  : initial tension
          EA  : rigidity of the cable
          rho : density of the cable

  bc is a np.array containing [[pox, poy, poz], [plx, ply, plz]]
          po  : first eyelet
          pl  : second eyelet

  cable_inc is a list containing [L, etaY, etaZ]
          L    : Unstretched length of the cable
          etaY : Initial y-slope
          etaZ : Initial z-slope

  nb_nodes is the number of nodes in the system

  Returns the positions in the shape of a 3*nb_nodes vector [x,y,z,...,x,y,z],
  the tension in the shape of a nb_nodes vector [T,...T]
  and the internal force vector in the shape of a 3*nb_nodes vector [H,V,B,...,H,V,B]
  */
  double xl = right_pylon_.x;
  double yl = right_pylon_.y;
  double zl = right_pylon_.z;
  double L = cable_inc.x;
  double etaY = cable_inc.y;
  double etaZ = cable_inc.z;
  double sq0 = sqrt(1 + std::pow(etaY, 2) + std::pow(etaZ, 2));
  double alpha = a_meca.get_alpha();
  double beta = a_meca.get_beta();
  double sqL = sqrt(1 + std::pow(etaY, 2) + std::pow((etaZ + sq0 * alpha * L), 2));
  double coef = a_meca.get_T0() / sq0;

  for (int i = 0; i < nb_nodes; i++) {
    int j = (a_reverse) ? nb_nodes - 1 - i : i;
    double S = L / (nb_nodes - 1) * j;
    double sqS = sqrt(1 + std::pow(etaY, 2) + std::pow((etaZ + sq0 * alpha * S), 2));
    double logS = log((etaZ + alpha * sq0 * L + sqL) / (etaZ + alpha * sq0 * S + sqS));

    // Elastic catenary positions
    if (i < nb_nodes - 1) {
      Point &p = a_q[i + q_offset];
      p.x = xl - beta / sq0 * (L - S) - 1 / (alpha * sq0) * logS;
      p.y = yl - etaY * beta / sq0 * (L - S) - etaY / (alpha * sq0) * logS;
      p.z = zl - etaZ * beta / sq0 * (L - S) -
            0.5 * alpha * beta * (std::pow(L, 2) - std::pow(S, 2)) -
            1 / (sq0 * alpha) * (sqL - sqS);
    }
    // Tension
    a_TS[i + q_offset] = coef * sqS;

    // internal forces
    Point &r = a_R[i + q_offset];
    r.x = coef;
    r.y = coef * etaY;
    r.z = coef * (etaZ + alpha * sq0 * S);
  }
}

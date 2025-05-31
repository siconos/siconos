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

#include "CableDS.hpp"

#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

siconos::fem::cable::CableDS::CableDS(Eigen::Ref<siconos::algebra::SiconosVector> q0,
                                      Eigen::Ref<siconos::algebra::SiconosVector> velocity0,
                                      Eigen::Ref<siconos::algebra::SiconosMatrix> mass,
                                      double a_EA, double a_elem_length)
    : LagrangianDS(q0, velocity0), EA_{a_EA}, l_e_{a_elem_length} {
  std::cout << " BUlD CABLE DS \n";
  setConstantMass(mass);

  TRNp_Np = TRNp_NpMatrix();

  siconos::algebra::print(*TRNp_Np);

  // Constructor with initial state and mass.
  // We assume that q0, v0 and mass are computed by the "cable model" based on mesh and other
  // inputs Maybe these things are to be computed inside the DS? To be discussed ...
  // The call to LagrangianDS base constructor  ensures a proper allocation of memories for q0,
  // v0 Mass is just a pointer link. Mass alloc : to be done in cable model if mass is a shared
  // pointer input. _mass =
  // std::make_shared<siconos::algebra::Matrix>(ndof_, ndof_,
  // UBLAS_TYPE::SPARSE);
  // We can deal with variable mass later.

  // What may happen in cable model: (see examples in testCableDS)
  // Case 1: no fext
  // cableDS{q0, v0, mass}

  // Case 2: constant fext
  // cableDS{q0, v0, mass}
  // call to setConstantFext()
  // cable

  totalForces_ = std::make_shared<siconos::algebra::SiconosVector>(ndof_);
  totalForces_->setZero();
  // ndof_ is given by the size of q0

  // We will use jacobianTotalForcesOver_q_ and _jacobianTotalForcesOver_velocity to save
  // tangent stiffness and damping matrices. Those are attributes of LagrangianDS class.
  jacobianTotalForcesOver_q_ = std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
  jacobianTotalForcesOver_q_->setZero();

  jacobianTotalForcesOver_velocity_ =
      std::make_shared<siconos::algebra::SiconosMatrix>(ndof_, ndof_);
  jacobianTotalForcesOver_velocity_->setZero();
}

void siconos::fem::cable::CableDS::computeTotalForces(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time) {
  // Update internal and external forces
  // Here you must:
  //  - compute internal forces and save them somewhere or directly in _forces
  //  - assume that jacobian are uptodate or maybe recompute them
  // *_forces = ...
  tangentStiffnessMatrix(q);

  //  - compute external forces in fExt (or just get them if they are constant)

  computeFext(time);
  if (fext_view_) *totalForces_ += *fext_view_;
}

// \f$ \nabla_q F \f$
void siconos::fem::cable::CableDS::computeJacobianTotalForcesOver_q(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time) {
  // Call a local routine which compute tangent stiffness and update a local operator

  // tangentStiffnessMatrix();
}

// \f$ \nabla_v F \f$
void siconos::fem::cable::CableDS::computeJacobianTotalForcesOver_velocity(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &velocity,
    const Eigen::Ref<const siconos::algebra::SiconosVector> &q, double time) {
  // Call a local routine to compute damping matrix and update a local operator
  dampingMatrix();
}

void siconos::fem::cable::CableDS::tangentStiffnessMatrix(
    const Eigen::Ref<const siconos::algebra::SiconosVector> q) {
  // must update jacobianTotalForcesOver_q
  // KT
  // intforces: K

  size_t nb_elem = ndof_ - 3;
  double k = EA_ / l_e_;
  siconos::algebra::SiconosVector Tq{6};
  siconos::algebra::SiconosMatrix TqqT{6, 6};

  // tous les points moins le dernier
  for (size_t i = 0; i < nb_elem - 3; i += 3) {
    double n_e = 0;
    for (size_t j = 0; j < 3; j++) {
      double d = (q(i + j) - q(i + 3 + j));
      n_e += d * d;
    }
    n_e = sqrt(n_e);

    double eps = n_e / l_e_ - 1;
    double f_e = fabs(eps);
    double kf = k / (1 + 1 / f_e);
    if (eps > 0) {
      // fi
      matmult(q, i, Tq);
      for (size_t j = 0; j < 6; j++) {
        (*totalForces_)(i + j) += kf * Tq(j);
      }

      // KT
      matmult2(Tq, TqqT);
      double kKT = 1 / (1 + f_e);
      kKT = k * kKT * kKT * (1 / (l_e_ * n_e));
      for (size_t j = 0; j < 6; j++) {
        for (size_t l = 0; l < 6; l++) {
          auto val = (*jacobianTotalForcesOver_q_)(i + j, i + l) + kKT * TqqT(j, l);
          jacobianTotalForcesOver_q_->setValue(i + j, i + l, val);
        }
      }
    }
    for (size_t j = 0; j < 6; j++) {
      for (size_t l = 0; l < 6; l++) {
        auto val = (*jacobianTotalForcesOver_q_)(i + j, i + l) + kf * (*TRNp_Np)(j, l);
        jacobianTotalForcesOver_q_->setValue(i + j, i + l, val);
      }
    }
  }
  // dernier point, à voir si on mets dans une seule boucle avec des if
  // dernier élément - premier élément
  double n_e = 0;
  for (size_t j = 0; j < 3; j++) {
    double d = (q(nb_elem + j) - q(j));
    n_e += d * d;
  }
  n_e = sqrt(n_e);

  double eps = n_e / l_e_ - 1;
  double f_e = fabs(eps);
  double kf = k / (1 + 1 / f_e);
  if (eps > 0) {
    // fi
    matmult(q, nb_elem, Tq);
    for (size_t j = 0; j < 3; j++) {
      (*totalForces_)(nb_elem + j) += kf * Tq(j);
      (*totalForces_)(j) += kf * Tq(j + 3);
    }

    // KT
    matmult2(Tq, TqqT);
    double kKT = 1 / (1 + f_e);
    kKT = k * kKT * kKT * (1 / (l_e_ * n_e));
    for (size_t j = 0; j < 3; j++) {
      for (size_t l = 0; l < 3; l++) {
        auto val = (*jacobianTotalForcesOver_q_)(nb_elem + j, nb_elem + l) + kKT * TqqT(j, l);
        jacobianTotalForcesOver_q_->setValue(nb_elem + j, nb_elem + l, val);
        val = (*jacobianTotalForcesOver_q_)(nb_elem + j, l) + kKT * TqqT(j, l + 3);
        jacobianTotalForcesOver_q_->setValue(nb_elem + j, l, val);
        val = (*jacobianTotalForcesOver_q_)(j, nb_elem + l) + kKT * TqqT(j + 3, l);
        jacobianTotalForcesOver_q_->setValue(j, nb_elem + l, val);
        val = (*jacobianTotalForcesOver_q_)(j, l) + kKT * TqqT(j + 3, l + 3);
        jacobianTotalForcesOver_q_->setValue(j, l, val);
      }
    }
  }
  for (size_t j = 0; j < 3; j++) {
    for (size_t l = 0; l < 3; l++) {
      auto val =
          (*jacobianTotalForcesOver_q_)(nb_elem + j, nb_elem + l) + kf * (*TRNp_Np)(j, l);
      jacobianTotalForcesOver_q_->setValue(nb_elem + j, nb_elem + l, val);
      val = (*jacobianTotalForcesOver_q_)(nb_elem + j, l) + kf * (*TRNp_Np)(j, l + 3);
      jacobianTotalForcesOver_q_->setValue(nb_elem + j, l, val);
      val = (*jacobianTotalForcesOver_q_)(j, nb_elem + l) + kf * (*TRNp_Np)(j + 3, l);
      jacobianTotalForcesOver_q_->setValue(j, nb_elem + l, val);
      val = (*jacobianTotalForcesOver_q_)(j, l) + kf * (*TRNp_Np)(j + 3, l + 3);
      jacobianTotalForcesOver_q_->setValue(j, l, val);
    }
  }
}

void siconos::fem::cable::CableDS::dampingMatrix(/** ...*/) {
  // must update jacobianqDotForces
  // constant ?
  // C = damp*M
}

void siconos::fem::cable::CableDS::matmult(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &V,
    siconos::algebra::SiconosSize_t a_startIdx,
    Eigen::Ref<siconos::algebra::SiconosVector> R) {
  R.setZero();
  assert(TRNp_Np);
  auto n = R.size();
  if (n + a_startIdx < V.size()) {
    for (auto i = 0; i < n; i++) {
      for (auto j = 0; j < n; j++) {
        auto val = (*TRNp_Np)(i, j);
        R(i) += val * V(j + a_startIdx);
      }
    }
  } else {
    for (auto i = 0; i < n; i++) {
      for (auto j = 0; j < 3; j++) {
        auto val = (*TRNp_Np)(i, j);
        R(i) += val * V(j + a_startIdx);
      }
      for (auto j = 3; j < 6; j++) {
        auto val = (*TRNp_Np)(i, j);
        R(i) += val * V(j - 3);
      }
    }
  }
}

void siconos::fem::cable::CableDS::matmult2(
    const Eigen::Ref<const siconos::algebra::SiconosVector> &V,
    Eigen::Ref<siconos::algebra::SiconosMatrix> R) {
  size_t n = V.size();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      R(i, j) = V(i) * V(j);
    }
  }
}
std::shared_ptr<siconos::algebra::SiconosMatrix>
siconos::fem::cable::CableDS::TRNp_NpMatrix() {
  auto vTRNp_Np =
      std::make_shared<siconos::algebra::SiconosMatrix>(6, 6);  // FP: HAS TO BE SPARSE
  /* vector<vector<double>> TRNp_Np = {{1, 0, 0, -1, 0, 0},
                                                                         {0, 1, 0, 0, -1, 0},
                                     {0, 0, 1, 0, 0, -1},
                                                                         {-1, 0, 0, 1, 0, 0},
                                     {0, -1, 0, 0, 1, 0},
                                                                         {0, 0, -1, 0, 0,
     1}};
                                                                         */
  vTRNp_Np->setIdentity();
  (*vTRNp_Np)(3, 1) = -1;
  (*vTRNp_Np)(4, 2) = -1;
  (*vTRNp_Np)(5, 3) = -1;

  return vTRNp_Np;
}

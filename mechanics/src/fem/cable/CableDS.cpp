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

siconos::fem::cable::CableDS::CableDS(
    std::shared_ptr<siconos::algebra::SiconosVector> q0,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity0, std::shared_ptr<Matrix> mass,
    double a_EA, double a_elem_length, ExternalForcesFunction fext)
    : LagrangianDS(q0, velocity0, mass), computefext_{fext}, _EA{a_EA}, _l_e{a_elem_length} {
  std::cout << " BUlD CABLE DS \n";

  TRNp_Np = TRNp_NpMatrix();

  TRNp_Np->display();

  // Constructor with initial state and mass.
  // We assume that q0, v0 and mass are computed by the "cable model" based on mesh and other
  // inputs Maybe these things are to be computed inside the DS? To be discussed ...
  // The call to LagrangianDS base constructor  ensures a proper allocation of memories for q0,
  // v0 Mass is just a pointer link. Mass alloc : to be done in cable model if mass is a shared
  // pointer input. _mass =
  // std::make_shared<siconos::algebra::Matrix>(_ndof, _ndof,
  // UBLAS_TYPE::SPARSE);
  // We can deal with variable mass later.
  _hasConstantMass = true;

  // What may happen in cable model: (see examples in testCableDS)
  // Case 1: no fext
  // cableDS{q0, v0, mass}

  // Case 2: constant fext
  // cableDS{q0, v0, mass}
  // call to setFExtPtr()
  // cable

  if (fext) {
    _hasConstantFExt = false;  // Indeed, this is the default for SecondOrderDS
    _fExt = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  } else
    _hasConstantFExt = true;
  // In that case _fExt = nullptr
  // setFextPtr is to be called later by cable model to set a constant fext.

  // _ndof is given by the size of q0 during SecondOrderDS build
  _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);

  // We will use _jacobianqForces and _jacobianvForces to save tangent stiffness and damping
  // matrices.
  // Those are attributes of LagrangianDS class.

  _jacobianqForces = std::make_shared<Matrix>(_ndof, _ndof);  // FP: has to be SPARSE !!

  _jacobianqDotForces = std::make_shared<Matrix>(_ndof, _ndof);  // FP: has to be SPARSE !!
}

void siconos::fem::cable::CableDS::computeFExt(double time) {
  assert(_fExt);
  assert(computefext_);
  // Call the std::function attribute that must be connected to some external function
  computefext_(time, _fExt);
  // ...
}

void siconos::fem::cable::CableDS::computeForces(
    double time, std::shared_ptr<siconos::algebra::SiconosVector> q,
    std::shared_ptr<siconos::algebra::SiconosVector> velocity) {
  assert(_forces);
  // if (!_forces) {
  //   _forces = std::make_shared<siconos::algebra::SiconosVector>(_ndof);
  // } // --> done during constructor call.
  // else
  _forces->zero();

  // Update internal and external forces
  // Here you must:
  //  - compute internal forces and save them somewhere or directly in _forces
  //  - assume that jacobian are uptodate or maybe recompute them
  // *_forces = ...
  tangentStiffnessMatrix(q);

  //  - compute external forces in fExt (or just get them if they are constant)
  if (!_hasConstantFExt) computeFExt(time);
  if (_fExt) *_forces += *_fExt;
}

// \f$ \nabla_q F \f$
void siconos::fem::cable::CableDS::computeJacobianqForces(double time) {
  // Call a local routine which compute tangent stiffness and update a local operator

  // tangentStiffnessMatrix();
}

// \f$ \nabla_v F \f$
void siconos::fem::cable::CableDS::computeJacobianvForces(double time) {
  // Call a local routine to compute damping matrix and update a local operator
  dampingMatrix();
}

void siconos::fem::cable::CableDS::tangentStiffnessMatrix(
    std::shared_ptr<siconos::algebra::SiconosVector> q) {
  // must update jacobianqForces
  // KT
  // intforces: K

  size_t nb_elem = _ndof - 3;
  double k = _EA / _l_e;
  auto Tq = std::make_shared<siconos::algebra::SiconosVector>(6);
  auto TqqT = std::make_shared<Matrix>(6, 6);

  // tous les points moins le dernier
  for (size_t i = 0; i < nb_elem - 3; i += 3) {
    double n_e = 0;
    for (size_t j = 0; j < 3; j++) {
      double d = ((*q)(i + j) - (*q)(i + 3 + j));
      n_e += d * d;
    }
    n_e = sqrt(n_e);

    double eps = n_e / _l_e - 1;
    double f_e = fabs(eps);
    double kf = k / (1 + 1 / f_e);
    if (eps > 0) {
      // fi
      matmult(q, i, Tq);
      for (size_t j = 0; j < 6; j++) {
        (*_forces)(i + j) += kf * (*Tq)(j);
      }

      // KT
      matmult2(Tq, TqqT);
      double kKT = 1 / (1 + f_e);
      kKT = k * kKT * kKT * (1 / (_l_e * n_e));
      for (size_t j = 0; j < 6; j++) {
        for (size_t l = 0; l < 6; l++) {
          auto val = _jacobianqForces->getValue(i + j, i + l) + kKT * (*TqqT)(j, l);
          _jacobianqForces->setValue(i + j, i + l, val);
        }
      }
    }
    for (size_t j = 0; j < 6; j++) {
      for (size_t l = 0; l < 6; l++) {
        auto val = _jacobianqForces->getValue(i + j, i + l) + kf * TRNp_Np->getValue(j, l);
        _jacobianqForces->setValue(i + j, i + l, val);
      }
    }
  }
  // dernier point, à voir si on mets dans une seule boucle avec des if
  // dernier élément - premier élément
  double n_e = 0;
  for (size_t j = 0; j < 3; j++) {
    double d = ((*q)(nb_elem + j) - (*q)(j));
    n_e += d * d;
  }
  n_e = sqrt(n_e);

  double eps = n_e / _l_e - 1;
  double f_e = fabs(eps);
  double kf = k / (1 + 1 / f_e);
  if (eps > 0) {
    // fi
    matmult(q, nb_elem, Tq);
    for (size_t j = 0; j < 3; j++) {
      (*_forces)(nb_elem + j) += kf * (*Tq)(j);
      (*_forces)(j) += kf * (*Tq)(j + 3);
    }

    // KT
    matmult2(Tq, TqqT);
    double kKT = 1 / (1 + f_e);
    kKT = k * kKT * kKT * (1 / (_l_e * n_e));
    for (size_t j = 0; j < 3; j++) {
      for (size_t l = 0; l < 3; l++) {
        auto val = _jacobianqForces->getValue(nb_elem + j, nb_elem + l) + kKT * (*TqqT)(j, l);
        _jacobianqForces->setValue(nb_elem + j, nb_elem + l, val);
        val = _jacobianqForces->getValue(nb_elem + j, l) + kKT * (*TqqT)(j, l + 3);
        _jacobianqForces->setValue(nb_elem + j, l, val);
        val = _jacobianqForces->getValue(j, nb_elem + l) + kKT * (*TqqT)(j + 3, l);
        _jacobianqForces->setValue(j, nb_elem + l, val);
        val = _jacobianqForces->getValue(j, l) + kKT * (*TqqT)(j + 3, l + 3);
        _jacobianqForces->setValue(j, l, val);
      }
    }
  }
  for (size_t j = 0; j < 3; j++) {
    for (size_t l = 0; l < 3; l++) {
      auto val =
          _jacobianqForces->getValue(nb_elem + j, nb_elem + l) + kf * TRNp_Np->getValue(j, l);
      _jacobianqForces->setValue(nb_elem + j, nb_elem + l, val);
      val = _jacobianqForces->getValue(nb_elem + j, l) + kf * TRNp_Np->getValue(j, l + 3);
      _jacobianqForces->setValue(nb_elem + j, l, val);
      val = _jacobianqForces->getValue(j, nb_elem + l) + kf * TRNp_Np->getValue(j + 3, l);
      _jacobianqForces->setValue(j, nb_elem + l, val);
      val = _jacobianqForces->getValue(j, l) + kf * TRNp_Np->getValue(j + 3, l + 3);
      _jacobianqForces->setValue(j, l, val);
    }
  }
}

void siconos::fem::cable::CableDS::dampingMatrix(/** ...*/) {
  // must update jacobianqDotForces
  // constant ?
  // C = damp*M
}

void siconos::fem::cable::CableDS::matmult(
    const std::shared_ptr<siconos::algebra::SiconosVector> &V, size_t a_startIdx,
    std::shared_ptr<siconos::algebra::SiconosVector> &R) {
  R->zero();
  assert(TRNp_Np);
  auto n = R->size();
  if (n + a_startIdx < V->size()) {
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
        auto val = TRNp_Np->getValue(i, j);
        (*R)(i) += val * (*V)(j + a_startIdx);
      }
    }
  } else {
    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < 3; j++) {
        auto val = TRNp_Np->getValue(i, j);
        (*R)(i) += val * (*V)(j + a_startIdx);
      }
      for (size_t j = 3; j < 6; j++) {
        auto val = TRNp_Np->getValue(i, j);
        (*R)(i) += val * (*V)(j - 3);
      }
    }
  }
}

void siconos::fem::cable::CableDS::matmult2(
    const std::shared_ptr<siconos::algebra::SiconosVector> &V, std::shared_ptr<Matrix> &R) {
  size_t n = V->size();
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      (*R)(i, j) = (*V)(i) * (*V)(j);
    }
  }
}
std::shared_ptr<siconos::modeling::SecondOrderDS::Matrix>
siconos::fem::cable::CableDS::TRNp_NpMatrix() {
  auto vTRNp_Np = std::make_shared<Matrix>(6, 6);  // FP: HAS TO BE SPARSE
  /* vector<vector<double>> TRNp_Np = {{1, 0, 0, -1, 0, 0},
                                                                         {0, 1, 0, 0, -1, 0},
                                     {0, 0, 1, 0, 0, -1},
                                                                         {-1, 0, 0, 1, 0, 0},
                                     {0, -1, 0, 0, 1, 0},
                                                                         {0, 0, -1, 0, 0, 1}};
                                                                         */
  vTRNp_Np->eye();
  vTRNp_Np->setValue(3, 1, -1);
  vTRNp_Np->setValue(4, 2, -1);
  vTRNp_Np->setValue(5, 3, -1);

  return vTRNp_Np;
}

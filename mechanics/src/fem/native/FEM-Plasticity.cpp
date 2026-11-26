// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "FEM-Plasticity.hpp"

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <algorithm>
#include <cstdio>
#include <iostream>  // DEBUG_EXPR
#include <memory>
#include <string>
#include <vector>

#include "CSparseMatrix.h"
#include "DruckerPragerPlasticityNSL.hpp"
#include "Interaction.hpp"
#include "NonSmoothDrivers.h"  // plasticity_2d_driver
#include "NumericsMatrix.h"
#include "OSNSMatrix.hpp"
#include "PlasticityProblem.h"
#include "Simulation.hpp"
#include "SimulationGraphs.hpp"
#include "SolverOptions.h"
#include "StressLinearTIR.hpp"
#include "VonMisesPlasticityNSL.hpp"
#include "siconos_debug.h"

namespace {
using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor, int>;
using Vec = Eigen::VectorXd;

/* NumericsMatrix (NM_SPARSE) -> Eigen sparse */
SpMat to_eigen(NumericsMatrix* nm) {
  CSparseMatrix* csc = NM_csc(nm);
  std::vector<Eigen::Triplet<double, int>> triplets;
  triplets.reserve(static_cast<size_t>(csc->p[csc->n]));
  for (CS_INT j = 0; j < csc->n; j++)
    for (CS_INT k = csc->p[j]; k < csc->p[j + 1]; k++)
      triplets.emplace_back(static_cast<int>(csc->i[k]), static_cast<int>(j), csc->x[k]);
  SpMat A(static_cast<int>(csc->m), static_cast<int>(csc->n));
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}

/* selection matrix (rows idx of an n-vector) */
SpMat selection(const std::vector<int>& idx, int n) {
  std::vector<Eigen::Triplet<double, int>> t;
  t.reserve(idx.size());
  for (size_t k = 0; k < idx.size(); k++) t.emplace_back(static_cast<int>(k), idx[k], 1.0);
  SpMat P(static_cast<int>(idx.size()), n);
  P.setFromTriplets(t.begin(), t.end());
  return P;
}

/* true if A is square diagonal with a non-zero diagonal */
bool is_invertible_diagonal(const SpMat& A) {
  if (A.rows() != A.cols()) return false;
  std::vector<bool> has_diag(static_cast<size_t>(A.rows()), false);
  for (int j = 0; j < A.outerSize(); ++j)
    for (SpMat::InnerIterator it(A, j); it; ++it) {
      if (it.row() != it.col()) {
        if (it.value() != 0.0) return false;
      } else if (it.value() != 0.0)
        has_diag[static_cast<size_t>(j)] = true;
    }
  for (bool d : has_diag)
    if (!d) return false;
  return true;
}
/* inverse of H_c, a permutation of a block-diagonal matrix with 3x3 blocks (each group of 3
   columns, one interaction, has its non-zeros in 3 rows not shared with the other groups):
   inverted block by block. Otherwise (general C), sparse LU with an identity right-hand
   side (expensive: dense columns). */
SpMat block_inverse(const SpMat& Hc) {
  int m = static_cast<int>(Hc.cols());
  std::vector<Eigen::Triplet<double, int>> t;
  t.reserve(static_cast<size_t>(3 * m));
  bool blocks = (m % 3 == 0);
  for (int k = 0; blocks && k < m / 3; k++) {
    std::vector<int> rows;
    for (int j = 3 * k; j < 3 * k + 3; j++)
      for (SpMat::InnerIterator it(Hc, j); it; ++it)
        if (it.value() != 0.0 && std::find(rows.begin(), rows.end(), it.row()) == rows.end())
          rows.push_back(static_cast<int>(it.row()));
    if (rows.size() != 3) {
      blocks = false;
      break;
    }
    Eigen::Matrix3d B = Eigen::Matrix3d::Zero();
    for (int j = 0; j < 3; j++)
      for (SpMat::InnerIterator it(Hc, 3 * k + j); it; ++it)
        for (int i = 0; i < 3; i++)
          if (it.row() == rows[static_cast<size_t>(i)]) B(i, j) = it.value();
    Eigen::FullPivLU<Eigen::Matrix3d> lu(B);
    if (!lu.isInvertible())
      THROW_EXCEPTION("fem::nonsmooth_formulations::Plasticity: H_c is not invertible");
    Eigen::Matrix3d Binv = lu.inverse();
    // Hc(rows, cols_k) = B  =>  Hc^{-1}(cols_k, rows) = B^{-1}
    for (int a = 0; a < 3; a++)
      for (int i = 0; i < 3; i++)
        if (Binv(a, i) != 0.0)
          t.emplace_back(3 * k + a, rows[static_cast<size_t>(i)], Binv(a, i));
  }
  if (blocks) {
    SpMat Hinv(m, m);
    Hinv.setFromTriplets(t.begin(), t.end());
    return Hinv;
  }
  Eigen::SparseLU<SpMat> lu_h(Hc);
  if (lu_h.info() != Eigen::Success)
    THROW_EXCEPTION("fem::nonsmooth_formulations::Plasticity: H_c is not invertible");
  SpMat I(m, m);
  I.setIdentity();
  SpMat Hinv = lu_h.solve(I);
  Hinv.prune(0.0);
  return Hinv;
}
}  // namespace

siconos::mechanics::fem::nonsmooth_formulations::Plasticity::Plasticity(int numericsSolverId)
    : GlobalFrictionContact(
          3, std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
                                            solver_options_delete)) {}

int siconos::mechanics::fem::nonsmooth_formulations::Plasticity::compute(double time) {
  if (!preCompute(time)) return 0;
  int info = solve_condensed();
  postCompute();
  return info;
}

int siconos::mechanics::fem::nonsmooth_formulations::Plasticity::solve_condensed() {
  auto& indexSet = *simulation()->indexSet(indexSetLevel());

  // --- Interactions: law and its parameters (one model for all the interactions)
  auto m_size = static_cast<int>(_sizeOutput);
  enum class Model { none, von_mises, drucker_prager } model = Model::none;
  std::vector<double> sigma_y, eta, theta, cohesion;  // per interaction
  std::vector<int> column_of_cone;  // first column of H (position in z) of each cone
  auto hypothesis = siconos::modeling::VonMisesPlasticityNSL::Hypothesis::plane_stress;
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui) {
    auto inter = indexSet.bundle(*ui);
    if (inter->dimension() != 3)
      THROW_EXCEPTION("fem::nonsmooth_formulations::Plasticity: interactions of size 3 only");
    if (inter->relation()->getSubType() != siconos::modeling::RelationSubType::StressLinearTIR)
      THROW_EXCEPTION(
          "fem::nonsmooth_formulations::Plasticity: only StressLinearTIR relations are "
          "supported");
    auto vm = std::dynamic_pointer_cast<siconos::modeling::VonMisesPlasticityNSL>(
        inter->nonSmoothLaw());
    auto dp = std::dynamic_pointer_cast<siconos::modeling::DruckerPragerPlasticityNSL>(
        inter->nonSmoothLaw());
    Model m = vm ? Model::von_mises : dp ? Model::drucker_prager : Model::none;
    if (m == Model::none)
      THROW_EXCEPTION(
          "fem::nonsmooth_formulations::Plasticity: only VonMisesPlasticityNSL and "
          "DruckerPragerPlasticityNSL are supported");
    if (model == Model::none) {
      model = m;
      if (vm) hypothesis = vm->hypothesis();
    } else if (m != model || (vm && vm->hypothesis() != hypothesis))
      THROW_EXCEPTION(
          "fem::nonsmooth_formulations::Plasticity: all the interactions must have the same "
          "model (and the same Von Mises hypothesis)");
    if (vm) {
      sigma_y.push_back(vm->sigma_y());
    } else {
      eta.push_back(dp->eta());
      theta.push_back(dp->theta());
      cohesion.push_back(dp->cohesion());
    }
    column_of_cone.push_back(static_cast<int>(indexSet.properties(*ui).absolute_position));
  }
  int nc = static_cast<int>(column_of_cone.size());
  if (3 * nc != m_size)
    THROW_EXCEPTION("fem::nonsmooth_formulations::Plasticity: inconsistent problem size");
  if (nc == 0) {
    _globalVelocities->setZero();
    return 0;
  }

  // --- Global problem W x = q + H z, y = H^T x + b
  SpMat W = to_eigen(_W->numericsMatrix().get());
  SpMat H = to_eigen(_H->numericsMatrix().get());
  int n = static_cast<int>(W.rows());
  Vec q = *_q;
  Vec b = _b->head(m_size);

  // rows of x touched by H (constrained) and the others
  std::vector<char> constrained(static_cast<size_t>(n), 0);
  for (int j = 0; j < H.outerSize(); ++j)
    for (SpMat::InnerIterator it(H, j); it; ++it)
      if (it.value() != 0.0) constrained[static_cast<size_t>(it.row())] = 1;
  std::vector<int> c_idx, r_idx;
  for (int i = 0; i < n; i++)
    (constrained[static_cast<size_t>(i)] ? c_idx : r_idx).push_back(i);
  if (static_cast<int>(c_idx.size()) != m_size)
    THROW_EXCEPTION(
        "fem::nonsmooth_formulations::Plasticity: the plasticity interactions must constrain "
        "as many unknowns as their dimension (one interaction per stress block): " +
        std::to_string(c_idx.size()) + " constrained unknowns for " + std::to_string(m_size) +
        " interaction components");

  SpMat Pc = selection(c_idx, n), Pr = selection(r_idx, n);
  SpMat Wcc = Pc * W * Pc.transpose(), Wcr = Pc * W * Pr.transpose();
  SpMat Wrc = Pr * W * Pc.transpose(), Wrr = Pr * W * Pr.transpose();
  SpMat Hc = Pc * H;
  Vec qc = Pc * q, qr = Pr * q;

  // --- Schur complement on the constrained stresses
  SpMat Wrr_inv_Wrc;
  Vec Wrr_inv_qr;
  Eigen::SparseLU<SpMat> lu_rr;
  bool diagonal_rr = is_invertible_diagonal(Wrr);
  if (diagonal_rr) {
    Vec dinv = Wrr.diagonal().cwiseInverse();
    Wrr_inv_Wrc = dinv.asDiagonal() * Wrc;
    Wrr_inv_qr = dinv.cwiseProduct(qr);
  } else if (Wrr.rows() > 0) {
    lu_rr.compute(Wrr);
    if (lu_rr.info() != Eigen::Success)
      THROW_EXCEPTION("fem::nonsmooth_formulations::Plasticity: LU of W_rr failed");
    Wrr_inv_Wrc = lu_rr.solve(Wrc);
    Wrr_inv_qr = lu_rr.solve(qr);
  } else {
    Wrr_inv_Wrc = SpMat(0, m_size);
    Wrr_inv_qr = Vec(0);
  }
  SpMat Wred = Wcc - Wcr * Wrr_inv_Wrc;
  Vec qt = qc - Wcr * Wrr_inv_qr;

  // --- Change of variables x_c = Hc^{-T} (y - b)
  SpMat Hc_inv = block_inverse(Hc);
  SpMat Mp = Hc_inv * Wred * SpMat(Hc_inv.transpose());
  Mp.prune(0.0);
  Vec qp = -(Mp * b) - Hc_inv * qt;

  // --- numerics problem 0 in Mp y + qp + N_K(y), y in Voigt coordinates
  NumericsMatrix* nm = NM_create(NM_SPARSE, m_size, m_size);
  NM_triplet_alloc(nm, Mp.nonZeros());
  for (int j = 0; j < Mp.outerSize(); ++j)
    for (SpMat::InnerIterator it(Mp, j); it; ++it)
      NM_zentry(nm, static_cast<int>(it.row()), j, it.value(), 0.0);

  // cones are ordered as the columns of H (z): reorder the parameters accordingly
  auto by_column = [&](const std::vector<double>& values) {
    std::vector<double> ordered(values.size());
    for (int k = 0; k < static_cast<int>(values.size()); k++) {
      if (column_of_cone[k] % 3 != 0 || column_of_cone[k] >= m_size)
        THROW_EXCEPTION(
            "fem::nonsmooth_formulations::Plasticity: unexpected interaction position");
      ordered[static_cast<size_t>(column_of_cone[k] / 3)] = values[static_cast<size_t>(k)];
    }
    return ordered;
  };
  std::vector<double> sigma_y_c = by_column(sigma_y), eta_c = by_column(eta),
                      theta_c = by_column(theta), cohesion_c = by_column(cohesion);
  bool has_cohesion =
      std::any_of(cohesion_c.begin(), cohesion_c.end(), [](double c) { return c != 0.0; });
  Vec qp_copy = qp;
  Plasticity_VonMises_model von_mises{};
  Plasticity_DruckerPrager_model drucker_prager{};
  PlasticityProblem problem{};
  problem.dimension = 3;
  problem.numberOfCones = nc;
  problem.M = nm;
  problem.q = qp_copy.data();
  if (model == Model::von_mises) {
    von_mises.sigma_y = sigma_y_c.data();
    von_mises.hypothesis =
        (hypothesis == siconos::modeling::VonMisesPlasticityNSL::Hypothesis::plane_stress)
            ? PLASTICITY_VONMISES_PLANE_STRESS
            : PLASTICITY_VONMISES_PLANE_STRAIN;
    problem.model_type = PLASTICITY_MODEL_VON_MISES;
    problem.model.von_mises = &von_mises;
  } else {
    drucker_prager.eta = eta_c.data();
    drucker_prager.theta = theta_c.data();
    drucker_prager.cohesion = has_cohesion ? cohesion_c.data() : nullptr;
    problem.model_type = PLASTICITY_MODEL_DRUCKER_PRAGER;
    problem.model.drucker_prager = &drucker_prager;
  }
  problem.stress_coordinates = PLASTICITY_STRESS_COORDINATES_CARTESIAN_VOIGT;

  _n_compute++;
  if (!_dump_prefix.empty() &&
      (_dump_steps.empty() ||
       std::find(_dump_steps.begin(), _dump_steps.end(), _n_compute) != _dump_steps.end())) {
    std::string filename = _dump_prefix + "_" + std::to_string(_n_compute) + ".dat";
    FILE* file = fopen(filename.c_str(), "w");
    if (!file)
      THROW_EXCEPTION("fem::nonsmooth_formulations::Plasticity: cannot open " + filename);
    plasticity_printInFile(&problem, file);
    fclose(file);
  }

  // warm start with the previous y
  Vec y = (_w->size() == m_size) ? Vec(*_w) : Vec::Zero(m_size);
  Vec z = Vec::Zero(m_size);
  int info = plasticity_2d_driver(&problem, y.data(), z.data(), &*_numerics_solver_options);
  NM_free(nm);

  // --- back to the global unknowns
  Vec xc = SpMat(Hc_inv.transpose()) * (y - b);
  Vec xr = diagonal_rr ? Vec(Wrr_inv_qr - Wrr_inv_Wrc * xc)
                       : (Wrr.rows() > 0 ? Vec(lu_rr.solve(Vec(qr - Wrc * xc))) : Vec(0));
  Vec x = Pc.transpose() * xc + Pr.transpose() * xr;
  if (_globalVelocities->size() != n) _globalVelocities->resize(n);
  *_globalVelocities = x;
  _z->resize(m_size);
  *_z = z;
  _w->resize(m_size);
  *_w = y;

  DEBUG_EXPR({
    Vec res = W * x - q - H * z;
    std::cout << "Plasticity: |W x - q - H z| = " << res.norm() << std::endl;
  });
  return info;
}

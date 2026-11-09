#pragma once

// Sparse factorization of a square matrix (LDL^T when symmetric, LU
// otherwise), with a lazily grown compliance block: G = A^{-1} restricted
// to the dofs that have been asked for (for contact problems: the dofs
// that appear in H).  G grows one multi-rhs batch per new set of dofs, so
// once every contact dof has been seen, H A^{-1} H^T needs no linear solve.

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseLU>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace siconos::algebra {

class sparse_factor {
 public:
  using matrix_t = Eigen::SparseMatrix<double>;
  using index_t = matrix_t::StorageIndex;

  // a must be compressed
  explicit sparse_factor(const matrix_t& a) : _a(a) {
    if (_a.rows() != _a.cols()) {
      throw std::runtime_error("sparse_factor: matrix is not square");
    }
    _symmetric = is_symmetric(_a);
    if (_symmetric) {
      _ldlt.compute(_a);
      _use_ldlt = (_ldlt.info() == Eigen::Success);
    }
    if (!_use_ldlt) {
      _lu.compute(_a);
      if (_lu.info() != Eigen::Success) {
        throw std::runtime_error("sparse_factor: factorization failed");
      }
    }
  }

  const matrix_t& matrix() const { return _a; }

  // the compliance block is only consistent for a symmetric matrix:
  // G(j, m) = (A^-1 e_m)_j = (A^-1 e_j)_m
  bool symmetric() const { return _use_ldlt; }

  template <typename B>
  Eigen::MatrixXd solve(const B& b) const {
    if (_use_ldlt) {
      return _ldlt.solve(b);
    } else {
      return _lu.solve(b);
    }
  }

  // G is dense over the tracked dofs: past this many, callers solve instead.
  static constexpr std::size_t max_compliance_dofs = 4096;

  bool can_track(const std::vector<index_t>& dofs) const {
    std::size_t missing = 0;
    for (auto j : dofs) missing += (_slot.find(j) == _slot.end());
    return _slot_dof.size() + missing <= max_compliance_dofs;
  }

  // Make sure G holds the dofs in `dofs` (local indices).  Missing dofs are
  // solved for in one multi-rhs batch.
  void ensure_compliance(const std::vector<index_t>& dofs) {
    std::vector<index_t> missing;
    for (auto j : dofs) {
      if (_slot.find(j) == _slot.end()) {
        _slot.emplace(j, static_cast<index_t>(_slot_dof.size() + missing.size()));
        missing.push_back(j);
      }
    }
    if (missing.empty()) return;

    auto old_size = static_cast<Eigen::Index>(_slot_dof.size());
    auto new_size = old_size + static_cast<Eigen::Index>(missing.size());
    _slot_dof.insert(_slot_dof.end(), missing.begin(), missing.end());

    Eigen::MatrixXd rhs = Eigen::MatrixXd::Zero(_a.rows(), missing.size());
    for (std::size_t c = 0; c < missing.size(); ++c) {
      rhs(missing[c], c) = 1.0;
    }
    Eigen::MatrixXd x = solve(rhs);

    _g.conservativeResize(new_size, new_size);
    // new columns: every tracked dof against each new dof
    for (std::size_t c = 0; c < missing.size(); ++c) {
      auto col = old_size + static_cast<Eigen::Index>(c);
      for (Eigen::Index s = 0; s < new_size; ++s) {
        _g(s, col) = x(_slot_dof[s], c);
      }
    }
    // new rows for old columns, by symmetry
    for (Eigen::Index s = old_size; s < new_size; ++s) {
      for (Eigen::Index col = 0; col < old_size; ++col) {
        _g(s, col) = _g(col, s);
      }
    }
  }

  index_t slot(index_t dof) const { return _slot.at(dof); }

  // (A^-1)_{j m} for tracked dofs j, m
  double compliance(index_t slot_j, index_t slot_m) const { return _g(slot_j, slot_m); }

 private:
  static bool is_symmetric(const matrix_t& a) {
    matrix_t at = a.transpose();
    double norm = a.norm();
    double diff = (a - at).norm();
    return diff <= 1e-12 * norm;
  }

  matrix_t _a;
  bool _symmetric = false;
  bool _use_ldlt = false;
  Eigen::SimplicialLDLT<matrix_t> _ldlt;
  Eigen::SparseLU<matrix_t> _lu;

  std::unordered_map<index_t, index_t> _slot;  // dof -> slot in _g
  std::vector<index_t> _slot_dof;              // slot -> dof
  Eigen::MatrixXd _g;
};

// Factorizations of a FEM system for the Moreau-Jean scheme, kept for the
// whole simulation: W = M + h^2 theta^2 K (free state, W matrix, input)
// and M (M^{-1} (f - K q)).  W is rebuilt when h or theta change.  M and K
// are assumed constant (linear elasticity): after changing them, the
// operator must be invalidated (one_step_integrator invalidate_fem_operators).
struct fem_operator {
  double h = 0.;
  double theta = 0.;
  std::shared_ptr<sparse_factor> w;
  std::shared_ptr<sparse_factor> m;
};

// A diagonal block of an assembled block diagonal matrix whose
// factorization is already known.
struct factored_block {
  std::ptrdiff_t offset;
  std::ptrdiff_t size;
  sparse_factor* factor;
};

}  // namespace siconos::algebra

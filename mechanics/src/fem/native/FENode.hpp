// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file FiniteElementLinearTIDS.hpp

 */
#ifndef FENODEH
#define FENODEH

#include <cstddef>  //size_t
#include <memory>
#include <span>
#include <vector>

#include "SiconosMatrix.hpp"

namespace siconos::mechanics::fem {

class MeshVertex;

class FENode {
 private:
  /* node number */
  std::size_t num_ = 0;

  /* associated mesh vertex */
  std::shared_ptr<MeshVertex> mVertex_{nullptr};

  /* associated dof numbers in the global dof vector*/
  std::vector<siconos::algebra::Index> global_dof_index_ = {};

  /** Rule of five */
  FENode() = delete;
  FENode(FENode&) = delete;
  FENode& operator=(const FENode&) = delete;
  FENode(FENode&&) = delete;
  FENode& operator=(FENode&&) = delete;

 public:
  /** Constructor
      \param num node number
      \param v vertex associated to the node
      \param dofIndex global dof index
   */
  FENode(size_t num, std::shared_ptr<MeshVertex> v,
         const std::vector<siconos::algebra::Index>& dofIndex)
      : num_(num), mVertex_(v), global_dof_index_(dofIndex) {};

  ~FENode() noexcept = default;
  auto num() { return num_; }

  /** \returns associated dof numbers in the global dof vector*/
  std::span<const siconos::algebra::Index> global_dof_index() const {
    return global_dof_index_;
  };

  /** \returns x coordinate */
  double x() const;

  /** \returns y coordinate */
  double y() const;

  /** \returns z coordinate */
  double z() const;

  /** Print node parameters and information */
  void display();
};

}  // namespace siconos::mechanics::fem

#endif

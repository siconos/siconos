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

/*! \file FiniteElementLinearTIDS.hpp

 */
#ifndef FENODEH
#define FENODEH

#include <cstddef>  //size_t
#include <memory>
#include <span>
#include <vector>

namespace siconos::mechanics::fem {

class MeshVertex;

class FENode {
 private:
  /* node number */
  std::size_t num_ = 0;

  /* associated mesh vertex */
  std::shared_ptr<MeshVertex> mVertex_{nullptr};

  /* associated dof numbers in the global dof vector*/
  std::vector<size_t> global_dof_index_ = {};

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
  FENode(size_t num, std::shared_ptr<MeshVertex> v, const std::vector<size_t>& dofIndex)
      : num_(num), mVertex_(v), global_dof_index_(dofIndex) {};

  ~FENode() noexcept = default;
  auto num() { return num_; }

  /** \returns associated dof numbers in the global dof vector*/
  std::span<const size_t> global_dof_index() const { return global_dof_index_; };

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

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
#include <vector>

namespace siconos::mechanics::fem {

class MVertex;

class FENode {
 private:
  /* node number */
  std::size_t _num = 0;

  /* associated Mvertex */
  std::shared_ptr<MVertex> _mVertex{nullptr};

  /* associated dof number  in the global dof vector*/
  std::shared_ptr<std::vector<size_t>> _dofIndex{nullptr};

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
  FENode(size_t num, std::shared_ptr<MVertex> v, std::shared_ptr<std::vector<size_t>> dofIndex)
      : _num(num), _mVertex(v), _dofIndex(dofIndex){};

  ~FENode() noexcept = default;
  auto num() { return _num; }

  /** \returns associated dof number in the global dof vector*/
  auto dofIndex() { return _dofIndex; };

  /** \returns x coordinate */
  double x();

  /** \returns y coordinate */
  double y();

  /** \returns z coordinate */
  double z();

  /** Print node parameters and information */
  void display();
};

}  // namespace siconos::mechanics::fem

#endif

/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#include "FENode.hpp"

#include <iostream>

#include "Mesh.hpp"

double siconos::mechanics::fem::FENode::x() const { return mVertex_->x(); }

double siconos::mechanics::fem::FENode::y() const { return mVertex_->y(); }

double siconos::mechanics::fem::FENode::z() const { return mVertex_->z(); }

void siconos::mechanics::fem::FENode::display() {
  std::cout << "     - Fe Node - number (element model/mesh): " << num_ << "/ "
            << mVertex_->num() << "               - ndof:" << global_dof_index_.size()
            << "               - dofIndex (first/last): " << global_dof_index_.front() << ":"
            << global_dof_index_.back() << "\n";
};

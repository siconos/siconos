/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * SPDX-FileCopyrightText: INRIA
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once
#include "FENode.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "FiniteElementModel.hpp"
#include "Material.hpp"
#include "MeshUtils.hpp"
#include "siconos/vkernel/model/model_head.hpp"
#include "siconos/vkernel/storage/data_holder.hpp"

namespace siconos::vkernel::model {

using finite_element_linear_tids =
    storage::data_holder<mechanics::fem::FiniteElementLinearTIDS>;

using material = storage::data_holder<mechanics::fem::Material>;

}  // namespace siconos::vkernel::model

#pragma once
#include "FENode.hpp"
#include "FiniteElementModel.hpp"
#include "FiniteElementLinearTIDS.hpp"
#include "Material.hpp"
#include "MeshUtils.hpp"

#include "siconos/model/model_head.hpp"
#include "siconos/storage/data_holder.hpp"

namespace siconos::model {

using finite_element_linear_tids =
    storage::data_holder<mechanics::fem::FiniteElementLinearTIDS>;

using material = storage::data_holder<mechanics::fem::Material>;

}  // namespace siconos::model

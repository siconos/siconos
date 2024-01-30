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
#include "Equality.hpp"

#include "NonSmoothLaw.hpp"
#include "NumericsSolversNamespace.h"  // solver_options stuff
#include "NumericsToolsNamespace.h"    // .
#include "OSNSMatrix.hpp"
#include "SiconosVector.hpp"
#include "Simulation.hpp"
#include "TypeName.hpp"  // check nslaw type, should be replaced by dynamic_cast or variant ?

siconos::nonsmooth_formulations::Equality::Equality(int numericsSolverId) : LinearOSNS()
//:
// Equality(std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
//                            solver_options_delete))
{
}

siconos::nonsmooth_formulations::Equality::Equality(
    std::shared_ptr<SolverOptions> options)
    : LinearOSNS(options)
{
}

bool siconos::nonsmooth_formulations::Equality::checkCompatibleNSLaw(
    siconos::modeling::NonSmoothLaw& nslaw)
{
  float type_number = (float)(siconos::types::type_value(nslaw));
  _nslawtype.insert(type_number);

  if (not(siconos::types::type_value(nslaw) ==
          siconos::modeling::Type::EqualityConditionNSL)) {
    THROW_EXCEPTION(
        "\nsiconos::nonsmooth_formulations::Equality::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with Equality one step nonsmooth problem. \n \
                      Compatible NonSmoothLaw are: EqualityConditionNSL\n");
    return false;
  }

  return true;
}
int siconos::nonsmooth_formulations::Equality::compute(double time)
{
  int info = 0;
  // --- Prepare data for EQUALITY computing ---
  bool cont = preCompute(time);
  if (!cont) return info;

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0) {
    auto* q_ = q()->data();
    auto* z_ = _z->data();
    for (decltype(_sizeOutput) i = 0; i < _sizeOutput; ++i) z_[i] = -q_[i];
    // info = NM_gesv(&*_M->numericsMatrix(), z_, true);
    // info =
    // NM_LU_solve(NM_preserve(&*_M->numericsMatrix()),
    // z_, 1);
    info = NM_LU_solve(&*_M->numericsMatrix(), z_, 1);

    // --- Recovering of the desired variables from EQUALITY output ---
    postCompute();
  }

  return info;
}

void siconos::nonsmooth_formulations::Equality::initialize(
    std::shared_ptr<siconos::simulation::Simulation> sim)
{
  // General initialize for LinearOSNS
  LinearOSNS::initialize(sim);
  // auto indexSet = simulation()->indexSet(levelMin());
  //_M = std::make_shared<OSNSMatrix>(indexSet,_numericsMatrixStorageType));
}

void siconos::nonsmooth_formulations::Equality::updateM()
{
  assert(0);
  // Get index set from Simulation
  auto& indexSet = *simulation()->indexSet(indexSetLevel());

  if (!_M) {
    // Creates and fills M using Interactionof indexSet
    _M = std::make_shared<OSNSMatrix>(indexSet, _numericsMatrixStorageType);
  }
  else {
    _M->setStorageType(_numericsMatrixStorageType);
    _M->fillM(indexSet);
  }
  _sizeOutput = _M->size();
}

void siconos::nonsmooth_formulations::Equality::display() const
{
  std::cout << "======= EQUALITY of size " << _sizeOutput << " with: \n";
  LinearOSNS::display();
}

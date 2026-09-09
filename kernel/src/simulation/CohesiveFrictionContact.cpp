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

/*! \file CohesiveFrictionContact.cpp
 *  \brief Implementation of CohesiveFrictionContact solver
 *
 *  This file implements the cohesive friction-contact solver using a coupled
 *  velocity-displacement formulation. The solver extends standard friction-contact
 *  to handle cohesive zone models (CZM) through a block matrix structure:
 *
 *  \f[
 *  \begin{bmatrix} W & V \\ U & X \end{bmatrix}
 *  \begin{bmatrix} r_v \\ r_u \end{bmatrix}
 *  +
 *  \begin{bmatrix} q_v \\ q_u \end{bmatrix}
 *  =
 *  \begin{bmatrix} u_v \\ u_u \end{bmatrix}
 *  \f]
 *
 *  See CohesiveFrictionContact.hpp for detailed documentation.
 */

#include "CohesiveFrictionContact.hpp"

#include <algorithm>

#include "FrictionContact.hpp"
// #include <boost/smart_ptr/shared_ptr.hpp>
// #include <memory>

#include "CohesiveZoneModelNIFNSL.hpp"
// #include "FrictionContactProblem.h"
#include "CohesiveFrictionContactProblem.h"
#include "Interaction.hpp"
#include "LinearOSNS.hpp"
#include "MoreauJeanOSI.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothDrivers.h"
#include "NonSmoothDynamicalSystem.hpp"
#include "NumericsMatrix.h"
#include "OSNSMatrix.hpp"
#include "Question.hpp"
#include "Simulation.hpp"

// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "SolverOptions.h"  // for solver_options_create, solver_options_delete
#include "siconos_debug.h"
// All laws compliant with FrictionContact must have a visitor function to return mu.
namespace siconos::nonsmooth_formulations::cohesive_friction_contact {

struct ForMu : public siconos::modeling::nonsmooth_laws::Question<double> {
  using Visitor::visit;

  void visit(const siconos::modeling::CohesiveZoneModelNIFNSL& nsl) override {
    DEBUG_EXPR(std::cout << "CohesiveZoneModelNIFNSL type " << std::endl;)
    answer = nsl.mu();
  }
  void visit(const siconos::modeling::NewtonImpactFrictionNSL& nsl) override {
    DEBUG_EXPR(std::cout << "NewtonImpactFrictionNSL type " << std::endl;);
    answer = nsl.mu();
  }
};

// struct ForC_n : public siconos::modeling::nonsmooth_laws::Question<double> {
//   using Visitor::visit;
//   siconos::modeling::Interaction& _inter;

//   ForC_n(siconos::modeling::Interaction& inter):_inter(inter) {};

//   void visit(const siconos::modeling::CohesiveZoneModelNIFNSL& nsl) override {
//     answer = nsl.r_cohesion(_inter)[0];
//   }
// };

// struct ForC_t : public siconos::modeling::nonsmooth_laws::Question<double> {
//   using Visitor::visit;
//   siconos::modeling::Interaction& _inter;

//   ForC_t(siconos::modeling::Interaction& inter):_inter(inter) {};

//   void visit(const siconos::modeling::CohesiveZoneModelNIFNSL& nsl) override {
//     answer = nsl.r_cohesion(_inter)[1];
//   }
// };

}  // namespace siconos::nonsmooth_formulations::cohesive_friction_contact

namespace siconos::nonsmooth_formulations {

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb, int numericsSolverId)
    : CohesiveFrictionContact(
          dimPb, std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
                                                solver_options_delete)) {}

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb,
                                                 std::shared_ptr<SolverOptions> options)
    : FrictionContact(dimPb, options),  _scaling_as_percussion(true) {
  _assemblyType = LinearOSNSAssemblyType::REDUCED_DIRECT;
  _numericsMatrixStorageType = NM_SPARSE;
  if (dimPb == 3) {
    _cohesiveFrictionContact_driver = &cohesive_friction_3d_driver;
  } else
    THROW_EXCEPTION(
        "Wrong dimension value (must be 3) for CohesiveFrictionContact constructor.");

  _c_n = std::make_shared<std::vector<double>>();
  _c_t = std::make_shared<std::vector<double>>();
}

void CohesiveFrictionContact::initialize(
    std::shared_ptr<siconos::simulation::Simulation> simulation) {
  DEBUG_BEGIN("CohesiveFrictionContact::initialize()\n");

  FrictionContact::initialize(simulation);

  // Initialize cohesion vector
  if (!_q_cohesion) {
    _q_cohesion = std::make_shared<siconos::algebra::SiconosVector>(LinearOSNS::maxSize());
    _q_cohesion->setZero();
  }

  // Reserve vectors for cohesion intensities friction coefficients
  auto size = simulation->nonSmoothDynamicalSystem()->topology()->indexSet(0)->size();
  _c_n->reserve(size);
  _c_t->reserve(size);

  // Initialize H0, V, U, X, matrices for cohesive contribution
  // Note: V matrix size will be determined by the number of cohesive interactions in indexset0

  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT) {
    switch (_numericsMatrixStorageType) {
      case NM_SPARSE: {
        if (!_H0) {
          _H0 = std::make_shared<OSNSMatrix>(
              simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
              simulation->indexSet(_indexSetLevel)->size(), _numericsMatrixStorageType);
        }
        if (!_V) {
          _V = std::make_shared<OSNSMatrix>(0, 0, _numericsMatrixStorageType);
        }
        if (!_U) {
          _U = std::make_shared<OSNSMatrix>(0, 0, _numericsMatrixStorageType);
        }
        if (!_X) {
          _X = std::make_shared<OSNSMatrix>(0, 0, _numericsMatrixStorageType);
        }
        break;
      }
        {
          default:
            THROW_EXCEPTION("CohesiveFrictionContact::initialize unknown _storageType");
        }
    }
  } else {
  }

  DEBUG_END("CohesiveFrictionContact::initialize()\n");
}
void CohesiveFrictionContact::updateCoefficients() {
  _mu->clear();
  auto indexSet = simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;

  for (std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui) {
    // auto nsl = std::dynamic_pointer_cast<siconos::modeling::CohesiveZoneModelNIFNSL>(
    //     indexSet->bundle(*ui)->nonSmoothLaw());
    // assert(nsl);
    auto mu_val = siconos::modeling::nonsmooth_laws::ask<
        siconos::nonsmooth_formulations::cohesive_friction_contact::ForMu>(
        *indexSet->bundle(*ui)->nonSmoothLaw());

    _mu->push_back(mu_val);
  }

  _c_n->clear();
  _c_t->clear();

  auto indexSet0 = simulation()->indexSet(0);
  for (std::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui) {
    siconos::modeling::Interaction& inter = *indexSet0->bundle(*ui);
    auto nslaw = (std::dynamic_pointer_cast<siconos::modeling::CohesiveZoneModelNIFNSL>(
        indexSet->bundle(*ui)->nonSmoothLaw()));
    if (nslaw) {
      auto c_n_val = nslaw->cohesion(inter)[0];
      if (_scaling_as_percussion) {
        _c_n->push_back(c_n_val * simulation()->currentTimeStep());
      } else {
	_c_n->push_back(c_n_val);
      }
      auto c_t_val = nslaw->cohesion(inter)[1];
      if (_scaling_as_percussion) {
        _c_t->push_back(c_t_val * simulation()->currentTimeStep());
      } else {
	_c_t->push_back(c_t_val);
      }
    }

    // An attempt with visitor Question to avoid the dynamic cast
    // --> does not work for the moment
    // cohesive_friction_contact::ForC_n for_c_n = cohesive_friction_contact::ForC_n(inter);

    // auto c_n_val = siconos::modeling::nonsmooth_laws::ask<for_c_n>(
    //     *indexSet->bundle(*ui)->nonSmoothLaw());

    //_c_n->push_back(c_n_val);
    // auto c_t_val =
    // siconos::modeling::nonsmooth_laws::ask<cohesive_friction_contact::ForC_t>(
    //     *indexSet->bundle(*ui)->nonSmoothLaw());

    // _c_t->push_back(c_t_val);
  }
  DEBUG_EXPR(
      std::cout << "_c_n = ["; bool first = true; for (double x : *_c_n) {
        if (!first) std::cout << ", ";
        std::cout << x;
        first = false;
      } std::cout << "]\n";
      std::cout << "_c_t = ["; for (double x : *_c_t) {
        if (!first) std::cout << ", ";
        std::cout << x;
        first = false;
      } std::cout << "]\n";);
}
std::shared_ptr<CohesiveFrictionContactProblem>
CohesiveFrictionContact::cohesiveFrictionContactProblem() {
  auto numerics_problem = std::make_shared<CohesiveFrictionContactProblem>();
  numerics_problem->dimension = _contactProblemDim;
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;

  numerics_problem->numberOfCohesivePoints = _sizeOutput_cohesion / _contactProblemDim;
  numerics_problem->M = NULL;
  numerics_problem->W = &*_M->numericsMatrix();
  numerics_problem->U = &*_U->numericsMatrix();
  numerics_problem->V = &*_V->numericsMatrix();
  numerics_problem->X = &*_X->numericsMatrix();
  numerics_problem->q_v = &*_q->data();
  numerics_problem->q_u = &*_q_cohesion->data();
  numerics_problem->mu = _mu->data();
  numerics_problem->c_n = _c_n->data();
  numerics_problem->c_t = _c_t->data();

  return numerics_problem;
}

void CohesiveFrictionContact::compute_q_cohesion_block(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::algebra::Index pos) {
  DEBUG_BEGIN("CohesiveFrictionContact::compute_q_cohesion_block()\n");

  auto indexSet = simulation()->indexSet(0);

  auto& osi1 = *indexSet->properties(vertex_inter).osi1;
  auto& osi2 = *indexSet->properties(vertex_inter).osi2;

  auto osi1_type = osi1.getType();
  auto osi2_type = osi2.getType();

  auto inter = indexSet->bundle(vertex_inter);
  auto nslaw = inter->nonSmoothLaw();
  auto nslaw_size = nslaw->size();

  // Check if this is a cohesive law
  auto cohesive_nslaw =
      std::dynamic_pointer_cast<siconos::modeling::CohesiveZoneModelNIFNSL>(nslaw);

  if (cohesive_nslaw) {
    // Only supported for MoreauJeanOSI currently
    using siconos::integrators::IntegratorType;
    if ((osi1_type == IntegratorType::MOREAUJEANOSI &&
         osi2_type == IntegratorType::MOREAUJEANOSI)) {
      // Compute free output contribution from cohesion
      osi1.computeFreeOutputPosition(vertex_inter, this);

      // Get the work vector containing the cohesion contribution
      auto& work_vecs = indexSet->properties(vertex_inter).workVectors;
      // Note: OSNSP_RHS_COHESION index needs to be defined in MoreauJeanOSI
      // For now, use OSNSP_RHS as a placeholder
      if (work_vecs && (*work_vecs).size() > 0) {
        // Use first work vector as placeholder for cohesion
        auto& osnsp_rhs_position = *(*work_vecs)[tools::enum_to_index(
            siconos::integrators::MoreauJeanOSI::wk_inter::osnsp_rhs_position)];
        // Copy to _q_cohesion at position pos
        _q_cohesion->segment(pos, nslaw_size) = osnsp_rhs_position.segment(0, nslaw_size);
      }
    } else {
      THROW_EXCEPTION(
          "CohesiveFrictionContact::compute_q_cohesion_block not yet implemented for OSI "
          "types " +
          std::to_string(static_cast<int>(osi1_type)) + " and " +
          std::to_string(static_cast<int>(osi2_type)));
    }
  }

  DEBUG_EXPR(siconos::algebra::print(*_q_cohesion););
  DEBUG_END("CohesiveFrictionContact::computeQCohesionBlock()\n");
}

void CohesiveFrictionContact::compute_q_cohesion(double time) {
  DEBUG_BEGIN("CohesiveFrictionContact::updateQWithQCohesion()\n");

  auto indexSet0 = simulation()->indexSet(0);

  // Resize and zero _q_cohesion
  if (_q_cohesion->size() != _sizeOutput_cohesion) {
    _q_cohesion->resize(_sizeOutput_cohesion);
  }
  _q_cohesion->setZero();

  // Loop through all interactions in indexSet0
  siconos::algebra::Index pos = 0;
  for (auto [ui, uiend] = indexSet0->vertices(); ui != uiend; ++ui) {
    pos = indexSet0->properties(*ui).absolute_position;
    auto inter = indexSet0->bundle(*ui);
    compute_q_cohesion_block(*ui, pos);
  }

  DEBUG_EXPR(siconos::algebra::print(*_q_cohesion););

  DEBUG_END("CohesiveFrictionContact::updateQWithQCohesion()\n");
}

void CohesiveFrictionContact::computeMatrices() {
  DEBUG_BEGIN("CohesiveFrictionContact::computeMatrices()\n");
  // Compute matrix V that maps cohesive forces from indexSet0 to the OSNS problem
  // This is similar to the M matrix computation but for indexSet0 interactions

  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT) {
    siconos::graphs::InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
    siconos::graphs::InteractionsGraph& indexSet0 = *simulation()->indexSet(0);
    siconos::graphs::DynamicalSystemsGraph& DSG0 =
        *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

    // fill _Winverse
    // _W_inverse->fillWinverse(DSG0); already done in computeM

    // fill H (transpose)
    //_H->fillHtrans(DSG0, indexSet);  already done in computeM

    // fill H0
    _H0->fillH(DSG0, indexSet0);

    // DEBUG_EXPR(NumericsMatrix* H0_NM = &*(_H0->numericsMatrix()); std::cout << "H0 :";
    //            NM_display(H0_NM););

    // ComputeV
    _V->computeV(_H->numericsMatrix(), _W_inverse->numericsMatrix(), _H0->numericsMatrix());

    // ComputeU
    _U->computeU(_H->numericsMatrix(), _W_inverse->numericsMatrix(), _H0->numericsMatrix());

    // ComputeX
    _X->computeX(_H0->numericsMatrix(), _W_inverse->numericsMatrix());

  } else
    THROW_EXCEPTION("CohesiveFrictionContact::computeMatrices unknown _assemblyTYPE");

  DEBUG_EXPR(_V->display(););
  // NumericsMatrix* V_NM = _V->numericsMatrix().get();
  // printf("V_NM : %p\n", V_NM );
  // if (V_NM)
  //   NM_display(V_NM);

  // getchar();
  DEBUG_END("CohesiveFrictionContact::computeMatrices()\n");
}

bool CohesiveFrictionContact::preCompute(double time) {
  DEBUG_BEGIN("CohesiveFrictionContact::preCompute()\n");

  DEBUG_EXPR(
      std::cout << "indexSet0 size : " << simulation()->indexSet(0)->size() << std::endl;
      std::cout << "indexSet1 size : " << simulation()->indexSet(1)->size() << std::endl;);

  // First do standard preCompute
  // _M and _q are computed on indexSet 1
  bool hasContactActive = LinearOSNS::preCompute(time);

  // In the case, that indexSet1 is empty, we compute M to fill en empty marix !!
  if (!hasContactActive) {
    LinearOSNS::computeM();
  }

  // Compute coupling matrices between cohesive points and contact points.
  computeMatrices();

  _sizeOutput = _M->cols();
  _sizeOutput_cohesion = _X->cols();

  // Add cohesive contribution to q
  compute_q_cohesion(time);

  // rescaling

  // Ugly Hack to get theta
  // we should consider a vector of theta

  siconos::graphs::InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  siconos::graphs::InteractionsGraph& indexSet0 = *simulation()->indexSet(0);
  double theta = 0.0;
  for (auto [ui, uiend] = indexSet0.vertices(); ui != uiend; ++ui) {
    auto inter = indexSet0.bundle(*ui);
    auto osi1 = indexSet0.properties(*ui).osi1;
    auto osi1_type = osi1->getType();
    auto& osi2 = *indexSet0.properties(*ui).osi1;
    auto osi2_type = osi2.getType();
    using siconos::integrators::IntegratorType;
    if ((osi1_type == IntegratorType::MOREAUJEANOSI &&
         osi2_type == IntegratorType::MOREAUJEANOSI)) {
      auto moreaujean_osi1 =
          std::dynamic_pointer_cast<siconos::integrators::MoreauJeanOSI>(osi1);
      theta = moreaujean_osi1->theta();
    }
    break;
  }
  if (_scaling_as_percussion) {
    NM_scal(theta, &*(_U->numericsMatrix()));
    NM_scal(theta, &*(_X->numericsMatrix()));
    *_q_cohesion = *_q_cohesion / simulation()->timeStep();

    // for (double val : *_c_n) {
    //   val = val * simulation()->timeStep();
    //   }
    // for (double val : *_c_t) {
    //   val = val * simulation()->timeStep();
    //   }

    // *_c_n->data() *= simulation()->timeStep();
    // *_c_t->data() *= simulation()->timeStep();

  } else {
    NM_scal(theta * simulation()->timeStep(), &*(_U->numericsMatrix()));
    NM_scal(simulation()->timeStep(), &*(_V->numericsMatrix()));
    NM_scal(theta * simulation()->timeStep() * simulation()->timeStep(),
            &*(_X->numericsMatrix()));
  }

  if (_z->size() != (_sizeOutput + _sizeOutput_cohesion)) {
    _z->resize(_sizeOutput + _sizeOutput_cohesion, Eigen::NoChange);
    _z->setZero();
  }

  if (_w->size() != (_sizeOutput + _sizeOutput_cohesion)) {
    _w->resize(_sizeOutput + _sizeOutput_cohesion);
    _w->setZero();
  }

  DEBUG_END("CohesiveFrictionContact::preCompute()\n");
  return true;
}

void CohesiveFrictionContact::postCompute() {
  DEBUG_BEGIN("CohesiveFrictionContact::postCompute()\n");

  // // DEBUG_EXPR(
  // std::cout << "w: " ;siconos::algebra::print(*_w);
  //            // );
  // // DEBUG_EXPR(
  // std::cout << "z: " ;   siconos::algebra::print(*_z);// );

  // Call parent postCompute
  FrictionContact::postCompute();

  // std::cout <<  "indexSetLevel() :" << indexSetLevel() <<std::endl;
  // std::cout <<  "inputOutputLevel() :" << inputOutputLevel() <<std::endl;

  // We store the cohesive forces in lambda[0] for all the cohesive points in
  // indexSet0
  // Warning, when mixing law with FrictionContact for instance.

  auto& indexSet0 = *simulation()->indexSet(0);
  siconos::graphs::InteractionsGraph::VIterator ui, uiend;
  for (std::tie(ui, uiend) = indexSet0.vertices(); ui != uiend; ++ui) {
    auto& inter = *indexSet0.bundle(*ui);
    // Get the  position of inter-interactionBlock in the vector w or z
    auto pos = indexSet0.properties(*ui).absolute_position;
    // Get lambda for the current Interaction
    auto lambda = inter.lambda(0);
    // Copy _z values, starting from index pos + _sizeOutput_cohesion into lambda[0].
    if (_scaling_as_percussion) {
      lambda->segment(0, lambda->size()) =
          _z->segment(pos + _sizeOutput, lambda->size()) / simulation()->currentTimeStep();
    } else {
      lambda->segment(0, lambda->size()) = _z->segment(pos + _sizeOutput, lambda->size());
    }

    // auto lambda_1 = inter.lambda(1);

    // std::cout << "contact percussion   : ";
    // siconos::algebra::print(_z->segment(pos, lambda->size()));
    // std::cout << "cohesion percussion  : ";
    // siconos::algebra::print(_z->segment(pos + _sizeOutput, lambda->size()));
    // std::cout << "contact velocity     : ";
    // siconos::algebra::print(_w->segment(pos, lambda->size()));
    // std::cout << "cohesion displacement: ";
    // siconos::algebra::print(_w->segment(pos + _sizeOutput, lambda->size()));

    DEBUG_EXPR(siconos::algebra::print(*lambda););
  }

  DEBUG_END("CohesiveFrictionContact::postCompute()\n");
}

int siconos::nonsmooth_formulations::CohesiveFrictionContact::solve()
// std::shared_ptr<FrictionContactProblem> problem) {
{
  //  if (!problem) {
  auto problem = cohesiveFrictionContactProblem();
  //}
  //cohesiveFrictionContact_display(&*problem);
  cohesiveFrictionContactProblem_build_M_q_from_blocks(&*problem);
  // getchar();
  return (*_cohesiveFrictionContact_driver)(&*problem, &*_z->data(), &*_w->data(),
                                            &*_numerics_solver_options);
}

bool CohesiveFrictionContact::checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) {
  // Accept both standard friction laws and cohesive laws
  auto cohesive_nslaw = dynamic_cast<siconos::modeling::CohesiveZoneModelNIFNSL*>(&nslaw);
  if (cohesive_nslaw) {
    return true;
  }
  // // Also accept standard NewtonImpactFrictionNSL
  // auto friction_nslaw = dynamic_cast<siconos::modeling::NewtonImpactFrictionNSL*>(&nslaw);
  // if (friction_nslaw) {
  //   return true;
  // }
  return false;
}

int siconos::nonsmooth_formulations::CohesiveFrictionContact::compute(double time) {
  int info = 0;
  // --- Prepare data for FrictionContact computing ---
  bool cont = preCompute(time);
  if (!cont) {
    return info;
  }
  // nothing to do
  if (indexSetLevel() == siconos::internal::LEVELMAX) {
    return info;
  }

  updateCoefficients();

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  DEBUG_EXPR(std::cout << "size_output :" << _sizeOutput
                       << " _sizeOutput_cohesion : " << _sizeOutput_cohesion << std::endl;);
  if (_sizeOutput + _sizeOutput_cohesion != 0) {
    // Call Numerics Driver for FrictionContact
    info = solve();
    postCompute();
  }
  // display();
  // getchar();

  return info;
}
void CohesiveFrictionContact::display() const {
  std::cout << "======= CohesiveFrictionContact display =======\n";
  FrictionContact::display();
  std::cout << "Cohesion contribution size: " << _sizeOutput_cohesion << "\n";
  std::cout << "_X  ";
  if (_X)
    _X->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "_U  ";
  if (_U)
    _U->display();
  else
    std::cout << "-> nullptr" << std::endl;
  std::cout << "_V  ";
  if (_V)
    _V->display();
  else
    std::cout << "-> nullptr" << std::endl;

  if (_q_cohesion) {
    std::cout << "q_cohesion:\n";
    siconos::algebra::print(*_q_cohesion);
  }
  std::cout << std::endl;
  std::cout << "The CohesiveFrictionContact works on the index set of level  "
            << _indexSetLevel << " for contacts points and 0 for cohesive points" << std::endl;

  std::cout << "================================================\n";
}

}  // namespace siconos::nonsmooth_formulations

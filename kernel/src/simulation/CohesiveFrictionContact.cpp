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

#include "CohesiveFrictionContact.hpp"

#include "CohesiveZoneModelNIFNSL.hpp"
#include "FrictionContactProblem.h"
#include "MoreauJeanOSI.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OSNSMatrix.hpp"
#include "Simulation.hpp"
#include "Topology.hpp"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"
#include "SolverOptions.h"  // for solver_options_create, solver_options_delete

namespace siconos::nonsmooth_formulations {

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb, int numericsSolverId)
    : CohesiveFrictionContact(
          dimPb, std::shared_ptr<SolverOptions>(solver_options_create(numericsSolverId),
                                                solver_options_delete)) {}

CohesiveFrictionContact::CohesiveFrictionContact(int dimPb,
                                                  std::shared_ptr<SolverOptions> options)
    : FrictionContact(dimPb, options) {}

void CohesiveFrictionContact::initialize(
    std::shared_ptr<siconos::simulation::Simulation> simulation) {
  DEBUG_BEGIN("CohesiveFrictionContact::initialize()\n");

  FrictionContact::initialize(simulation);

  // Initialize cohesion vector
  if (!_q_cohesion) {
    _q_cohesion = std::make_shared<siconos::algebra::SiconosVector>(LinearOSNS::maxSize());
    _q_cohesion->setZero();
  }

  // Initialize V matrix for cohesive contribution
  // Note: V matrix size will be determined by the number of cohesive interactions

  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_BLOCK or
      _assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT) {
    if (!_V) {
      switch (_numericsMatrixStorageType) {
        case NM_DENSE:
        case NM_SPARSE: {
          _V = std::make_shared<OSNSMatrix>(0 , 0, _numericsMatrixStorageType);
          break;
        }
        case NM_SPARSE_BLOCK: {
          // = number of Interactionin the largest considered indexSet
          if (indexSetLevel() != siconos::internal::LEVELMAX &&
              simulation->nonSmoothDynamicalSystem()->topology()->indexSetsSize() >
                  indexSetLevel()) {
            _V = std::make_shared<OSNSMatrix>(0 ,simulation->indexSet(indexSetLevel())->size(),
                                              _numericsMatrixStorageType);
          } else {
            _V = std::make_shared<OSNSMatrix>(0 , 1, _numericsMatrixStorageType);
          }
          break;
        }
          {
            default:
              THROW_EXCEPTION("LinearOSNS::initOSNSMatrix unknown _storageType");
          }
      }
    }
  }

  // Initialize H0 for direct assembly if needed  
  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_DIRECT)
    {
      if(!_H0)
	{

	  switch(_numericsMatrixStorageType)
	    {
	    case NM_DENSE:
	      {
		_H0 = std::make_shared<OSNSMatrix>(LinearOSNS::maxSize(), LinearOSNS::maxSize(), NM_DENSE);
		break;
	      }
	    case NM_SPARSE:
	      {
              _H0 = std::make_shared<OSNSMatrix>(
                  simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
                  simulation->indexSet(_indexSetLevel)->size(), NM_SPARSE);
		break;
	      }
	    case NM_SPARSE_BLOCK:
	      {
              _H0 = std::make_shared<OSNSMatrix>(
                  simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(),
                  simulation->indexSet(_indexSetLevel)->size(), NM_SPARSE_BLOCK);
		break;
	      }
	      {
		default:
		  THROW_EXCEPTION("LinearOSNS::initOSNSMatrix unknown _storageType");
	      }
	    }
	}
    }  

  DEBUG_END("CohesiveFrictionContact::initialize()\n");
}

void CohesiveFrictionContact::computeQCohesionBlock(
    siconos::graphs::InteractionsGraph::VDescriptor& vertex_inter,
    siconos::algebra::Index pos) {
  DEBUG_BEGIN("CohesiveFrictionContact::computeQCohesionBlock()\n");

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
      osi1.computeFreeOutput(vertex_inter, this);

      // Get the work vector containing the cohesion contribution
      auto& work_vecs = indexSet->properties(vertex_inter).workVectors;
      // Note: OSNSP_RHS_COHESION index needs to be defined in MoreauJeanOSI
      // For now, use OSNSP_RHS as a placeholder
      if (work_vecs && (*work_vecs).size() > 0) {
        // Use first work vector as placeholder for cohesion
        auto& osnsp_rhs_cohesion = *(*work_vecs)[tools::enum_to_index(
            siconos::integrators::MoreauJeanOSI::wk_inter::osnsp_rhs_cohesion)];        
        // Copy to _q_cohesion at position pos
        _q_cohesion->segment(pos, nslaw_size) = osnsp_rhs_cohesion.segment(0, nslaw_size);
      }
    } else {
      THROW_EXCEPTION(
          "CohesiveFrictionContact::computeQCohesionBlock not yet implemented for OSI types " +
          std::to_string(static_cast<int>(osi1_type)) + " and " +
          std::to_string(static_cast<int>(osi2_type)));
    }
  }

  DEBUG_EXPR(siconos::algebra::print(*_q_cohesion););
  DEBUG_END("CohesiveFrictionContact::computeQCohesionBlock()\n");
}

void CohesiveFrictionContact::updateQWithQCohesion(double time) {
  DEBUG_BEGIN("CohesiveFrictionContact::updateQWithQCohesion()\n");

  auto indexSet0 = simulation()->indexSet(0);

  // Note: sizeColumn() might not exist, use alternative
  _sizeOutput_cohesion = _V->cols();

  DEBUG_PRINTF("_sizeOutput_cohesion = %d\n", _sizeOutput_cohesion);

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
    computeQCohesionBlock(*ui, pos);
  }

  DEBUG_EXPR(siconos::algebra::print(*_q_cohesion););

  // Add cohesion contribution to q: q = q + V * q_cohesion
  if (_V && _V->numericsMatrix()) {
    DEBUG_PRINT("Before adding cohesion:");
    DEBUG_EXPR(siconos::algebra::print(*_q););

    // Use Numerics matrix-vector product
    NM_gemv(1.0, _V->numericsMatrix().get(), _q_cohesion->data(), 1.0, _q->data());

    DEBUG_PRINT("After adding cohesion:");
    DEBUG_EXPR(siconos::algebra::print(*_q););
  }

  DEBUG_END("CohesiveFrictionContact::updateQWithQCohesion()\n");
}

void CohesiveFrictionContact::computeV() {
  DEBUG_BEGIN("CohesiveFrictionContact::computeV()\n");
  // Compute matrix V that maps cohesive forces from indexSet0 to the OSNS problem
  // This is similar to the M matrix computation but for indexSet0 interactions
  
  if (_assemblyType == LinearOSNSAssemblyType::REDUCED_BLOCK)
  {

    siconos::graphs::InteractionsGraph& indexSet0 = *simulation()->indexSet(0);
    siconos::graphs::InteractionsGraph& indexSet1 = *simulation()->indexSet(1);
    indexSet0.update_vertices_indices();
    indexSet0.update_edges_indices();
    // Computes new _interactionBlocks if required
    updateInteractionBlocks(indexSet0);

    _V->fillV(indexSet1, indexSet0, !_hasBeenUpdated);
    DEBUG_PRINT("partial V");
    DEBUG_EXPR( _V->display(););


  }
    else if (_assemblyType ==LinearOSNSAssemblyType::REDUCED_DIRECT)
  {
     siconos::graphs::InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
     siconos::graphs::InteractionsGraph& indexSet0 = *simulation()->indexSet(0);
     siconos::graphs::DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

    // fill _Winverse
    // _W_inverse->fillWinverse(DSG0); already done in computeM

    // fill H (transpose)
    //_H->fillHtrans(DSG0, indexSet);  already done in computeM

    // fill H0
    _H0->fillH(DSG0, indexSet0);

    // ComputeV
    _V->computeV(_H->numericsMatrix(), _W_inverse->numericsMatrix(), _H0->numericsMatrix());

  }
  else
    THROW_EXCEPTION("CohesiveFrictionContact::computeV unknown _assemblyTYPE");


  DEBUG_EXPR(_V->display(););
  // NumericsMatrix* V_NM = _V->numericsMatrix().get();
  // printf("V_NM : %p\n", V_NM );
  // if (V_NM)
  //   NM_display(V_NM);

  // getchar();
  DEBUG_END("CohesiveFrictionContact::computeV()\n");
}

bool CohesiveFrictionContact::preCompute(double time) {
  DEBUG_BEGIN("CohesiveFrictionContact::preCompute()\n");

  // First do standard preCompute
  // _M and _q are computed on indexSet 1
  bool hasContactActive = LinearOSNS::preCompute(time);

  if (!hasContactActive) return false;

  // Update cohesion contribution
  computeV();

  //Add cohesive contribution to q
  updateQWithQCohesion(time);

  siconos::graphs::InteractionsGraph& indexSet = *simulation()->indexSet(indexSetLevel());
  if(_keepLambdaAndYState)
  {
    siconos::graphs::InteractionsGraph::VIterator ui, uiend;
    for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      auto& inter = *indexSet.bundle(*ui);
      auto nslaw = inter.nonSmoothLaw();
      auto nslaw_CohesiveZoneModelNIFNSL(std::dynamic_pointer_cast<siconos::modeling::CohesiveZoneModelNIFNSL>(nslaw));
      if (nslaw_CohesiveZoneModelNIFNSL)
      {
        // Get the position of inter-interactionBlock in the vector w
        // or z
        unsigned int pos = indexSet.properties(*ui).absolute_position;
	//auto osnsp_rhs_cohesion = *(*indexSet.properties(*ui).workVectors)[siconos::integrators::MoreauJeanOSI::OSNSP_RHS_COHESION];
	auto& osnsp_rhs_cohesion = *(*indexSet.properties(*ui).workVectors)[tools::enum_to_index(siconos::integrators::MoreauJeanOSI::wk_inter::osnsp_rhs_cohesion)];

        for (int k =0; k < osnsp_rhs_cohesion.size(); k++)
        {
          (*_z)(pos+k) -= osnsp_rhs_cohesion(k);
        }
      }
    }
  }

  
  DEBUG_END("CohesiveFrictionContact::preCompute()\n");
  return true;
}

void CohesiveFrictionContact::postCompute() {
  DEBUG_BEGIN("CohesiveFrictionContact::postCompute()\n");

  DEBUG_EXPR(siconos::algebra::print(*_w););
  DEBUG_EXPR(siconos::algebra::print(*_z););  

  // Call parent postCompute
  FrictionContact::postCompute();
  std::cout <<  "indexSetLevel() :" << indexSetLevel() <<std::endl;
  std::cout <<  "inputOutputLevel() :" << inputOutputLevel() <<std::endl;
  
  DEBUG_END("CohesiveFrictionContact::postCompute()\n");
}

bool CohesiveFrictionContact::checkCompatibleNSLaw(siconos::modeling::NonSmoothLaw& nslaw) {
  // Accept both standard friction laws and cohesive laws
  auto cohesive_nslaw =
      dynamic_cast<siconos::modeling::CohesiveZoneModelNIFNSL*>(&nslaw);
  if (cohesive_nslaw) {
    return true;
  }
  // Also accept standard NewtonImpactFrictionNSL
  auto friction_nslaw =
      dynamic_cast<siconos::modeling::NewtonImpactFrictionNSL*>(&nslaw);
  if (friction_nslaw) {
    return true;
  }
  return false;
}

void CohesiveFrictionContact::display() const {
  std::cout << "======= CohesiveFrictionContact display =======\n";
  FrictionContact::display();
  std::cout << "Cohesion contribution size: " << _sizeOutput_cohesion << "\n";
  if (_q_cohesion) {
    std::cout << "q_cohesion:\n";
    siconos::algebra::print(*_q_cohesion);
  }
  std::cout << "================================================\n";
}

}  // namespace siconos::nonsmooth_formulations

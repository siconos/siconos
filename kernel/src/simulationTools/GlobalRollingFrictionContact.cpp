/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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
#include "SiconosPointers.hpp"
#include "NumericsMatrix.h"
#include "GlobalRollingFrictionContact.hpp"
#include "Simulation.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relation.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "MoreauJeanGOSI.hpp" // Numerics Header
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonImpactRollingFrictionNSL.hpp"
#include "OSNSMatrix.hpp"

#include "TypeName.hpp"

// --- Numerics headers ---
#include "NonSmoothDrivers.h"
#include "gfc3d_Solvers.h"
#include "NumericsSparseMatrix.h"
// #define DEBUG_NOCOLOR
// #define DEBUG_STDOUT
// #define DEBUG_MESSAGES
#include "siconos_debug.h"

// Constructor from solver id - Uses delegated constructor
GlobalRollingFrictionContact::GlobalRollingFrictionContact(int dimPb, const int numericsSolverId):
  GlobalRollingFrictionContact(dimPb, SP::SolverOptions(solver_options_create(numericsSolverId),
                        solver_options_delete))
{}

// Constructor based on a pre-defined solver options set.
GlobalRollingFrictionContact::GlobalRollingFrictionContact(int dimPb, SP::SolverOptions options):
  GlobalFrictionContact(dimPb,options), _g_rolling_driver(&g_rolling_fc3d_driver)
{
  // Only rolling fc3d for the moment.
  if(_contactProblemDim != 5)
    THROW_EXCEPTION("GlobalRollingFrictionContact No solver for 2 dimensional problems");

  //Reset default storage type for numerics matrices.
  _numericsMatrixStorageType = NM_SPARSE;
}


void GlobalRollingFrictionContact::initialize(SP::Simulation sim)
{

  GlobalFrictionContact::initialize(sim);
  // get topology
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Fill vector of rolling friction coefficients
  SP::InteractionsGraph I0 = topology->indexSet0();
  _mu_r.reset(new MuStorage());
  _mu_r->reserve(I0->size());

}

SP::GlobalRollingFrictionContactProblem GlobalRollingFrictionContact::globalRollingFrictionContactProblem()
{
  SP::GlobalRollingFrictionContactProblem numerics_problem(globalRollingFrictionContactProblem_new());
  numerics_problem->M = &*_W->numericsMatrix();
  numerics_problem->H = &*_H->numericsMatrix();
  numerics_problem->q = _q->getArray();
  numerics_problem->b = _b->getArray();
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->mu = _mu->data();
  numerics_problem->mu_r = _mu_r->data();
  numerics_problem->dimension = 5;
  return numerics_problem;
}

GlobalRollingFrictionContactProblem *GlobalRollingFrictionContact::globalRollingFrictionContactProblemPtr()
{
  GlobalRollingFrictionContactProblem *numerics_problem = &_numerics_problem;
  numerics_problem->M = &*_W->numericsMatrix();
  numerics_problem->H = &*_H->numericsMatrix();
  numerics_problem->q = _q->getArray();
  numerics_problem->b = _b->getArray();
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->mu = _mu->data();
  numerics_problem->mu_r = _mu_r->data();
  numerics_problem->dimension = 5;
  return numerics_problem;
}


bool GlobalRollingFrictionContact::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{

  float type_number= (float) (Type::value(nslaw) + 0.1 * nslaw.size());
  _nslawtype.insert(type_number);

  if (Type::value(nslaw) != Type::NewtonImpactRollingFrictionNSL)
  {
    THROW_EXCEPTION("\nGlobalRollingFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with GlobalRollingFrictionalContact one step nonsmooth problem. \n\
                      Compatible NonSmoothLaw is NewtonImpactRollingFrictionNSL (3D) \n");
    return false;
  }
  if (_nslawtype.size() > 1)
  {
    THROW_EXCEPTION("\nFrictionContact::checkCompatibleNSLaw -  \n\
                     Compatible NonSmoothLaw is : NewtonImpactRollingFrictionNSL (3D), but you cannot mix them \n");
    return false;
  }

  return true;
}



bool GlobalRollingFrictionContact::preCompute(double time)
{
  DEBUG_BEGIN("GlobalRollingFrictionContact::preCompute(double time)\n");
  // This function is used to prepare data for the GlobalRollingFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();
  DEBUG_PRINTF("indexSetLevel = %i\n", indexSetLevel());
  if(indexSetLevel() == LEVELMAX)
  {
    DEBUG_END("GlobalRollingFrictionContact::preCompute(double time)\n");
    return false;
  }
  if(!_hasBeenUpdated)
  {
    InteractionsGraph& indexSet = *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel);
    DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();



    // compute size and nnz of M and collect all matrices
    // compute nnz of H and collect H blocks
    // fill  mu

    // if (_sizeOutput == 0)
    // {
    //   DEBUG_END("GlobalRollingFrictionContact::preCompute(double time)\n");
    //   return false; }

    _mu->clear();
    _mu_r->clear();
//    _mu.reserve(indexSet.size())

#if !defined(SICONOS_USE_MAP_FOR_HASH)
    typedef std::unordered_map<SP::DynamicalSystem, SiconosMatrix*> dsMatMap;
    typedef std::unordered_map<SP::DynamicalSystem, size_t> dsPosMap;
#else
    typedef std::map<SP::DynamicalSystem, SiconosMatrix*> dsMatMap;
    typedef std::map<SP::DynamicalSystem, size_t> dsPosMap;
#endif
    dsMatMap dsMat;
    dsPosMap absPosDS;

    size_t sizeM = 0;


    // fill _W
    _W->fillW(DSG0);
    sizeM = _W->size();
    _sizeGlobalOutput = sizeM;
    DEBUG_PRINTF("sizeM = %lu \n", sizeM);


    // fill _q
    if(_q->size() != _sizeGlobalOutput)
      _q->resize(_sizeGlobalOutput);

    size_t offset = 0;
    DynamicalSystemsGraph::VIterator dsi, dsend;
    for(std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi)
    {
      SP::DynamicalSystem ds = DSG0.bundle(*dsi);
      Type::Siconos dsType = Type::value(*ds);
      size_t dss = ds->dimension();
      DEBUG_PRINTF("offset = %lu \n", offset);

      OneStepIntegrator& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
      OSI::TYPES osiType = Osi.getType();
      if(osiType == OSI::MOREAUJEANGOSI)
      {
        VectorOfVectors& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;

        if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
        {
          SiconosVector& vfree = *ds_work_vectors[MoreauJeanGOSI::FREE];
          setBlock(vfree, _q, dss, 0, offset);
        }
        else  if(dsType == Type::NewtonEulerDS)
        {
          SiconosVector& vfree = *ds_work_vectors[MoreauJeanGOSI::FREE];
          setBlock(vfree, _q, dss, 0, offset);
        }
      }
      else
      {
        THROW_EXCEPTION("GlobalRollingFrictionContact::computeq. Not yet implemented for Integrator type : " + std::to_string(osiType));
      }
      offset += dss;
    }
    DEBUG_EXPR(_q->display(););

    /************************************/


    // fill H
    _H->fillH(DSG0, indexSet);
    DEBUG_EXPR(NM_display(_H->numericsMatrix().get()););

    _sizeOutput =_H->sizeColumn();
    DEBUG_PRINTF("_sizeOutput = %i\n ", _sizeOutput);


    //fill _b
    if(_b->size() != _sizeOutput)
      _b->resize(_sizeOutput);

    size_t pos = 0;
    InteractionsGraph::VIterator ui, uiend;
    for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet.bundle(*ui);

      assert(Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactRollingFrictionNSL);
      _mu->push_back(std::static_pointer_cast<NewtonImpactRollingFrictionNSL>(inter->nonSmoothLaw())->mu());
      _mu_r->push_back(std::static_pointer_cast<NewtonImpactRollingFrictionNSL>(inter->nonSmoothLaw())->muR());

      SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
      SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;
      OneStepIntegrator& Osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
      OneStepIntegrator& Osi2 = *DSG0.properties(DSG0.descriptor(ds2)).osi;

      OSI::TYPES osi1Type = Osi1.getType();
      OSI::TYPES osi2Type = Osi2.getType();
      if(osi1Type == OSI::MOREAUJEANGOSI  && osi2Type == OSI::MOREAUJEANGOSI)
      {
        static_cast<MoreauJeanGOSI&>(Osi1).NSLcontrib(inter, *this);
      }
      else
      {
        THROW_EXCEPTION("GlobalRollingFrictionContact::computeq. Not yet implemented for Integrator type : " + std::to_string(osi1Type));
      }
      SiconosVector& osnsp_rhs = *(*indexSet.properties(*ui).workVectors)[MoreauJeanGOSI::OSNSP_RHS];
      pos =  indexSet.properties(*ui).absolute_position;
      size_t sizeY = inter->dimension();
      setBlock(osnsp_rhs, _b, sizeY, 0, pos);
    }
    DEBUG_EXPR(_b->display(););
    // Checks z and w sizes and reset if necessary
    if(_z->size() != _sizeOutput)
    {
      _z->resize(_sizeOutput, false);
      _z->zero();
    }

    if(_w->size() != _sizeOutput)
    {
      _w->resize(_sizeOutput);
      _w->zero();
    }

    if(_globalVelocities->size() != _sizeGlobalOutput)
    {
      _globalVelocities->resize(_sizeGlobalOutput);
      _globalVelocities->zero();
    }

  }
  DEBUG_END("GlobalRollingFrictionContact::preCompute(double time)\n");
  return true;
}

int GlobalRollingFrictionContact::compute(double time)
{
  DEBUG_BEGIN("GlobalRollingFrictionContact::compute(double time)\n")
  int info = 0;
  // --- Prepare data for GlobalRollingFrictionContact computing ---
  bool cont = preCompute(time);
  if(!cont)
    return info;
  updateMu();
  updateMur();
  // --- Call Numerics solver ---
  info= solve();
  DEBUG_EXPR(display(););
  postCompute();
  DEBUG_END("GlobalRollingFrictionContact::compute(double time)\n")
  return info;
}



int GlobalRollingFrictionContact::solve(SP::GlobalRollingFrictionContactProblem problem)
{
  if(!problem)
  {
    problem = globalRollingFrictionContactProblem();
  }
  return (*_g_rolling_driver)(&*problem,
                         _z->getArray(),
                         _w->getArray(),
                         _globalVelocities->getArray(),
                         &*_numerics_solver_options);
}




void GlobalRollingFrictionContact::updateMur()
{
  _mu_r->clear();
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    _mu_r->push_back(std::static_pointer_cast<NewtonImpactRollingFrictionNSL>
                     (indexSet->bundle(*ui)->nonSmoothLaw())->muR());
  }
}
void GlobalRollingFrictionContact::updateMu()
{
  _mu_r->clear();
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    _mu_r->push_back(std::static_pointer_cast<NewtonImpactRollingFrictionNSL>
                     (indexSet->bundle(*ui)->nonSmoothLaw())->mu());
  }
}
void GlobalRollingFrictionContact::display() const
{
  GlobalFrictionContact::display();


  std::cout << " - Vector mu : " <<std::endl;
  if(_mu)
  {
    for (int i = 0; i < _mu->size(); i++) {
      std::cout << (*_mu)[i] << " ";
    }
    std::cout << std::endl;
  }
  else std::cout << "-> nullptr" <<std::endl;

  std::cout << " - Vector mu_r : " <<std::endl;
  if(_mu_r)
  {
    for (int i = 0; i < _mu_r->size(); i++) {
      std::cout << (*_mu_r)[i] << " ";
    }
    std::cout  << std::endl;
  }
  else std::cout << "-> nullptr" <<std::endl;


}

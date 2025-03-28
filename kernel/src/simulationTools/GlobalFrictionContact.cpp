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
#include "SiconosPointers.hpp"
#include "NumericsMatrix.h"
#include "GlobalFrictionContact.hpp"
#include "Simulation.hpp"
//#include "Interaction.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relation.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "MoreauJeanOSI.hpp" // Numerics Header
#include "MoreauJeanGOSI.hpp" // Numerics Header
#include "LagrangianDS.hpp"
#include "NewtonEulerDS.hpp"
#include "NewtonImpactNSL.hpp"
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
GlobalFrictionContact::GlobalFrictionContact(int dimPb, const int numericsSolverId):
  GlobalFrictionContact(dimPb, SP::SolverOptions(solver_options_create(numericsSolverId),
                        solver_options_delete))
{}

// Constructor based on a pre-defined solver options set.
GlobalFrictionContact::GlobalFrictionContact(int dimPb, SP::SolverOptions options):
  LinearOSNS(options,GLOBAL), _contactProblemDim(dimPb), _gfc_driver(&gfc3d_driver)
{
  if(dimPb == 2)
  {
    _gfc_driver = &gfc2d_driver;
  }
  else if(dimPb == 3)
  {
    _gfc_driver = &gfc3d_driver;
  }

  //Reset default storage type for numerics matrices.
  _numericsMatrixStorageType = NM_SPARSE;
}


void GlobalFrictionContact::initVectorsMemory()
{
  // Memory allocation for reaction, and velocity
  LinearOSNS::initVectorsMemory();

  if(!_globalVelocities)
    _globalVelocities.reset(new SiconosVector(_maxSize));
  else
  {
    if(_globalVelocities->size() != _maxSize)
      _globalVelocities->resize(_maxSize);
  }

  if(!_b)
    _b.reset(new SiconosVector(_maxSize));
  else
  {
    if(_b->size() != _maxSize)
      _b->resize(_maxSize);
  }
}

void GlobalFrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  initVectorsMemory();

  // get topology
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  SP::InteractionsGraph I0 = topology->indexSet0();
  _mu.reset(new MuStorage());
  _mu->reserve(I0->size());


  initOSNSMatrix();


}

SP::GlobalFrictionContactProblem GlobalFrictionContact::globalFrictionContactProblem()
{
  SP::GlobalFrictionContactProblem numerics_problem(globalFrictionContactProblem_new());
  numerics_problem->M = &*_W->numericsMatrix();
  numerics_problem->H = &*_H->numericsMatrix();
  numerics_problem->q = _q->getArray();
  numerics_problem->b = _b->getArray();
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->mu = _mu->data();
  numerics_problem->dimension = _contactProblemDim;
  if (_assemblyType == GLOBAL_REDUCED)
  {
    numerics_problem->M_inverse = &*_W_inverse->numericsMatrix();
  }
  return numerics_problem;
}

GlobalFrictionContactProblem *GlobalFrictionContact::globalFrictionContactProblemPtr()
{
  GlobalFrictionContactProblem *numerics_problem = &_numerics_problem;
  numerics_problem->M = &*_W->numericsMatrix();
  numerics_problem->H = &*_H->numericsMatrix();
  numerics_problem->q = _q->getArray();
  numerics_problem->b = _b->getArray();
  numerics_problem->numberOfContacts = _sizeOutput / _contactProblemDim;
  numerics_problem->mu = _mu->data();
  numerics_problem->dimension = _contactProblemDim;
  if (_assemblyType == GLOBAL_REDUCED)
  {
    numerics_problem->M_inverse = &*_W_inverse->numericsMatrix();
  }
  return numerics_problem;
}


bool GlobalFrictionContact::checkCompatibleNSLaw(NonSmoothLaw& nslaw)
{

  float type_number= (float) (Type::value(nslaw) + 0.1 * nslaw.size());
  _nslawtype.insert(type_number);

  if (Type::value(nslaw) != Type::NewtonImpactFrictionNSL)
  {
    THROW_EXCEPTION("\nGlobalFrictionContact::checkCompatibleNSLaw -  \n\
                      The chosen nonsmooth law is not compatible with FrictionalContact one step nonsmooth problem. \n\
                      Compatible NonSmoothLaw is NewtonImpactFrictionNSL (2D or 3D) \n");
    return false;
  }
  if (_nslawtype.size() > 1)
  {
    THROW_EXCEPTION("\nFrictionContact::checkCompatibleNSLaw -  \n\
                     Compatible NonSmoothLaw is : NewtonImpactFrictionNSL (2D or 3D), but you cannot mix them \n");
    return false;
  }

  return true;
}

//#define WITH_TIMER


bool GlobalFrictionContact::preCompute(double time)
{
  DEBUG_BEGIN("GlobalFrictionContact::preCompute(double time)\n");
  // This function is used to prepare data for the GlobalFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();
  DEBUG_PRINTF("indexSetLevel = %i\n", indexSetLevel());
  if(indexSetLevel() == siconos::internal::LEVELMAX)
  {
    DEBUG_END("GlobalFrictionContact::preCompute(double time)\n");
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
    //   DEBUG_END("GlobalFrictionContact::preCompute(double time)\n");
    //   return false; }

    _mu->clear();
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

#ifdef WITH_TIMER
    std::chrono::time_point<std::chrono::system_clock> start, end, end_old;
    start = std::chrono::system_clock::now();
#endif
    // fill _W
    _W->fillW(DSG0);
    sizeM = _W->size();
    _sizeGlobalOutput = sizeM;
    DEBUG_PRINTF("sizeM = %lu \n", sizeM);
#ifdef WITH_TIMER
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::microseconds> (end-start).count();
    std::cout << "\nGlobalFrictionContact: fill W  " << elapsed << " ms" << std::endl;
#endif    

    
    if (_assemblyType == GLOBAL_REDUCED)
    {
      // fill _W_inverse
      _W_inverse->fillWinverse(DSG0);
    }
#ifdef WITH_TIMER
    end_old=end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>
      (end-end_old).count();
    std::cout << "GlobalFrictionContact: fillW inverse " << elapsed << " ms" << std::endl;
#endif
 
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

      //OneStepIntegrator& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
      //OSI::TYPES osiType = Osi.getType();
      SP::MoreauJeanGOSI mjgosi =  std::dynamic_pointer_cast<MoreauJeanGOSI>(DSG0.properties(DSG0.descriptor(ds)).osi);
      if(mjgosi)
      {
        VectorOfVectors& ds_work_vectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;
        if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::NewtonEulerDS)
        {
          SiconosVector& vfree = *ds_work_vectors[MoreauJeanGOSI::FREE];
          setBlock(vfree, _q, dss, 0, offset);
        }
      }
      else
      {
        THROW_EXCEPTION("GlobalFrictionContact::computeq. Not yet implemented for the given Integrator type : ");
      }
      offset += dss;
    }
    DEBUG_EXPR(_q->display(););
#ifdef WITH_TIMER
    end_old=end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>
      (end-end_old).count();
    std::cout << "GlobalFrictionContact: fill q " << elapsed << " ms" << std::endl;
#endif
 
    /************************************/


    // fill H
    _H->fillH(DSG0, indexSet);
    DEBUG_EXPR(NM_display(_H->numericsMatrix().get()););

    _sizeOutput =_H->sizeColumn();
    DEBUG_PRINTF("_sizeOutput = %i\n ", _sizeOutput);
#ifdef WITH_TIMER
    end_old=end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>
      (end-end_old).count();
    std::cout << "GlobalFrictionContact: fill H " << elapsed << " ms" << std::endl;
#endif

    //fill _b
    if(_b->size() != _sizeOutput)
      _b->resize(_sizeOutput);

    size_t pos = 0;
    InteractionsGraph::VIterator ui, uiend;
    for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet.bundle(*ui);

      assert(Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL);
      _mu->push_back(std::static_pointer_cast<NewtonImpactFrictionNSL>(inter->nonSmoothLaw())->mu()); //curious !!

      SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
      SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;
      OneStepIntegrator& osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
      OneStepIntegrator& osi2 = *DSG0.properties(DSG0.descriptor(ds2)).osi;

      if (typeid(osi1) == typeid(MoreauJeanGOSI) and typeid(osi2) == typeid(MoreauJeanGOSI))
      {
        //std::cout << "MoreauJeanGOSI case" << std::endl;
        SP::MoreauJeanGOSI mjgosi1 =  std::dynamic_pointer_cast<MoreauJeanGOSI>(DSG0.properties(DSG0.descriptor(ds1)).osi);
        mjgosi1->NonSmoothLawContributionToOutput(inter, *this);
      }
      else
      {
        THROW_EXCEPTION("GlobalFrictionContact::computeq. Not yet implemented for the given Integrator type ");
      }
      SiconosVector& osnsp_rhs = *(*indexSet.properties(*ui).workVectors)[MoreauJeanGOSI::OSNSP_RHS];
      pos =  indexSet.properties(*ui).absolute_position;
      size_t sizeY = inter->dimension();
      setBlock(osnsp_rhs, _b, sizeY, 0, pos);
    }
    DEBUG_EXPR(_b->display(););
#ifdef WITH_TIMER
    end_old=end;
    end = std::chrono::system_clock::now();
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>
      (end-end_old).count();
    std::cout << "GlobalFrictionContact: fill b " << elapsed << " ms" << std::endl;
#endif
    
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
  // nothing to do (IsLinear and not changed)
#ifdef WITH_TIMER
  end_old=end;
  end = std::chrono::system_clock::now();
  elapsed = std::chrono::duration_cast<std::chrono::microseconds>
    (end-end_old).count();
  std::cout << "GlobalFrictionContact: init w and z and V " << elapsed << " ms" << std::endl;
#endif
  }
  DEBUG_END("GlobalFrictionContact::preCompute(double time)\n");
  return true;
}

int GlobalFrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for GlobalFrictionContact computing ---
  bool cont = preCompute(time);
  if(!cont)
    return info;
  updateMu();
  
  // --- Call Numerics solver ---
  //if(_sizeGlobalOutput != 0)
  {
    info= solve();
    DEBUG_EXPR(display(););
    postCompute();
  }
  return info;
}



int GlobalFrictionContact::solve(SP::GlobalFrictionContactProblem problem)
{
  if(!problem)
  {
    problem = globalFrictionContactProblem();
  }
  return (*_gfc_driver)(&*problem,
                        _z->getArray(),
                        _w->getArray(),
                        _globalVelocities->getArray(),
                        &*_numerics_solver_options);
}



void GlobalFrictionContact::postCompute()
{
  DEBUG_BEGIN("GlobalFrictionContact::postCompute(double time)\n");

  // This function is used to set y/lambda values using output from primalfrictioncontact_driver
  // Only Interactions (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  InteractionsGraph& indexSet = *simulation()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel);
  // y and lambda vectors
  SP::SiconosVector  y, lambda;

  //   // === Loop through "active" Interactions (ie present in indexSets[1]) ===

  size_t pos = 0;

  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui, pos += _contactProblemDim)
  {
    Interaction& inter = *indexSet.bundle(*ui);
    // Get Y and Lambda for the current Interaction
    y = inter.y(inputOutputLevel());
    lambda = inter.lambda(inputOutputLevel());
    // Copy _w/_z values, starting from index pos into y/lambda.

    //setBlock(*_w, y, y->size(), pos, 0);// Warning: yEquivalent is
    // saved in y !!
    setBlock(*_z, lambda, lambda->size(), pos, 0);
    DEBUG_EXPR(lambda->display(););
  }
  DynamicalSystemsGraph& DSG0 = *simulation()->nonSmoothDynamicalSystem()->dynamicalSystems();

  unsigned int sizeDS;
  SP::OneStepIntegrator  Osi; 
  DynamicalSystemsGraph::VIterator dsi, dsend;
  pos=0;
  for(std::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi)
  {
    DynamicalSystem& ds = *DSG0.bundle(*dsi);
    Type::Siconos dsType = Type::value(ds);

    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&>(ds);
      sizeDS = d.dimension();
      SP::SiconosVector velocity = d.velocity();
      DEBUG_PRINTF("ds.number() : %i \n",ds.number());
      DEBUG_EXPR(velocity->display(););
      DEBUG_EXPR(_globalVelocities->display(););
      pos = DSG0.properties(*dsi).absolute_position;
      setBlock(*_globalVelocities, velocity, sizeDS, pos, 0);
      DEBUG_EXPR(velocity->display(););
    }
    else if(dsType == Type::NewtonEulerDS)
    {
      NewtonEulerDS& neds = static_cast<NewtonEulerDS&>(ds);;
      sizeDS = neds.dimension();
      SP::SiconosVector twist = neds.twist();
      DEBUG_PRINTF("ds.number() : %i \n",ds.number());
      DEBUG_EXPR(twist->display(););
      DEBUG_EXPR(_globalVelocities->display(););
      pos = DSG0.properties(*dsi).absolute_position;
      setBlock(*_globalVelocities, twist, sizeDS, pos, 0);
      DEBUG_EXPR(twist->display(););
    }
    else THROW_EXCEPTION("GlobalFrictionContact::postCompute() - not yet implemented for Dynamical system of type: " +  Type::name(ds));

  }

  DEBUG_END("GlobalFrictionContact::postCompute(double time)\n");

}

void GlobalFrictionContact::display() const
{

  std::cout << "===== " << _contactProblemDim << "D Global Friction Contact Problem " <<std::endl;
  std::cout << "size (_sizeOutput) " << _sizeOutput << "(ie " << _sizeOutput / _contactProblemDim << " contacts)."<<std::endl;
  std::cout << "and  size (_sizeGlobalOutput) " << _sizeGlobalOutput  <<std::endl;
  std::cout << "_numericsMatrixStorageType" << _numericsMatrixStorageType<< std::endl;
  std::cout << " - Matrix M  : " <<std::endl;
  // if (_W) _W->display();
  // else std::cout << "-> nullptr" <<std::endl;
  NumericsMatrix* W_NM = _W->numericsMatrix().get();
  if(W_NM)
  {
    NM_display(W_NM);
  }
  std::cout << " - Matrix H : " <<std::endl;
  // if (_H) _H->display();
  // else std::cout << "-> nullptr" <<std::endl;
  NumericsMatrix* H_NM = _H->numericsMatrix().get();
  if(H_NM)
  {
    NM_display(H_NM);
  }

  std::cout << " - Vector q : " <<std::endl;
  if(_q) _q->display();
  else std::cout << "-> nullptr" <<std::endl;
  std::cout << " - Vector b : " <<std::endl;
  if(_b) _b->display();
  else std::cout << "-> nullptr" <<std::endl;

  std::cout << " - Vector z (reaction) : " <<std::endl;
  if(_z) _z->display();
  else std::cout << "-> nullptr" <<std::endl;

  std::cout << " - Vector w (local velocities): " <<std::endl;
  if(_w) _w->display();
  else std::cout << "-> nullptr" <<std::endl;

  std::cout << " - Vector globalVelocities : " <<std::endl;
  if(_globalVelocities) _globalVelocities->display();
  else std::cout << "-> nullptr" <<std::endl;

  std::cout << "============================================================" <<std::endl;
}

void GlobalFrictionContact::updateMu()
{
  _mu->clear();
  SP::InteractionsGraph indexSet = simulation()->indexSet(indexSetLevel());
  InteractionsGraph::VIterator ui, uiend;
  for(std::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    _mu->push_back(std::static_pointer_cast<NewtonImpactFrictionNSL>
                   (indexSet->bundle(*ui)->nonSmoothLaw())->mu());
  }
}

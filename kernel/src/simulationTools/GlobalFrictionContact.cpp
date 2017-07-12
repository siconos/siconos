/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include "GlobalFrictionContact.hpp"
#include "Simulation.hpp"
//#include "Interaction.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relation.hpp"
#include "NewtonImpactFrictionNSL.hpp"
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
#include "debug.h"

// Constructor from a set of data
// Required input: simulation
// Optional: newNumericsSolverName
GlobalFrictionContact::GlobalFrictionContact(int dimPb, const int numericsSolverId):
  LinearOSNS(numericsSolverId), _contactProblemDim(dimPb)
{
  // Connect to the right function according to dim. of the problem
  if (_contactProblemDim == 2)
  {
    RuntimeException::selfThrow("GlobalFrictionContact No solver for 2 dimensional problems");
  }
  else if(_contactProblemDim == 3)
  {
    gfc3d_setDefaultSolverOptions(&*_numerics_solver_options, _numerics_solver_id);
    _gfc_driver = &gfc3d_driver;
  }
  else
  {
     RuntimeException::selfThrow("GlobalFrictionContact size not supported");
  }
  //default storage
  _numericsMatrixStorageType = NM_SPARSE;


}
GlobalFrictionContact::~GlobalFrictionContact()
{
  solver_options_delete(&*_numerics_solver_options);
}

void GlobalFrictionContact::initVectorsMemory()
{
  // Memory allocation for reaction, and velocity
  LinearOSNS::initVectorsMemory();

  if (!_globalVelocities)
    _globalVelocities.reset(new SiconosVector(_maxSize));
  else
  {
    if (_globalVelocities->size() != _maxSize)
      _globalVelocities->resize(_maxSize);
  }

  if (!_b)
    _b.reset(new SiconosVector(_maxSize));
  else
  {
    if (_b->size() != _maxSize)
      _b->resize(_maxSize);
  }
}


void GlobalFrictionContact::initOSNSMatrix()
{
  // Default size for M = _maxSize
  if (!_M)
  {
    // if (_numericsMatrixStorageType == NM_DENSE)
    //   _M.reset(new OSNSMatrix(_maxSize, NM_DENSE));
    // else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered graph of ds
    //   _M.reset(new OSNSMatrix(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), 1));

    switch (_numericsMatrixStorageType)
    {
    case NM_DENSE:
    {
      _M.reset(new OSNSMatrix(_maxSize, NM_DENSE));
      break;
    }
    case NM_SPARSE:
    {
      _M.reset(new OSNSMatrix(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), NM_SPARSE));
      break;
    }
    case NM_SPARSE_BLOCK:
    {
      _M.reset(new OSNSMatrix(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), NM_SPARSE_BLOCK));
      break;
    }
    {
    default:
      RuntimeException::selfThrow("GlobalFrictionContact::initOSNSMatrix unknown _storageType");
    }
    }
  }


  if (!_H)
  {

    switch (_numericsMatrixStorageType)
    {
    case NM_DENSE:
    {
      _H.reset(new OSNSMatrix(_maxSize, NM_DENSE));
      break;
    }
    case NM_SPARSE:
    {
      _H.reset(new OSNSMatrix(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), simulation()->indexSet(_indexSetLevel)->size()   , NM_SPARSE));
      break;
    }
    case NM_SPARSE_BLOCK:
    {
      _H.reset(new OSNSMatrix(simulation()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), simulation()->indexSet(_indexSetLevel)->size()   , NM_SPARSE_BLOCK));
      break;
    }
    {
    default:
      RuntimeException::selfThrow("GlobalFrictionContact::initOSNSMatrix unknown _storageType");
    }
    }
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
  DEBUG_PRINTF( "indexSetLevel = %i\n", indexSetLevel() );
  if (indexSetLevel() == LEVELMAX)
  {
    DEBUG_END("GlobalFrictionContact::preCompute(double time)\n");
    return false;
  }
  if (!_hasBeenUpdated)
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

#if defined(SICONOS_STD_UNORDERED_MAP) && !defined(SICONOS_USE_MAP_FOR_HASH)
    typedef std::unordered_map<SP::DynamicalSystem, SiconosMatrix*> dsMatMap;
    typedef std::unordered_map<SP::DynamicalSystem, size_t> dsPosMap;
#else
    typedef std::map<SP::DynamicalSystem, SiconosMatrix*> dsMatMap;
    typedef std::map<SP::DynamicalSystem, size_t> dsPosMap;
#endif
    dsMatMap dsMat;
    dsPosMap absPosDS;

    size_t sizeM = 0;


    // fill _M
    _M->fillM(DSG0);
    sizeM = _M->size();
    _sizeGlobalOutput = sizeM;
    DEBUG_PRINTF("sizeM = %lu \n", sizeM);


    // fill _q
    if (_q->size() != _sizeGlobalOutput)
      _q->resize(_sizeGlobalOutput);

    size_t offset = 0;
    DynamicalSystemsGraph::VIterator dsi, dsend;
    for(std11::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi)
    {
      SP::DynamicalSystem ds = DSG0.bundle(*dsi);
      Type::Siconos dsType = Type::value(*ds);
      size_t dss = ds->dimension();
      DEBUG_PRINTF("offset = %lu \n", offset);

      OneStepIntegrator& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
      OSI::TYPES osiType = Osi.getType();
      if (osiType == OSI::MOREAUJEANGOSI)
      {
        VectorOfVectors& workVectors = *DSG0.properties(DSG0.descriptor(ds)).workVectors;

        if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
        {
          SiconosVector& vfree = *workVectors[OneStepIntegrator::free];
          setBlock(vfree, _q, dss, 0, offset);
        }
        else  if (dsType == Type::NewtonEulerDS)
        {
          SiconosVector& vfree = *workVectors[OneStepIntegrator::free];
          setBlock(vfree, _q, dss, 0, offset);
        }
      }
      else
      {
        RuntimeException::selfThrow("GlobalFrictionContact::computeq. Not yet implemented for Integrator type : " + osiType);
      }
      offset += dss;
    }
    DEBUG_EXPR(_q->display(););

    /************************************/


     // fill H
    _H->fillH(DSG0, indexSet);
    DEBUG_EXPR(NM_display(_H->numericsMatrix().get() ););

    _sizeOutput =_H->sizeColumn();
    DEBUG_PRINTF( "_sizeOutput = %i\n ", _sizeOutput );


    //fill _b
    if (_b->size() != _sizeOutput)
      _b->resize(_sizeOutput);

    size_t pos = 0;
    InteractionsGraph::VIterator ui, uiend;
    for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      SP::Interaction inter = indexSet.bundle(*ui);

      assert(Type::value(*(inter->nonSmoothLaw())) == Type::NewtonImpactFrictionNSL);
      _mu->push_back(std11::static_pointer_cast<NewtonImpactFrictionNSL>(inter->nonSmoothLaw())->mu());

      SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
      SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;
      OneStepIntegrator& Osi1 = *DSG0.properties(DSG0.descriptor(ds1)).osi;
      OneStepIntegrator& Osi2 = *DSG0.properties(DSG0.descriptor(ds2)).osi;

      OSI::TYPES osi1Type = Osi1.getType();
      OSI::TYPES osi2Type = Osi2.getType();
      if (osi1Type == OSI::MOREAUJEANGOSI  && osi2Type == OSI::MOREAUJEANGOSI)
      {
        static_cast<MoreauJeanGOSI&>(Osi1).NSLcontrib(inter, *this);
      }
      else
      {
        RuntimeException::selfThrow("GlobalFrictionContact::computeq. Not yet implemented for Integrator type : " + osi1Type);
      }
      SiconosVector& osnsp_rhs = *(*indexSet.properties(*ui).workVectors)[MoreauJeanGOSI::OSNSP_RHS];
      pos =  indexSet.properties(*ui).absolute_position;
      setBlock(osnsp_rhs, _b, 3, 0, pos);
    }
    DEBUG_EXPR(_b->display(););
    // Checks z and w sizes and reset if necessary
    if (_z->size() != _sizeOutput)
    {
      _z->resize(_sizeOutput, false);
      _z->zero();
    }

    if (_w->size() != _sizeOutput)
    {
      _w->resize(_sizeOutput);
      _w->zero();
    }

    if (_globalVelocities->size() != _sizeGlobalOutput)
    {
      _globalVelocities->resize(_sizeGlobalOutput);
      _globalVelocities->zero();
    }

  }
  DEBUG_END("GlobalFrictionContact::preCompute(double time)\n");
  return true;
}

int GlobalFrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for GlobalFrictionContact computing ---
  bool cont = preCompute(time);
  if (!cont)
    return info;


  // --- Call Numerics solver ---
  //if (_sizeOutput != 0)
  {
    // The GlobalFrictionContact Problem in Numerics format
    GlobalFrictionContactProblem numerics_problem;
    globalFrictionContact_null(&numerics_problem);
    numerics_problem.M = &*_M->numericsMatrix();
    numerics_problem.H = &*_H->numericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.b = _b->getArray();
    numerics_problem.numberOfContacts = _sizeOutput / _contactProblemDim;
    numerics_problem.mu = _mu->data();
    numerics_problem.dimension = 3;
    DEBUG_EXPR(display(););
    info = (*_gfc_driver)(&numerics_problem,
                          _z->getArray(),
                           _w->getArray(),
                           _globalVelocities->getArray(),
			  &*_numerics_solver_options);
    DEBUG_EXPR(display(););
    postCompute();

  }

  return info;
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
  for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui, pos += 3)
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
  for(std11::tie(dsi, dsend) = DSG0.vertices(); dsi != dsend; ++dsi)
  {
    DynamicalSystem& ds = *DSG0.bundle(*dsi);
    Type::Siconos dsType = Type::value(ds);

    if(dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianLinearDiagonalDS)
    {
      LagrangianDS& d = static_cast<LagrangianDS&> (ds);
      sizeDS = d.dimension();
      SP::SiconosVector velocity = d.velocity();
      DEBUG_PRINTF("ds.number() : %i \n",ds.number());
      DEBUG_EXPR(velocity->display(););
      DEBUG_EXPR(_globalVelocities->display(););
      pos = DSG0.properties(*dsi).absolute_position;
      setBlock(*_globalVelocities, velocity, sizeDS, pos, 0 );
      DEBUG_EXPR(velocity->display(););
    }

    else RuntimeException::selfThrow("GlobalFrictionContact::postCompute() - not yet implemented for Dynamical system of type: " +  Type::name(ds));

  }

  DEBUG_END("GlobalFrictionContact::postCompute(double time)\n");

}

void GlobalFrictionContact::display() const
{

  std::cout << "===== " << _contactProblemDim << "D Primal Friction Contact Problem " <<std::endl;
  std::cout << "size (_sizeOutput) " << _sizeOutput << "(ie " << _sizeOutput / _contactProblemDim << " contacts)."<<std::endl;
  std::cout << "and  size (_sizeGlobalOutput) " << _sizeGlobalOutput  <<std::endl;
  std::cout << "_numericsMatrixStorageType" << _numericsMatrixStorageType<< std::endl;
  std::cout << " - Matrix M  : " <<std::endl;
  // if (_M) _M->display();
  // else std::cout << "-> NULL" <<std::endl;
  NumericsMatrix* M_NM = _M->numericsMatrix().get();
  if (M_NM)
  {
    NM_display(M_NM);
  }
  std::cout << " - Matrix H : " <<std::endl;
  // if (_H) _H->display();
  // else std::cout << "-> NULL" <<std::endl;
  NumericsMatrix* H_NM = _H->numericsMatrix().get();
  if (H_NM)
  {
    NM_display(H_NM);
  }

  std::cout << " - Vector q : " <<std::endl;
  if (_q) _q->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " - Vector b : " <<std::endl;
  if (_b) _b->display();
  else std::cout << "-> NULL" <<std::endl;

  std::cout << " - Vector z (reaction) : " <<std::endl;
  if (_z) _z->display();
  else std::cout << "-> NULL" <<std::endl;

  std::cout << " - Vector w (local velocities): " <<std::endl;
  if (_w) _w->display();
  else std::cout << "-> NULL" <<std::endl;

  std::cout << " - Vector globalVelocities : " <<std::endl;
  if (_globalVelocities) _globalVelocities->display();
  else std::cout << "-> NULL" <<std::endl;

  std::cout << "============================================================" <<std::endl;
}

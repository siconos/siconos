/* Siconos-Kernel, Copyright INRIA 2005-2012.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY ory FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
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
#include "NewtonImpactNSL.hpp"
#include "OSNSMatrix.hpp"

#include "gfc3d_Solvers.h"

 #define DEBUG_STDOUT
 #define DEBUG_MESSAGES
#include "debug.h"

// Constructor from a set of data
// Required input: simulation
// Optional: newNumericsSolverName
GlobalFrictionContact::GlobalFrictionContact(int dimPb, const int numericsSolverId):
  LinearOSNS(numericsSolverId), _contactProblemDim(dimPb)
{
}

GlobalFrictionContact::~GlobalFrictionContact()
{
  deleteSolverOptions(&*_numerics_solver_options);
}

void GlobalFrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

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

  // Memory allocation for reaction, and velocity
  initVectorsMemory();

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

  // get topology
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  SP::InteractionsGraph I0 = topology->indexSet0();
  _mu.reset(new MuStorage());
  _mu->reserve(I0->size());

  // Default size for M = _maxSize
  if (!_M)
  {
    if (MStorageType == 0)
      M.reset(new OSNSMatrix(_maxSize, 0));
    else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered DynamicalSystemsSet
      M.reset(new OSNSMatrix(simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), 1));
  }
  if (!_H)
  {
    if (_MStorageType == 0)
      _H.reset(new OSNSMatrix(_maxSize, 0));
    else // if(_MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered DynamicalSystemsSet
      _H.reset(new OSNSMatrix(simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), simulation->indexSet(_indexSetLevel)->size()   , 1));
  }


}

bool GlobalFrictionContact::preCompute(double time)
{
  // This function is used to prepare data for the GlobalFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  SP::Topology topology = simulation()->model()->nonSmoothDynamicalSystem()->topology();

  if (indexSetLevel() == LEVELMAX)
    return false;

  if (!_hasBeenUpdated)
  {
    InteractionsGraph& indexSet = *simulation()->model()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel);
    DynamicalSystemsGraph& DSG0 = *simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();

    _sizeOutput = 3*indexSet.size();

    if (_sizeOutput == 0) { return false; }

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

    size_t nnzM = 0;
    size_t nnzH = 0;
    size_t sizeM = 0;
    size_t nbDS = 0;



    // compute size and nnz of M and collect all matrices
    // compute nnz of H and collect H blocks
    // fill _b and mu
    if (_b->size() != _sizeOutput)
      _b->resize(_sizeOutput);

    size_t pos = 0;
    InteractionsGraph::VIterator ui, uiend;
    for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      Interaction& inter = *indexSet.bundle(*ui);
      VectorOfSMatrices& workMInter = *indexSet.properties(*ui).workMatrices;

      _mu->push_back(std11::static_pointer_cast<NewtonImpactFrictionNSL>(inter.nonSmoothLaw())->mu());

      OneStepIntegrator& Osi = *indexSet.properties(*ui).osi;
      OSI::TYPES osiType = Osi.getType();
      if (osiType == OSI::MOREAUJEANOSI2)
      {
        static_cast<MoreauJeanGOSI&>(Osi).NSLcontrib(inter, *this);
      }
      else
      {
        RuntimeException::selfThrow("GlobalFrictionContact::computeq. Not yet implemented for Integrator type : " + osiType);
      }
      setBlock(*inter.yForNSsolver(), _b, 3, 0, pos);
      nnzH += inter.getLeftInteractionBlock(workMInter).nnz();
      pos += 3;

      SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
      SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;
      SiconosMatrix* W = DSG0.properties(DSG0.descriptor(ds1)).W.get();
      assert(W);
      bool inserted = dsMat.insert(std::make_pair(ds1, W)).second;

      if (inserted) // first time we see this DS
      {
        absPosDS.insert(std::make_pair(ds1, sizeM));

        // update sizes
        sizeM += W->size(0);
        nnzM += W->nnz();
        ++nbDS;
      }
      if (ds1 != ds2)
      {
        SiconosMatrix* W = DSG0.properties(DSG0.descriptor(ds2)).W.get();
        assert(W);
        bool inserted = dsMat.insert(std::make_pair(ds2, W)).second;

        if (inserted) // first time we see this DS
        {
          absPosDS.insert(std::make_pair(ds2, sizeM));

          // update sizes
          sizeM += W->size(0);
          nnzM += W->nnz();
          ++nbDS;
        }
      }
    }

    // fill M and _q
    _sizeGlobalOutput = sizeM;
    if (_q->size() != _sizeGlobalOutput)
      _q->resize(_sizeGlobalOutput);

    NumericsMatrix& M_NM = *_M->getNumericsMatrix();

    M_NM.storageType = NM_SPARSE;
    M_NM.size0 = sizeM;
    M_NM.size1 = sizeM;
    NM_csc_alloc(&M_NM, nnzM);
    M_NM.matrix2->origin = NS_CSC;
    CSparseMatrix* Mcsc = NM_csc(&M_NM);
    Mcsc->p[0] = 0.;

    size_t offset = 0;
    for (dsMatMap::iterator it = dsMat.begin(); it != dsMat.end(); ++it)
    {
      SP::DynamicalSystem ds = (*it).first;
      SiconosMatrix& mat = *(*it).second;

      size_t dss = ds->getDim();

      // compute q (aka free velocity) = v^k + contribution from forces
      OneStepIntegrator& Osi = *DSG0.properties(DSG0.descriptor(ds)).osi;
      OSI::TYPES osiType = Osi.getType();
      if (osiType == OSI::MOREAUJEANOSI2)
      {
        setBlock(*ds->workspace(DynamicalSystem::free), _q, dss, 0, offset);
      }
      else
      {
          RuntimeException::selfThrow("GlobalFrictionContact::computeq. Not yet implemented for Integrator type : " + osiType);
      }

      // fill matrix
      mat.fillCSC(Mcsc, offset, offset);
      offset += dss;
    }


    // fill H
    NumericsMatrix& H_NM = *_H->getNumericsMatrix();

    H_NM.storageType = NM_SPARSE;
    H_NM.size0 = sizeM;
    H_NM.size1 = _sizeOutput;
    NM_csc_alloc(&H_NM, nnzH);
    H_NM.matrix2->origin = NS_CSC;
    CSparseMatrix* Hcsc = NM_csc(&H_NM);
    Hcsc->p[0] = 0;

    pos = 0;
    offset = 0;


    SP::SiconosMatrix leftInteractionBlock;
    for (std11::tie(ui, uiend) = indexSet.vertices(); ui != uiend; ++ui)
    {
      Interaction& inter = *indexSet.bundle(*ui);
      VectorOfSMatrices& workMInter = *indexSet.properties(*ui).workMatrices;

      SP::DynamicalSystem ds1 = indexSet.properties(*ui).source;
      SP::DynamicalSystem ds2 = indexSet.properties(*ui).target;

      bool endl = false;
      size_t posBlock = indexSet.properties(*ui).source_pos;
      size_t pos2 = indexSet.properties(*ui).target_pos;
      for (SP::DynamicalSystem ds = ds1; !endl; ds = ds2, posBlock = pos2)
      {
        endl = (ds == ds2);
        size_t sizeDS = ds->getDim();
        // this whole part is a hack. Just should just get the rightblock
        leftInteractionBlock.reset(new SimpleMatrix(3, sizeDS));
        inter.getLeftInteractionBlockForDS(posBlock, leftInteractionBlock, workMInter);
        leftInteractionBlock->trans();
        leftInteractionBlock->fillCSC(Hcsc, absPosDS[ds], pos);
      }
      pos += 3;


    }
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
  if (_sizeOutput != 0)
  {
    // The GlobalFrictionContact Problem in Numerics format
    GlobalFrictionContactProblem numerics_problem;
    globalFrictionContact_null(&numerics_problem);
    numerics_problem.M = &*_M->getNumericsMatrix();
    numerics_problem.H = &*_H->getNumericsMatrix();
    numerics_problem.q = _q->getArray();
    numerics_problem.b = _b->getArray();
    numerics_problem.numberOfContacts = _sizeOutput / _contactProblemDim;
    numerics_problem.mu = &(_mu->at(0));
    numerics_problem.dimension = 3;
    info = (*_gfc_driver)(&numerics_problem,
                           _z->getArray(),
                           _w->getArray(),
                           _globalVelocities->getArray(),
                           &*_numerics_solver_options,
                           &*_numerics_options);
    postCompute();

  }

  return info;
}

void GlobalFrictionContact::postCompute()
{
  // This function is used to set y/lambda values using output from primalfrictioncontact_driver
  // Only Interactions (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  InteractionsGraph& indexSet = *simulation()->model()->nonSmoothDynamicalSystem()->topology()->indexSet(_indexSetLevel);
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

//  SP::DynamicalSystemsGraph DSG = simulation()->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
//  DSIterator itDS;
//  unsigned int sizeDS;
//  SP::OneStepIntegrator  Osi;
//  std::string osiType; // type of the current one step integrator
//  std::string dsType; // type of the current Dynamical System
  //   for(itDS = allDS->begin(); itDS!=  allDS->end(); ++itDS)
  //     {
  //       dsType = (*itDS) -> getType();
  //       if(dsType!=Type::LagrangianDS && dsType!=Type::LagrangianLinearTIDS)
  //      RuntimeException::selfThrow("GlobalFrictionContact::postCompute not yet implemented for dynamical system of types "+dsType);

  //       pos = M->getPositionOfDSBlock(*itDS);
  //       sizeDS = (*itDS)->getDim();
  //       setBlock((static_cast<LagrangianDS*>(*itDS))->velocity(),velocity,sizeDS, pos, 0 );

  //     }



}

void GlobalFrictionContact::display() const
{
  std::cout << "===== " << _contactProblemDim << "D Primal Friction Contact Problem " <<std::endl;
  std::cout << "size (_sizeOutput) " << _sizeOutput <<std::endl;
  std::cout << "and  size (_sizeGlobalOutput) " << _sizeGlobalOutput << "(ie " << _sizeGlobalOutput / _contactProblemDim << " contacts)." <<std::endl;
  std::cout << " - Matrix M  : " <<std::endl;
  if (_M) _M->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " - Matrix H : " <<std::endl;
  if (_H) _H->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " - Vector q : " <<std::endl;
  if (_q) _q->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "============================================================" <<std::endl;
}

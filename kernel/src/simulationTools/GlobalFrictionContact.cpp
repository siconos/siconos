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
#include "MoreauJeanOSI.hpp" // Numerics Header
#include "LagrangianDS.hpp"
#include "NewtonImpactNSL.hpp"


// Constructor from a set of data
// Required input: simulation
// Optional: newNumericsSolverName
GlobalFrictionContact::GlobalFrictionContact(int dimPb,
    const int numericsSolverId,
    const std::string& newId):
  LinearOSNS(numericsSolverId, "GlobalFrictionContact", newId), _contactProblemDim(dimPb)
{}

void GlobalFrictionContact::setLocalVelocity(const SiconosVector& newValue)
{
  assert(0);
}

void GlobalFrictionContact::setLocalReaction(const SiconosVector& newValue)
{
  assert(0);
}

void GlobalFrictionContact::setTildeLocalVelocity(const SiconosVector& newValue)
{
  assert(0);
}

void GlobalFrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  // updateInteractionDSBlocks(); This computation is not made because, we that InteractionDSBlocks =H^T
  updateInteractionDSBlocks(); // blocks of H

  // Connect to the right function according to dim. of the problem
  if (_contactProblemDim == 2)
    ;
  //    primalFrictionContact_driver = &pfc_2D_driver;
  else // if(_contactProblemDim == 3)
    //primalFrictionContact_driver = &frictionContact3D_driver;
    ;

  // Memory allocation for reaction, and velocity
  initVectorsMemory();

  if (!_localVelocity)
    _localVelocity.reset(new SiconosVector(_maxSize));
  else
  {
    if (_localVelocity->size() != _maxSize)
      _localVelocity->resize(_maxSize);
  }

  if (!_localReaction)
    _localReaction.reset(new SiconosVector(_maxSize));
  else
  {
    if (_localReaction->size() != _maxSize)
      _localReaction->resize(_maxSize);
  }

  if (!_tildeLocalVelocity)
    _tildeLocalVelocity.reset(new SiconosVector(_maxSize));
  else
  {
    if (_tildeLocalVelocity->size() != _maxSize)
      _tildeLocalVelocity->resize(_maxSize);
  }

  // get topology
  SP::Topology topology = simulation()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  SP::InteractionsGraph I0 = topology->indexSet0();
  _mu.reset(new MuStorage());
  _mu->reserve(I0->size());

  SP::DynamicalSystemsSet allDS = simulation->nonSmoothDynamicalSystem()->dynamicalSystems();;


  // Default size for M = _maxSize
  if (!M)
  {
    if (MStorageType == 0)
      M.reset(new OSNSMatrix(_maxSize, 0));
    else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered DynamicalSystemsSet
      M.reset(new OSNSMatrix(simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), 1));
  }
  if (!H)
  {
    if (MStorageType == 0)
      H.reset(new OSNSMatrix(_maxSize, 0));
    else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered DynamicalSystemsSet


      H.reset(new OSNSMatrix(simulation->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), simulation->indexSet(levelMin)->size()   , 1));
  }


}

void GlobalFrictionContact::computeInteractionBlock(SP::Interaction inter1, SP::Interaction inter2)
{
  // Computes matrix interactionBlocks[inter1][inter2] (and allocates memory if necessary) if inter1 and inter2 have commond DynamicalSystem.
  // How interactionBlocks are computed depends explicitely on the type of Relation of each Interaction.

}

void GlobalFrictionContact::computeDSBlock(SP::DynamicalSystem DS)
{
  // Computes matrix DSBlocks[DS1](and allocates memory if necessary)

  DSIterator itDS;
  SP::OneStepIntegrator  Osi;
  std::string osiType; // type of the current one step integrator
  std::string dsType; // type of the current Dynamical System;

  Osi = simulation->integratorOfDS(DS); // get OneStepIntegrator of current dynamical system
  osiType = Osi->getType();
  if (osiType == MOREAU || osiType == MOREAU2)
  {
    DSBlocks[DS] = (std11::static_pointer_cast<MoreauJeanOSI> (Osi))->W(DS); // get its W matrix ( pointer link!)
    //       std::cout << "GlobalFrictionContact::computeDSBlock(SP::DynamicalSystem DS) " <<std::endl;
    //       DSBlocks[DS]->display();
  }
  else if (osiType == LSODAR) // Warning: LagrangianDS only at the time !!!
  {
    RuntimeException::selfThrow("GlobalFrictionContact::computeDSBlocks. Not yet implemented for LsodarOSI Integrator");
  }
  else
    RuntimeException::selfThrow("GlobalFrictionContact::computeDSBlocks. nNot yet implemented for Integrator of type " + osiType);
}

/** computes  DSInteractionBlock-matrix that corresponds to inter1 and DS2
 *  Move this to Interaction class?
 *  \param a pointer to Interaction inter1
 *  \param a pointer to DynamicalSystems DS2
 */
void GlobalFrictionContact::computeInteractionDSBlock(SP::Interaction inter, SP::DynamicalSystem DS)
{
  unsigned int sizeDS = (DS)->getDim();
  unsigned int nslawSize = inter->nonSmoothLaw()->size();
  RELATION::TYPES relationType = inter->relation()->getType();

  if (relationType == Lagrangian)
  {
    interactionDSBlocks[inter][DS].reset(new SimpleMatrix(nslawSize, sizeDS));
    inter->getLeftInteractionBlockForDS(DS, interactionDSBlocks[inter][DS]);
  }
  else RuntimeException::selfThrow("GlobalFrictionContact::computeInteractionDSBlock not yet implemented for relation of type " + relationType);
}

void GlobalFrictionContact::computeq(double time)
{
  if (q->size() != _sizeOutput)
    q->resize(_sizeOutput);
  q->zero();

  SP::DynamicalSystemsSet allDS = simulation->nonSmoothDynamicalSystem()->dynamicalSystems();;
  DSIterator itDS;

  // type of the current one step integrator
  unsigned int pos = 0;
  for (itDS = allDS->begin(); itDS !=  allDS->end(); ++itDS)
  {
    pos = M->getPositionOfDSBlock(*itDS);
    computeqBlockDS((*itDS), pos);
  }
}

void GlobalFrictionContact::computeqBlockDS(SP::DynamicalSystem ds, unsigned int pos)
{
  SP::OneStepIntegrator  Osi = simulation->integratorOfDS(ds);
  std::string osiType = Osi->getType();
  unsigned int sizeDS;
  if (osiType == MOREAU2)
  {
    sizeDS = (ds)->getDim();
    setBlock(Osi->getWorkX(ds), q, sizeDS, 0, pos);
  }
  else
  {
    RuntimeException::selfThrow("GlobalFrictionContact::computeq. Not yet implemented for Integrator type : " + osiType);
  }
}

void GlobalFrictionContact::computeTildeLocalVelocityBlock(SP::Interaction inter, unsigned int pos)
{
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();
  std::string nslawType = inter->getNonSmoothLawType();

  std::string simulationType = simulation->getType();

  SP::DynamicalSystem ds = *(inter->dynamicalSystemsBegin());
  std::string osiType = simulation->integratorOfDS(ds)->getType();

  unsigned int sizeY = inter->nonSmoothLaw()->size();
  Index coord(8);

  unsigned int relativePosition = 0;
  SP::Interaction mainInteraction = inter;
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = pos;
  coord[7] = pos + sizeY;

  SP::SiconosMatrix  H;
  SP::SiconosVector workX = inter->workX;
  if (osiType == MOREAU2)
  {
  }
  else
    RuntimeException::selfThrow("GlobalFrictionContact::computeTildeLocalVelocityBlock not yet implemented for OSI of type " + osiType);

  // Add "non-smooth law effect" on q
  if (inter->relation()->getType() == Lagrangian)
  {
    double e;
    if (nslawType == NEWTONIMPACTNSLAW)
    {

      e = (std11::static_pointer_cast<NewtonImpactNSL>(mainInteraction->nonSmoothLaw()))->getE();

      Index subCoord(4);
      if (simulationType == "TimeStepping")
      {
        subCoord[0] = 0;
        subCoord[1] = inter->nonSmoothLaw()->size();
        subCoord[2] = pos;
        subCoord[3] = pos + subCoord[1];
        subscal(e, *inter->yOld(levelMin), *q, subCoord, false);
      }
      else if (simulationType == "EventDriven")
      {
        subCoord[0] = pos;
        subCoord[1] = pos + inter->nonSmoothLaw()->size();
        subCoord[2] = pos;
        subCoord[3] = subCoord[1];
        subscal(e, *_tildeLocalVelocity, *_tildeLocalVelocity, subCoord, false); // _tildeLocalVelocity = _tildeLocalVelocity + e * _tildeLocalVelocity
      }
      else
        RuntimeException::selfThrow("FrictionContact::computetildeLocalVelocityBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {

      e = (std11::static_pointer_cast<NewtonImpactFrictionNSL>(mainInteraction->nonSmoothLaw()))->getEn();

      // Only the normal part is multiplied by e
      if (simulationType == "TimeStepping")
        (*_tildeLocalVelocity)(pos) +=  e * (*inter->yOld(levelMin))(0);

      else RuntimeException::selfThrow("FrictionContact::computetildeLocalVelocityBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);

    }
    else
      RuntimeException::selfThrow("FrictionContact::computeTILDELOCALVELOCITYBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType);
  }


}
void GlobalFrictionContact::computeTildeLocalVelocity(double time)
{
  if (_tildeLocalVelocity->size() != _sizeLocalOutput)
    _tildeLocalVelocity->resize(_sizeLocalOutput);
  _tildeLocalVelocity->zero();
  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(levelMin);
  SP::Interaction inter;
  for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
  {
    inter = indexSet->bundle(*ui);
    pos = H->getPositionOfInteractionBlock(*inter);
    computeTildeLocalVelocityBlock(inter, pos);
  }

  // // === Loop through "active" Interactions (ie present in indexSets[level]) ===

  // unsigned int pos = 0;
  // InteractionsIterator itCurrent;
  // for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  // {
  //   pos = H->getPositionOfInteractionBlock(*itCurrent);
  //   computeTildeLocalVelocityBlock((*itCurrent), pos);
  // }

}

bool GlobalFrictionContact::preCompute(double time)
{
  // This function is used to prepare data for the GlobalFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  SP::Topology topology = simulation->nonSmoothDynamicalSystem()->topology();
  SP::DynamicalSystemsSet allDS = simulation->nonSmoothDynamicalSystem()->dynamicalSystems();;
  
  if (!_hasBeenUpdated)
  {
    // Computes new interactionBlocks if required
    updateInteractionBlocks();
    // updateInteractionDSBlocks(); This computation is not made because, we that InteractionDSBlocks =H^T
    updateInteractionDSBlocks(); // blocks of H



    // Updates matrix M
    SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(levelMin);
    M->fill(allDS, DSBlocks);
    // Note FP: old version ==
    // H->fill(indexSet, allDS, interactionDSBlocks);
    // fill function with this signature does not exists. How can it works????
    H->fill(indexSet); 
    _sizeOutput = M->size();
    _sizeLocalOutput = H->size();

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
    // Checks z and w sizes and reset if necessary
    if (_localReaction->size() != _sizeLocalOutput)
    {
      _localReaction->resize(_sizeLocalOutput, false);
      _localReaction->zero();
    }

    if (_localVelocity->size() != _sizeLocalOutput)
    {
      _localVelocity->resize(_sizeLocalOutput);
      _localVelocity->zero();
    }
 
    // Update mu
    _mu->clear();
    InteractionsGraph::VIterator ui, uiend;
    for (std11::tie(ui, uiend) = indexSet->vertices(); ui != uiend; ++ui)
    {
      _mu->push_back(std11::static_pointer_cast<NewtonImpactFrictionNSL>((indexSet->bundle(*ui))->nonSmoothLaw())->getMu());
    }

    //for (ConstInteractionsIterator itI = indexSet->begin(); itI != indexSet->end(); ++itI)
    // _mu->push_back(std11::static_pointer_cast<NewtonImpactFrictionNSL>((*itI)->nonSmoothLaw())->getMu());
  }

  // Computes q
  computeq(time);
  computeTildeLocalVelocity(time);

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
    numerics_problem.M = M->getNumericsMatrix();
    numerics_problem.q = q->getArray();
    numerics_problem.numberOfContacts = _sizeOutput / _contactProblemDim;
    numerics_problem.mu = &(_mu->at(0));
    // Call Numerics Driver for GlobalFrictionContact
    //  {
    //    int info2 = -1;
    //    int iparam2[7];
    //    iparam2[0] = 100000;
    //    iparam2[1] = 1;
    //    iparam2[3] = iparam2[0];
    //    iparam2[5] = 0;
    //    iparam2[6] = 0;// 0: proj 1: AC/FB

    //    double dparam2[5];
    //    dparam2[0] = 1e-6;
    //    dparam2[2] = 1e-6;
    //    SiconosVector * z= new SiconosVector(_sizeOutput);
    //    SiconosVector * w= new SiconosVector(_sizeOutput);
    //    M->display();
    //    q->display();

    //    pfc_3D_nsgs(numerics_problem.numberOfContacts, numerics_problem.M->matrix0, numerics_problem.q , z->getArray() , w->getArray(), mu->data(), &info2, iparam2 , dparam2 );
    //    std::cout << " ... " << info2 << " " << dparam2[1] <<std::endl;
    //    z->display();
    //    w->display();


    //  }
    //       M->display();
    //       q->display();
    //      info = (*primalFrictionContact_driver)(&numerics_problem, reaction->getArray() , velocity->getArray() , solver->numericsSolverOptions(), numerics_options);
    //  std::cout << " step 2 "<< info <<std::endl;
    //  _z->display();
    //  velocity->display();
    // --- Recovering of the desired variables from LCP output ---
    postCompute();

  }

  return info;
}

void GlobalFrictionContact::postCompute()
{
  // This function is used to set y/lambda values using output from primalfrictioncontact_driver
  // Only Interactions (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::InteractionsGraph indexSet = model()->nonSmoothDynamicalSystem()->topology()->indexSet(levelMin);
  // y and lambda vectors
  SP::SiconosVector  y, lambda;

  //   // === Loop through "active" Interactions (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;

  //   for(InteractionsIterator itCurrent = indexSet->begin(); itCurrent!=  indexSet->end(); ++itCurrent)
  //     {
  //       // size of the interactionBlock that corresponds to the current Interaction
  //       nsLawSize = (*itCurrent)->nonSmoothLaw()->size();
  //       // Get the relative position of inter-interactionBlock in the vector velocity or reaction
  //       pos = H->getPositionOfInteractionBlock(*itCurrent);

  //       // Get Y and Lambda for the current Interaction
  //       y = (*itCurrent)->y(levelMin);
  //       lambda = (*itCurrent)->lambda(levelMin);

  //       // Copy velocity/_z values, starting from index pos into y/lambda.

  //       setBlock(_localVelocity, y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
  //       setBlock(_localReaction, lambda, lambda->size(), pos, 0);

  //     }

  SP::DynamicalSystemsSet allDS = simulation->nonSmoothDynamicalSystem()->dynamicalSystems();;
  DSIterator itDS;
  unsigned int sizeDS;
  SP::OneStepIntegrator  Osi;
  std::string osiType; // type of the current one step integrator
  std::string dsType; // type of the current Dynamical System
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
  std::cout << "and  size (_sizeLocalOutput) " << _sizeLocalOutput << "(ie " << _sizeLocalOutput / _contactProblemDim << " contacts)." <<std::endl;
  std::cout << " - Matrix M  : " <<std::endl;
  if (M) M->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " - Matrix H : " <<std::endl;
  if (H) H->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " - Vector q : " <<std::endl;
  if (q) q->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " - Vector _tildeLocalVelocity : " <<std::endl;
  if (_tildeLocalVelocity) _tildeLocalVelocity->display();
  else std::cout << "-> NULL" <<std::endl;
  std::cout << " Friction coefficients: " <<std::endl;
  if (_mu) print(_mu->begin(), _mu->end());
  else std::cout << "-> NULL" <<std::endl;
  std::cout << "============================================================" <<std::endl;
}

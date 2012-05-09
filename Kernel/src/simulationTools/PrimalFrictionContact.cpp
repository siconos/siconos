/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
#include "PrimalFrictionContact.hpp"
#include "FrictionContactXML.hpp"
#include "Simulation.hpp"
//#include "Interaction.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Relation.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "Moreau.hpp" // Numerics Header
#include "LagrangianDS.hpp"
#include "NewtonImpactNSL.hpp"
using namespace std;

// xml constructor
PrimalFrictionContact::PrimalFrictionContact(SP::OneStepNSProblemXML osNsPbXml):
  LinearOSNS(osNsPbXml, "PrimalFrictionContact"), _contactProblemDim(3)
{
  SP::FrictionContactXML xmlFC = boost::static_pointer_cast<FrictionContactXML>(osNsPbXml);

  // Read dimension of the problem (required parameter)
  if (!xmlFC->hasProblemDim())
    RuntimeException::selfThrow(
      "PrimalFrictionContact: xml constructor failed, attribute for dimension of the problem (2D or 3D) is missing.");

  _contactProblemDim = xmlFC->getProblemDim();
}

// Constructor from a set of data
// Required input: simulation
// Optional: newNumericsSolverName
PrimalFrictionContact::PrimalFrictionContact(int dimPb,
    const int newNumericsSolverId,
    const string& newId):
  LinearOSNS(newNumericsSolverId, "PrimalFrictionContact", newId), _contactProblemDim(dimPb)
{}

void PrimalFrictionContact::setLocalVelocity(const SiconosVector& newValue)
{
  assert(0);
}

void PrimalFrictionContact::setLocalReaction(const SiconosVector& newValue)
{
  assert(0);
}

void PrimalFrictionContact::setTildeLocalVelocity(const SiconosVector& newValue)
{
  assert(0);
}

void PrimalFrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  updateDSBlocks(); //blocks of M
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
    _localVelocity.reset(new SimpleVector(_maxSize));
  else
  {
    if (_localVelocity->size() != _maxSize)
      _localVelocity->resize(_maxSize);
  }

  if (!_localReaction)
    _localReaction.reset(new SimpleVector(_maxSize));
  else
  {
    if (_localReaction->size() != _maxSize)
      _localReaction->resize(_maxSize);
  }

  if (!_tildeLocalVelocity)
    _tildeLocalVelocity.reset(new SimpleVector(_maxSize));
  else
  {
    if (_tildeLocalVelocity->size() != _maxSize)
      _tildeLocalVelocity->resize(_maxSize);
  }

  // get topology
  SP::Topology topology = simulation()->model()->nonSmoothDynamicalSystem()->topology();

  // Note that interactionBlocks is up to date since updateInteractionBlocks has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  SP::InteractionsGraph I0 = topology->indexSet0();
  _mu.reset(new MuStorage());
  _mu->reserve(I0->size());

  SP::DynamicalSystemsSet allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();;


  // Default size for M = _maxSize
  if (!M)
  {
    if (MStorageType == 0)
      M.reset(new OSNSMatrix(_maxSize, 0));
    else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered DynamicalSystemsSet
      M.reset(new OSNSMatrix(simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), 1));
  }
  if (!H)
  {
    if (MStorageType == 0)
      H.reset(new OSNSMatrix(_maxSize, 0));
    else // if(MStorageType == 1) size = number of DSBlocks = number of DS in the largest considered DynamicalSystemsSet


      H.reset(new OSNSMatrix(simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems()->size(), simulation->indexSet(levelMin)->size()   , 1));
  }


}

void PrimalFrictionContact::computeInteractionBlock(SP::Interaction inter1, SP::Interaction inter2)
{
  // Computes matrix interactionBlocks[inter1][inter2] (and allocates memory if necessary) if inter1 and inter2 have commond DynamicalSystem.
  // How interactionBlocks are computed depends explicitely on the type of Relation of each Interaction.

}

void PrimalFrictionContact::computeDSBlock(SP::DynamicalSystem DS)
{
  // Computes matrix DSBlocks[DS1](and allocates memory if necessary)

  DSIterator itDS;
  SP::OneStepIntegrator  Osi;
  string osiType; // type of the current one step integrator
  string dsType; // type of the current Dynamical System;

  Osi = simulation->integratorOfDS(DS); // get OneStepIntegrator of current dynamical system
  osiType = Osi->getType();
  if (osiType == MOREAU || osiType == MOREAU2)
  {
    DSBlocks[DS] = (boost::static_pointer_cast<Moreau> (Osi))->W(DS); // get its W matrix ( pointer link!)
    //       cout << "PrimalFrictionContact::computeDSBlock(SP::DynamicalSystem DS) " << endl;
    //       DSBlocks[DS]->display();
  }
  else if (osiType == LSODAR) // Warning: LagrangianDS only at the time !!!
  {
    RuntimeException::selfThrow("PrimalFrictionContact::computeDSBlocks. Not yet implemented for Lsodar Integrator");
  }
  else
    RuntimeException::selfThrow("PrimalFrictionContact::computeDSBlocks. nNot yet implemented for Integrator of type " + osiType);
}

/** computes  DSInteractionBlock-matrix that corresponds to inter1 and DS2
 *  Move this to Interaction class?
 *  \param a pointer to Interaction inter1
 *  \param a pointer to DynamicalSystems DS2
 */
void PrimalFrictionContact::computeInteractionDSBlock(SP::Interaction inter, SP::DynamicalSystem DS)
{
  unsigned int sizeDS = (DS)->getDim();
  unsigned int nslawSize = inter->getNonSmoothLawSize();
  RELATION::TYPES relationType = inter->getRelationType();

  if (relationType == Lagrangian)
  {
    interactionDSBlocks[inter][DS].reset(new SimpleMatrix(nslawSize, sizeDS));
    inter->getLeftInteractionBlockForDS(DS, interactionDSBlocks[inter][DS]);
  }
  else RuntimeException::selfThrow("PrimalFrictionContact::computeInteractionDSBlock not yet implemented for relation of type " + relationType);
}

void PrimalFrictionContact::computeq(const double time)
{
  if (q->size() != _sizeOutput)
    q->resize(_sizeOutput);
  q->zero();

  SP::DynamicalSystemsSet allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();;
  DSIterator itDS;

  // type of the current one step integrator
  unsigned int pos = 0;
  for (itDS = allDS->begin(); itDS !=  allDS->end(); ++itDS)
  {
    pos = M->getPositionOfDSBlock(*itDS);
    computeqBlock((*itDS), pos);
  }
}

void PrimalFrictionContact::computeqBlock(SP::DynamicalSystem DS, unsigned int pos)
{
  SP::OneStepIntegrator  Osi = simulation->integratorOfDS(DS);
  string osiType = Osi->getType();
  unsigned int sizeDS;
  if (osiType == MOREAU2)
  {
    sizeDS = (DS)->getDim();
    setBlock(Osi->getWorkX(DS), q, sizeDS, 0, pos);
  }
  else
  {
    RuntimeException::selfThrow("PrimalFrictionContact::computeq. Not yet implemented for Integrator type : " + osiType);
  }
}
void PrimalFrictionContact::computeTildeLocalVelocityBlock(SP::Interaction inter, unsigned int pos)
{
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->getRelationType();
  RELATION::SUBTYPES relationSubType = inter->getRelationSubType();
  string nslawType = inter->getNonSmoothLawType();

  string simulationType = simulation->getType();

  SP::DynamicalSystem ds = *(inter->dynamicalSystemsBegin());
  string osiType = simulation->integratorOfDS(ds)->getType();

  unsigned int sizeY = inter->getNonSmoothLawSize();
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
    RuntimeException::selfThrow("PrimalFrictionContact::computeTildeLocalVelocityBlock not yet implemented for OSI of type " + osiType);

  // Add "non-smooth law effect" on q
  if (inter->getRelationType() == Lagrangian)
  {
    double e;
    if (nslawType == NEWTONIMPACTNSLAW)
    {

      e = (boost::static_pointer_cast<NewtonImpactNSL>(mainInteraction->nonSmoothLaw()))->getE();

      Index subCoord(4);
      if (simulationType == "TimeStepping")
      {
        subCoord[0] = 0;
        subCoord[1] = inter->getNonSmoothLawSize();
        subCoord[2] = pos;
        subCoord[3] = pos + subCoord[1];
        subscal(e, *inter->yOld(levelMin), *q, subCoord, false);
      }
      else if (simulationType == "EventDriven")
      {
        subCoord[0] = pos;
        subCoord[1] = pos + inter->getNonSmoothLawSize();
        subCoord[2] = pos;
        subCoord[3] = subCoord[1];
        subscal(e, *_tildeLocalVelocity, *_tildeLocalVelocity, subCoord, false); // _tildeLocalVelocity = _tildeLocalVelocity + e * _tildeLocalVelocity
      }
      else
        RuntimeException::selfThrow("FrictionContact::computetildeLocalVelocityBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {

      e = (boost::static_pointer_cast<NewtonImpactFrictionNSL>(mainInteraction->nonSmoothLaw()))->getEn();

      // Only the normal part is multiplied by e
      if (simulationType == "TimeStepping")
        (*_tildeLocalVelocity)(pos) +=  e * (*inter->yOld(levelMin))(0);

      else RuntimeException::selfThrow("FrictionContact::computetildeLocalVelocityBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);

    }
    else
      RuntimeException::selfThrow("FrictionContact::computeTILDELOCALVELOCITYBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType);
  }


}
void PrimalFrictionContact::computeTildeLocalVelocity(const double time)
{
  if (_tildeLocalVelocity->size() != _sizeLocalOutput)
    _tildeLocalVelocity->resize(_sizeLocalOutput);
  _tildeLocalVelocity->zero();
  // === Get index set from Simulation ===
  SP::InteractionsSet indexSet = simulation->indexSet(levelMin);
  // === Loop through "active" Interactions (ie present in indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsIterator itCurrent;
  for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    pos = H->getPositionOfInteractionBlock(*itCurrent);
    computeTildeLocalVelocityBlock((*itCurrent), pos);
  }

}

void PrimalFrictionContact::preCompute(const double time)
{
  // This function is used to prepare data for the PrimalFrictionContact problem
  // - computation of M, H _tildeLocalVelocity and q
  // - set _sizeOutput, sizeLocalOutput

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, _sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  SP::Topology topology = simulation->model()->nonSmoothDynamicalSystem()->topology();
  SP::DynamicalSystemsSet allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();;
  SP::InteractionsSet indexSet = simulation->indexSet(levelMin);

  if (!_hasBeenUpdated)
  {
    // Computes new interactionBlocks if required
    updateInteractionBlocks();
    updateDSBlocks(); //blocks of M
    // updateInteractionDSBlocks(); This computation is not made because, we that InteractionDSBlocks =H^T
    updateInteractionDSBlocks(); // blocks of H



    // Updates matrix M
    SP::InteractionsSet indexSet = simulation->indexSet(levelMin);
    M->fill(allDS, DSBlocks);
    H->fill(indexSet, allDS, interactionDSBlocks);
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
    for (ConstInteractionsIterator itI = indexSet->begin(); itI != indexSet->end(); ++itI)
      _mu->push_back(boost::static_pointer_cast<NewtonImpactFrictionNSL>((*itI)->nonSmoothLaw())->getMu());
  }

  // Computes q
  computeq(time);
  computeTildeLocalVelocity(time);

}

int PrimalFrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for PrimalFrictionContact computing ---
  preCompute(time);

  // --- Call Numerics solver ---
  if (_sizeOutput != 0)
  {
    // The PrimalFrictionContact Problem in Numerics format
    PrimalFrictionContactProblem numerics_problem;
    numerics_problem.M = M->getNumericsMatrix();
    numerics_problem.q = q->getArray();
    numerics_problem.numberOfContacts = _sizeOutput / _contactProblemDim;
    numerics_problem.mu = &((*_mu)[0]);
    // Call Numerics Driver for PrimalFrictionContact
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
    //    SimpleVector * z= new SimpleVector(_sizeOutput);
    //    SimpleVector * w= new SimpleVector(_sizeOutput);
    //    M->display();
    //    q->display();

    //    pfc_3D_nsgs(numerics_problem.numberOfContacts, numerics_problem.M->matrix0, numerics_problem.q , z->getArray() , w->getArray(), mu->data(), &info2, iparam2 , dparam2 );
    //    cout << " ... " << info2 << " " << dparam2[1] << endl;
    //    z->display();
    //    w->display();


    //  }
    //       M->display();
    //       q->display();
    //      info = (*primalFrictionContact_driver)(&numerics_problem, reaction->getArray() , velocity->getArray() , solver->numericsSolverOptions(), numerics_options);
    //  cout << " step 2 "<< info << endl;
    //  _z->display();
    //  velocity->display();
    // --- Recovering of the desired variables from LCP output ---
    postCompute();

  }

  return info;
}

void PrimalFrictionContact::postCompute()
{
  // This function is used to set y/lambda values using output from primalfrictioncontact_driver
  // Only Interactions (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::InteractionsSet indexSet = simulation->indexSet(levelMin);

  // y and lambda vectors
  SP::SiconosVector  y, lambda;

  //   // === Loop through "active" Interactions (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;

  //   for(InteractionsIterator itCurrent = indexSet->begin(); itCurrent!=  indexSet->end(); ++itCurrent)
  //     {
  //       // size of the interactionBlock that corresponds to the current Interaction
  //       nsLawSize = (*itCurrent)->getNonSmoothLawSize();
  //       // Get the relative position of inter-interactionBlock in the vector velocity or reaction
  //       pos = H->getPositionOfInteractionBlock(*itCurrent);

  //       // Get Y and Lambda for the current Interaction
  //       y = (*itCurrent)->y(levelMin);
  //       lambda = (*itCurrent)->lambda(levelMin);

  //       // Copy velocity/_z values, starting from index pos into y/lambda.

  //       setBlock(_localVelocity, y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
  //       setBlock(_localReaction, lambda, lambda->size(), pos, 0);

  //     }

  SP::DynamicalSystemsSet allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();;
  DSIterator itDS;
  unsigned int sizeDS;
  SP::OneStepIntegrator  Osi;
  string osiType; // type of the current one step integrator
  string dsType; // type of the current Dynamical System
  //   for(itDS = allDS->begin(); itDS!=  allDS->end(); ++itDS)
  //     {
  //       dsType = (*itDS) -> getType();
  //       if(dsType!=Type::LagrangianDS && dsType!=Type::LagrangianLinearTIDS)
  //      RuntimeException::selfThrow("PrimalFrictionContact::postCompute not yet implemented for dynamical system of types "+dsType);

  //       pos = M->getPositionOfDSBlock(*itDS);
  //       sizeDS = (*itDS)->getDim();
  //       setBlock((static_cast<LagrangianDS*>(*itDS))->velocity(),velocity,sizeDS, pos, 0 );

  //     }



}

void PrimalFrictionContact::display() const
{
  cout << "===== " << _contactProblemDim << "D Primal Friction Contact Problem " << endl;
  cout << "size (_sizeOutput) " << _sizeOutput << endl;
  cout << "and  size (_sizeLocalOutput) " << _sizeLocalOutput << "(ie " << _sizeLocalOutput / _contactProblemDim << " contacts)." << endl;
  cout << " - Matrix M  : " << endl;
  if (M) M->display();
  else cout << "-> NULL" << endl;
  cout << " - Matrix H : " << endl;
  if (H) H->display();
  else cout << "-> NULL" << endl;
  cout << " - Vector q : " << endl;
  if (q) q->display();
  else cout << "-> NULL" << endl;
  cout << " - Vector _tildeLocalVelocity : " << endl;
  if (_tildeLocalVelocity) _tildeLocalVelocity->display();
  else cout << "-> NULL" << endl;
  cout << " Friction coefficients: " << endl;
  if (_mu) print(_mu->begin(), _mu->end());
  else cout << "-> NULL" << endl;
  cout << "============================================================" << endl;
}

PrimalFrictionContact* PrimalFrictionContact::convert(OneStepNSProblem* osnsp)
{
  PrimalFrictionContact* fc2d = dynamic_cast<PrimalFrictionContact*>(osnsp);
  return fc2d;
}



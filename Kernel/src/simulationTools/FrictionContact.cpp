/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "FrictionContact.h"
#include "FrictionContactXML.h"
#include "Topology.h"
#include "Interaction.h"
#include "Simulation.h"
#include "UnitaryRelation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "DynamicalSystem.h"
#include "Relation.h"
#include "NewtonImpactFrictionNSL.h"
#include "FrictionContact_Problem.h" // Numerics Header
#include "RelationTypes.h"
#include "NewtonImpactNSL.h"
#include "OneStepIntegrator.h"
using namespace std;
using namespace RELATION;
// xml constructor
FrictionContact::FrictionContact(SP::OneStepNSProblemXML osNsPbXml):
  OneStepNSProblem("FrictionContact", osNsPbXml), contactProblemDim(3)
{
  SP::FrictionContactXML  xmlFC = (boost::static_pointer_cast<FrictionContactXML>(osNsPbXml));

  // Read dimension of the problem (required parameter)
  if (!xmlFC->hasProblemDim())
    RuntimeException::selfThrow("FrictionContact: xml constructor failed, attribute for dimension of the problem (2D or 3D) is missing.");

  contactProblemDim = xmlFC->getProblemDim();

  // Read storage type if given (optional , default = dense)
  if (onestepnspbxml->hasStorageType())
    MStorageType = onestepnspbxml->getStorageType();

}

// Constructor from a set of data
// Required input: simulation
// Optional: NonSmoothSolver and id
FrictionContact::FrictionContact(int dimPb, SP::NonSmoothSolver  newSolver, const string& newId):
  OneStepNSProblem("FrictionContact", newId, newSolver), contactProblemDim(dimPb)
{}

// Setters

void FrictionContact::setVelocity(const SimpleVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("FrictionContact: setVelocity, inconsistent size between given velocity size and problem size. You should set sizeOutput before");

  if (! velocity)
    velocity.reset(new SimpleVector(sizeOutput));
  else if (velocity->size() != sizeOutput)
    RuntimeException::selfThrow("FrictionContact: setVelocity, velocity size differs from sizeOutput");

  *velocity = newValue;
}

void FrictionContact::setVelocityPtr(SP::SimpleVector newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("FrictionContact: setVelocityPtr, inconsistent size between given velocity size and problem size. You should set sizeOutput before");

  velocity = newPtr;

}


void FrictionContact::setReaction(const SimpleVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("FrictionContact: setReaction, inconsistent size between given reaction size and problem size. You should set sizeOutput before");

  if (! reaction)
    reaction.reset(new SimpleVector(sizeOutput));

  *reaction = newValue;
}

void FrictionContact::setReactionPtr(SP::SimpleVector newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("FrictionContact: setReactionPtr, inconsistent size between given reaction size and problem size. You should set sizeOutput before");
  reaction = newPtr;
}

void FrictionContact::setM(const OSNSMatrix& newValue)
{
  // Useless ?
  RuntimeException::selfThrow("FrictionContact::setM, forbidden operation. Try setMPtr().");
}

void FrictionContact::setMPtr(SP::OSNSMatrix newPtr)
{
  // Note we do not test if newPtr and M sizes are equal. Not necessary?
  M = newPtr;
}

void FrictionContact::setQ(const SimpleVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("FrictionContact: setQ, inconsistent size between given q size and problem size. You should set sizeOutput before");

  if (! q)
    q.reset(new SimpleVector(sizeOutput));

  *q = newValue;
}

void FrictionContact::setQPtr(SP::SimpleVector newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("FrictionContact: setQPtr, inconsistent size between given q size and problem size. You should set sizeOutput before");
  q = newPtr;

}

void FrictionContact::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);


  // Connect to the right function according to dim. of the problem
  if (contactProblemDim == 2)
    frictionContact_driver = &pfc_2D_driver;
  else // if(contactProblemDim == 3)
    frictionContact_driver = &frictionContact3D_driver;

  // Memory allocation for reaction, and velocity
  // If one of them has already been allocated, nothing is done.
  // We suppose that user has chosen a correct size.

  if (! velocity)
  {
    velocity.reset(new SimpleVector(maxSize));
  }
  else
  {
    if (velocity->size() != maxSize)
      velocity->resize(maxSize);
  }
  if (! reaction)
  {
    reaction.reset(new SimpleVector(maxSize));
  }
  else
  {
    if (reaction->size() != maxSize)
      reaction->resize(maxSize);
  }

  if (! q)
  {
    q.reset(new SimpleVector(maxSize));
  }
  else
  {
    if (q->size() != maxSize)
      q->resize(maxSize);
  }

  // get topology
  SP::Topology topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // Note that unitaryBlocks is up to date since updateUnitaryBlocks has been called during OneStepNSProblem::initialize()

  // Fill vector of friction coefficients
  SP::UnitaryRelationsSet I0 = topology->getIndexSet0Ptr();

  mu.reset(new std::vector<double>());

  mu->reserve(I0->size());

  // If the topology is TimeInvariant ie if M structure does not change during simulation:
  if (topology->isTimeInvariant() &&   !OSNSInteractions->isEmpty())
  {
    // Get index set from Simulation
    SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);
    if (! M)
    {
      // Creates and fills M using UR of indexSet

      M.reset(new OSNSMatrix(indexSet, unitaryBlocks, MStorageType));

    }
    else
    {
      M->setStorageType(MStorageType);
      M->fill(indexSet, unitaryBlocks);
    }

    for (ConstUnitaryRelationsIterator itUR = indexSet->begin(); itUR != indexSet->end(); ++itUR)
    {

      mu->push_back(boost::static_pointer_cast<NewtonImpactFrictionNSL> ((*itUR)->getInteractionPtr()->getNonSmoothLawPtr())->getMu());

    }
  }
  else // in that case, M and mu will be updated during preCompute
  {
    // Default size for M = maxSize
    if (! M)
    {

      if (MStorageType == 0)
        M.reset(new OSNSMatrix(maxSize, 0));
      else // if(MStorageType == 1) size = number of blocks = number of UR in the largest considered indexSet
        M.reset(new OSNSMatrix(simulation->getIndexSetPtr(levelMin)->size(), 1));
    }
  }
}

void FrictionContact::computeUnitaryBlock(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{

  // Computes matrix unitaryBlocks[UR1][UR2] (and allocates memory if necessary) if UR1 and UR2 have commond DynamicalSystem.
  // How unitaryBlocks are computed depends explicitely on the type of Relation of each UR.

  // Warning: we suppose that at this point, all non linear operators (G for lagrangian relation for example) have been computed through plug-in mechanism.
  // Get DS common between UR1 and UR2
  DynamicalSystemsSet commonDS;
  intersection(*UR1->getDynamicalSystemsPtr(), *UR2->getDynamicalSystemsPtr(), commonDS);

  if (!commonDS.isEmpty()) // Nothing to be done if there are no common DS between the two UR.
  {

    // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
    unsigned int nslawSize1 = UR1->getNonSmoothLawSize();
    unsigned int nslawSize2 = UR2->getNonSmoothLawSize();
    // Check allocation
    if (! unitaryBlocks[UR1][UR2])
    {

#ifndef WithSmartPtr
      unitaryBlocks[UR1][UR2] = new SimpleMatrix(nslawSize1, nslawSize2);
#else
      unitaryBlocks[UR1][UR2].reset(new SimpleMatrix(nslawSize1, nslawSize2));
#endif
    }

    // Get DS common between UR1 and UR2
    DSIterator itDS;

    // Get the W and Theta maps of one of the Unitary Relation - Warning: in the current version, if OSI!=Moreau, this fails.
    MapOfDSMatrices W;
    MapOfDouble Theta;
    getOSIMaps(UR1, W, Theta);

    SP::SiconosMatrix currentUnitaryBlock = unitaryBlocks[UR1][UR2];

    currentUnitaryBlock->zero();

    SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock;
    unsigned int sizeDS;
    RELATION::TYPES relationType1, relationType2;
    bool flagRightUnitaryBlock = false;

    // General form of the unitaryBlock is :   unitaryBlock = a*extraUnitaryBlock + b * leftUnitaryBlock * OP * rightUnitaryBlock
    // a and b are scalars, OP a matrix depending on the integrator, the simulation type ...
    // left, right and extra depend on the relation type and the non smooth law.

    // loop over the common DS
    for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
    {
      sizeDS = (*itDS)->getDim();

      // get unitaryBlocks corresponding to the current DS
      // These unitaryBlocks depends on the relation type.
      leftUnitaryBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));
      UR1->getLeftUnitaryBlockForDS(*itDS, leftUnitaryBlock);

      relationType1 = UR1->getRelationType();
      relationType2 = UR2->getRelationType();
      // Computing depends on relation type -> move this in UnitaryRelation method?
      if (relationType1 == Lagrangian || relationType2 == Lagrangian)
      {
        if (UR1 == UR2)
          rightUnitaryBlock = leftUnitaryBlock ;
        else
        {

          rightUnitaryBlock.reset(new SimpleMatrix(nslawSize2, sizeDS));
          UR2->getLeftUnitaryBlockForDS(*itDS, rightUnitaryBlock);
          // Warning: we use getLeft for Right unitaryBlock because right = transpose(left) and because of size checking inside the getBlock function,
          // a getRight call will fail.
          flagRightUnitaryBlock = true;
        }
        SP::SimpleMatrix work(new SimpleMatrix(*rightUnitaryBlock));
        work->trans();
        // W contains a lu-factorized matrix and we solve
        // W * X = rightUnitaryBlock with PLU
        // Work is a temporary matrix.
        W[*itDS]->PLUForwardBackwardInPlace(*work);
        //*currentUnitaryBlock +=  *leftUnitaryBlock * *work; // left = right = G or H
        gemm(CblasNoTrans, CblasNoTrans, 1.0, *leftUnitaryBlock, *work, 1.0, *currentUnitaryBlock);
      }
      else RuntimeException::selfThrow("FrictionContact::computeUnitaryBlock not yet implemented for relation of type " + relationType1);
    }
  }
}

void FrictionContact::computeQBlock(SP::UnitaryRelation UR, unsigned int pos)
{

  // Get relation and non smooth law types
  RELATION::TYPES relationType = UR->getRelationType();
  RELATION::SUBTYPES relationSubType = UR->getRelationSubType();
  string nslawType = UR->getNonSmoothLawType();

  string simulationType = simulation->getType();

  SP::DynamicalSystem ds = *(UR->dynamicalSystemsBegin());
  string osiType = simulation->getIntegratorOfDSPtr(ds)->getType();

  unsigned int sizeY = UR->getNonSmoothLawSize();
  std::vector<unsigned int> coord(8);

  unsigned int relativePosition = UR->getRelativePosition();
  SP::Interaction mainInteraction = UR->getInteractionPtr();
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = pos;
  coord[7] = pos + sizeY;

  SP::SiconosMatrix  H;
  SP::SiconosVector workX = UR->getWorkXPtr();
  if (osiType == "Moreau" || osiType == "Lsodar")
  {
    if (relationType == FirstOrder)
    {
      if (relationSubType == Type1R)//|| relationType =="FirstOrderType2R" || relationType =="FirstOrderType3R")
      {
        H = boost::static_pointer_cast<FirstOrderR>(mainInteraction->getRelationPtr())->getJacobianHPtr(0);
        if (H)
        {
          coord[3] = H->size(1);
          coord[5] = H->size(1);
          subprod(*H, *workX, *q, coord, true);
        }
      }

      else if (relationSubType == LinearTIR || relationSubType == LinearR)
      {
        // q = HXfree + e + Fz
        H = boost::static_pointer_cast<FirstOrderLinearR>(mainInteraction->getRelationPtr())->getCPtr();
        if (H)
        {
          coord[3] = H->size(1);
          coord[5] = H->size(1);
          subprod(*H, *workX, (*q), coord, true);
        }
        SP::SiconosVector  e = boost::static_pointer_cast<FirstOrderLinearR>(mainInteraction->getRelationPtr())->getEPtr();
        if (e)
          boost::static_pointer_cast<SimpleVector>(q)->addBlock(pos, *e);

        H = boost::static_pointer_cast<FirstOrderLinearR>(mainInteraction->getRelationPtr())->getFPtr();
        if (H)
        {
          SP::SiconosVector  workZ = UR->getWorkZPtr();
          coord[3] = H->size(1);
          coord[5] = H->size(1);
          subprod(*H, *workZ, *q, coord, false);
        }
      }
    }
    else if (relationType == Lagrangian)
    {
      if (relationSubType == CompliantR || relationSubType == ScleronomousR || relationSubType == RheonomousR)
      {
        // q = jacobian_q h().v_free
        H = boost::static_pointer_cast<LagrangianR>(mainInteraction->getRelationPtr())->getGPtr(0);
        if (H)
        {
          coord[3] = H->size(1);
          coord[5] = H->size(1);
          subprod(*H, *workX, *q, coord, true);
        }
      }
      else if (relationSubType == LinearR || relationSubType == LinearTIR)
      {
        // q = H.v_free
        H = boost::static_pointer_cast<LagrangianLinearR>(mainInteraction->getRelationPtr())->getHPtr();
        if (H)
        {
          coord[3] = H->size(1);
          coord[5] = H->size(1);
          subprod(*H, *workX, *q, coord, true);
        }
      }
    }
    else
      RuntimeException::selfThrow("FrictionContact::computeQBlock, not yet implemented for first order relations of subtype " + relationType);

  }
  else if (osiType == "Moreau2")
  {
  }
  else
    RuntimeException::selfThrow("FrictionContact::computeQBlock not yet implemented for OSI of type " + osiType);

  // Add "non-smooth law effect" on q
  if (UR->getRelationType() == Lagrangian)
  {
    double e;
    if (nslawType == NEWTONIMPACTNSLAW)
    {
      e = (boost::static_pointer_cast<NewtonImpactNSL>(mainInteraction->getNonSmoothLawPtr()))->getE();

      std::vector<unsigned int> subCoord(4);
      if (simulationType == "TimeStepping")
      {
        subCoord[0] = 0;
        subCoord[1] = UR->getNonSmoothLawSize();
        subCoord[2] = pos;
        subCoord[3] = pos + subCoord[1];
        subscal(e, *UR->getYOldPtr(levelMin), *q, subCoord, false);
      }
      else if (simulationType == "EventDriven")
      {
        subCoord[0] = pos;
        subCoord[1] = pos + UR->getNonSmoothLawSize();
        subCoord[2] = pos;
        subCoord[3] = subCoord[1];
        subscal(e, *q, *q, subCoord, false); // q = q + e * q
      }
      else
        RuntimeException::selfThrow("FrictionContact::computeQBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {

      e = (boost::static_pointer_cast<NewtonImpactFrictionNSL>(mainInteraction->getNonSmoothLawPtr()))->getEn();

      // Only the normal part is multiplied by e
      if (simulationType == "TimeStepping")
        (*q)(pos) +=  e * (*UR->getYOldPtr(levelMin))(0);

      else RuntimeException::selfThrow("FrictionContact::computeQBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);

    }
    else
      RuntimeException::selfThrow("FrictionContact::computeQBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType);
  }
}


void FrictionContact::computeQ(const double time)
{
  if (q->size() != sizeOutput)
    q->resize(sizeOutput);
  q->zero();

  // === Get index set from Simulation ===
  SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);

  // === Loop through "active" Unitary Relations (ie present in indexSets[level]) ===

  unsigned int pos = 0;
  UnitaryRelationsIterator itCurrent, itLiked;
  string simulationType = simulation->getType();
  for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // *itCurrent is a SP::UnitaryRelation.

    // Compute q, this depends on the type of non smooth problem, on the relation type and on the non smooth law
    pos = M->getPositionOfUnitaryBlock(*itCurrent);
    computeQBlock((*itCurrent), pos); // free output is saved in y
  }
}

void FrictionContact::preCompute(const double time)
{
  // This function is used to prepare data for the FrictionContact problem
  // - computation of M and q
  // - set sizeOutput
  // - check dim. for z,w

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, sizeOutput have been computed in initialize and are uptodate.

  // Get topology

  assert(M || !"M is NULL");


  SP::Topology topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  if (!topology->isTimeInvariant())
  {
    // Computes new unitaryBlocks if required
    updateUnitaryBlocks();

    // Updates matrix M
    SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);
    M->fill(indexSet, unitaryBlocks);
    sizeOutput = M->size();

    // Checks z and w sizes and reset if necessary
    if (reaction->size() != sizeOutput)
    {
      reaction->resize(sizeOutput, false);
      reaction->zero();
    }

    if (velocity->size() != sizeOutput)
    {
      velocity->resize(sizeOutput);
      velocity->zero();
    }

    // Update mu
    mu->clear();
    for (ConstUnitaryRelationsIterator itUR = indexSet->begin(); itUR != indexSet->end(); ++itUR)
    {

      mu->push_back(boost::static_pointer_cast<NewtonImpactFrictionNSL>((*itUR)->getInteractionPtr()->getNonSmoothLawPtr())->getMu());
    }

  }

  // Computes q
  computeQ(time);
}

int FrictionContact::compute(double time)
{
  int info = 0;
  // --- Prepare data for FrictionContact computing ---
  preCompute(time);

  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)
  if (sizeOutput != 0)
  {
    // The FrictionContact Problem in Numerics format
    FrictionContact_Problem numerics_problem;
    numerics_problem.M = &*M->getNumericsMatrix();
    numerics_problem.q = &*q->getArray();
    numerics_problem.numberOfContacts = sizeOutput / contactProblemDim;
    numerics_problem.isComplete = 1;
    numerics_problem.mu = &((*mu)[0]);
    // Call Numerics Driver for FrictionContact
    info = (*frictionContact_driver)(&numerics_problem, &*reaction->getArray() , &*velocity->getArray() , &*solver->getNumericsSolverOptionsPtr(), &*numerics_options);
    //  cout << " step 2 "<< info << endl;
    //  reaction->display();
    //  velocity->display();
    // --- Recovering of the desired variables from LCP output ---
    postCompute();

  }

  return info;
}

void FrictionContact::postCompute()
{
  // This function is used to set y/lambda values using output from lcp_driver (w,z).
  // Only UnitaryRelations (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);

  // y and lambda vectors
  SP::SiconosVector y;
  SP::SiconosVector lambda;

  //   // === Loop through "active" Unitary Relations (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;

  for (UnitaryRelationsIterator itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // size of the unitaryBlock that corresponds to the current UnitaryRelation
    nsLawSize = (*itCurrent)->getNonSmoothLawSize();
    // Get the relative position of UR-unitaryBlock in the vector velocity or reaction
    pos = M->getPositionOfUnitaryBlock(*itCurrent);

    // Get Y and Lambda for the current Unitary Relation
    y = (*itCurrent)->getYPtr(levelMin);
    lambda = (*itCurrent)->getLambdaPtr(levelMin);

    // Copy velocity/reaction values, starting from index pos into y/lambda.

    setBlock(velocity, y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
    setBlock(reaction, lambda, lambda->size(), pos, 0);

  }
}

void FrictionContact::display() const
{
  cout << "===== " << contactProblemDim << "D Friction Contact Problem " << endl;
  cout << "of size " << sizeOutput << "(ie " << sizeOutput / contactProblemDim << " contacts)." << endl;
  cout << " - Matrix M  : " << endl;
  if (M) M->display();
  else cout << "-> NULL" << endl;
  cout << " - Vector q : " << endl;
  if (q) q->display();
  else cout << "-> NULL" << endl;
  cout << " Friction coefficients: " << endl;
  if (mu) print(mu->begin(), mu->end());
  else cout << "-> NULL" << endl;
  cout << "============================================================" << endl;
}

void FrictionContact::saveNSProblemToXML()
{
  //   OneStepNSProblem::saveNSProblemToXML();
  //   if(onestepnspbxml != NULL)
  //     {
  //       (boost::static_pointer_cast<FrictionContactXML>(onestepnspbxml))->setM(*M);
  //       (boost::static_pointer_cast<FrictionContactXML>(onestepnspbxml))->setQ(*q);
  //     }
  //   else RuntimeException::selfThrow("FrictionContact::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  RuntimeException::selfThrow("FrictionContact::saveNSProblemToXML - Not yet implemented.");
}

FrictionContact* FrictionContact::convert(OneStepNSProblem* osnsp)
{
  FrictionContact* fc2d = dynamic_cast<FrictionContact*>(osnsp);
  return fc2d;
}



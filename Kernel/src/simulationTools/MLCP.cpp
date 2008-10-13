/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "MLCP.h"
#include "Topology.h"
#include "UnitaryRelation.h"
#include "MixedComplementarityConditionNSL.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "Relation.h"
#include "DynamicalSystem.h"
#include "TimeDiscretisation.h"
#include "MixedLinearComplementarity_Problem.h" // Numerics structure
#include "NumericsMatrix.h"
#include "RelationTypes.h"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"

#include "OneStepIntegrator.h"
//#include "mlcpDefaultSolver.h"
#include <stdio.h>
#include <stdlib.h>

using namespace std;
using namespace RELATION;

// xml constructor
MLCP::MLCP(SP::OneStepNSProblemXML onestepnspbxml):
  OneStepNSProblem("MLCP", onestepnspbxml), MStorageType(0)
{}

// Constructor from a set of data
MLCP::MLCP(SP::NonSmoothSolver newSolver, const string& newId):
  OneStepNSProblem("MLCP", newId, newSolver), MStorageType(0)
{}

// Setters

void MLCP::setW(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("MLCP: setW, inconsistent size between given w size and problem size. You should set sizeOutput before");

  if (! w)
    w.reset(new SimpleVector(sizeOutput));

  else if (w->size() != sizeOutput)
    RuntimeException::selfThrow("MLCP: setW, w size differs from sizeOutput");

  *w = newValue;
}

void MLCP::setWPtr(SP::SiconosVector newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("MLCP: setWPtr, inconsistent size between given w size and problem size. You should set sizeOutput before");
  w = newPtr;

}


void MLCP::setZ(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("MLCP: setZ, inconsistent size between given z size and problem size. You should set sizeOutput before");

  if (! z)
  {
#ifndef WithSmartPtr
    z = new SimpleVector(sizeOutput);
    isZAllocatedIn = true;
#else
    z.reset(new SimpleVector(sizeOutput));
#endif
  }

  *z = newValue;
}

void MLCP::setZPtr(SP::SiconosVector newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("MLCP: setZPtr, inconsistent size between given z size and problem size. You should set sizeOutput before");

#ifndef WithSmartPtr
  if (isZAllocatedIn) delete z;
  isZAllocatedIn = false;
#endif

  z = newPtr;

}

void MLCP::setM(const OSNSMatrix& newValue)
{
  // Useless ?
  RuntimeException::selfThrow("MLCP: setM, forbidden operation. Try setMPtr().");
}

void MLCP::setMPtr(SP::OSNSMatrix newPtr)
{

  // Note we do not test if newPtr and M sizes are equal. Not necessary?
#ifndef WithSmartPtr
  if (isMAllocatedIn)
    delete M;
  isMAllocatedIn = false;
#endif

  M = newPtr;

}

void MLCP::setQ(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("MLCP: setQ, inconsistent size between given q size and problem size. You should set sizeOutput before");

  if (! q)
  {

#ifndef WithSmartPtr
    q = new SimpleVector(sizeOutput);
    isQAllocatedIn = true;
#else
    q.reset(new SimpleVector(sizeOutput));
#endif

  }

  *q = newValue;
}

void MLCP::setQPtr(SP::SiconosVector newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("MLCP: setQPtr, inconsistent size between given q size and problem size. You should set sizeOutput before");

#ifndef WithSmartPtr
  if (isQAllocatedIn) delete q;
  isQAllocatedIn = false;
#endif

  q = newPtr;


}

void MLCP::initialize(SP::Simulation sim)
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize(sim);

  // Memory allocation for w, M, z and q.
  // If one of them has already been allocated, nothing is done.
  // We suppose that user has chosen a correct size.

  if (!w)
  {
    w.reset(new SimpleVector(maxSize));
  }
  else
  {
    if (w->size() != maxSize)
      w->resize(maxSize);
  }
  if (! z)
  {
    z.reset(new SimpleVector(maxSize));
  }
  else
  {
    if (z->size() != maxSize)
      z->resize(maxSize);
  }

  if (! q)
  {
    q.reset(new SimpleVector(maxSize));
  }

  // get topology
  SP::Topology topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // Note that unitaryBlocks is up to date since updateUnitaryBlocks has been called during OneStepNSProblem::initialize()

  // If the topology is TimeInvariant ie if M structure does not change during simulation:
  if (topology->isTimeInvariant() &&   !OSNSInteractions->isEmpty())
  {
    updateM();
  }
  else // in that case, M will be updated during preCompute
  {
    // Default size for M = maxSize
    if (!M)
    {
      if (MStorageType == 0)
        M.reset(new OSNSMatrix(maxSize, 0));
      else // if(MStorageType == 1) size = number of unitaryBlocks = number of UR in the largest considered indexSet
        M.reset(new OSNSMatrix(simulation->getIndexSetPtr(levelMin)->size(), 1));
    }
  }
}

void MLCP::updateM()
{
  // Get index set from Simulation
  SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);

  if (!M)
  {
    // Creates and fills M using UR of indexSet
    M.reset(new OSNSMatrix(indexSet, unitaryBlocks, MStorageType));
    numerics_problem.M = &*M->getNumericsMatrix();
    numerics_problem.A = 0;
    numerics_problem.B = 0;
    numerics_problem.C = 0;
    numerics_problem.D = 0;
    numerics_problem.a = 0;
    numerics_problem.b = 0;
    numerics_problem.problemType = 0;
    numerics_problem.n = n;
    numerics_problem.m = m;
  }
  else
  {
    M->setStorageType(MStorageType);
    M->fill(indexSet, unitaryBlocks);

  }
  sizeOutput = M->size();
}
void  MLCP::reset()
{
  mlcp_driver_reset(&numerics_problem, solver->getNumericsSolverOptionsPtr());
}

void MLCP::computeUnitaryBlock(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{

  // Computes matrix unitaryBlocks[UR1][UR2] (and allocates memory if necessary) if UR1 and UR2 have commond DynamicalSystem.
  // How unitaryBlocks are computed depends explicitely on the type of Relation of each UR.

  // Get DS common between UR1 and UR2
  DynamicalSystemsSet commonDS;
  intersection(*UR1->getDynamicalSystemsPtr(), *UR2->getDynamicalSystemsPtr(), commonDS);
  m = 0;
  n = 0;
  if (!commonDS.isEmpty()) // Nothing to be done if there are no common DS between the two UR.
  {
    DSIterator itDS;
    // Warning: we suppose that at this point, all non linear operators (G for lagrangian relation for example) have been computed through plug-in mechanism.

    // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
    unsigned int nslawSize1 = UR1->getNonSmoothLawSize();
    unsigned int nslawSize2 = UR2->getNonSmoothLawSize();
    unsigned int equalitySize1 =  MixedComplementarityConditionNSL::convert(UR1->getInteractionPtr()->getNonSmoothLawPtr())->getEqualitySize();
    unsigned int equalitySize2 = MixedComplementarityConditionNSL::convert(UR1->getInteractionPtr()->getNonSmoothLawPtr())->getEqualitySize();



    // Check allocation
    if (! unitaryBlocks[UR1][UR2])
    {
#ifndef WithSmartPtr
      unitaryBlocks[UR1][UR2] = new SimpleMatrix(nslawSize1, nslawSize2);
#else
      (unitaryBlocks[UR1][UR2]).reset(new SimpleMatrix(nslawSize1, nslawSize2));
#endif

    }
    // Get the W and Theta maps of one of the Unitary Relation - Warning: in the current version, if OSI!=Moreau, this fails.
    // If OSI = MOREAU, centralUnitaryBlocks = W
    // if OSI = LSODAR, centralUnitaryBlocks = M (mass matrices)
    MapOfDSMatrices centralUnitaryBlocks;
    MapOfDouble Theta; // If OSI = LSODAR, Theta remains empty
    getOSIMaps(UR1, centralUnitaryBlocks, Theta);

    SP::SiconosMatrix currentUnitaryBlock = unitaryBlocks[UR1][UR2];

    SP::SiconosMatrix leftUnitaryBlock, rightUnitaryBlock;

    unsigned int sizeDS;
    RELATION::TYPES relationType1, relationType2;

    double h = simulation->getTimeDiscretisationPtr()->getCurrentTimeStep();
    printf("h : %f \n", h);

    // General form of the unitaryBlock is :   unitaryBlock = a*extraUnitaryBlock + b * leftUnitaryBlock * centralUnitaryBlocks * rightUnitaryBlock
    // a and b are scalars, centralUnitaryBlocks a matrix depending on the integrator (and on the DS), the simulation type ...
    // left, right and extra depend on the relation type and the non smooth law.
    relationType1 = UR1->getRelationType();
    relationType2 = UR2->getRelationType();
    // ==== First Order Relations - Specific treatment for diagonal unitaryBlocks ===
    if (UR1 == UR2)
    {
      UR1->getExtraUnitaryBlock(currentUnitaryBlock);
      m += nslawSize1 - equalitySize1;
      n += equalitySize1;
    }
    else
      currentUnitaryBlock->zero();


    // loop over the common DS
    for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
    {
      sizeDS = (*itDS)->getDim();
      // get blocks corresponding to the current DS
      // These blocks depends on the relation type.
      leftUnitaryBlock.reset(new SimpleMatrix(nslawSize1, sizeDS));

      UR1->getLeftUnitaryBlockForDS(*itDS, leftUnitaryBlock);
      // Computing depends on relation type -> move this in UnitaryRelation method?
      if (relationType1 == FirstOrder && relationType2 == FirstOrder)
      {
        rightUnitaryBlock.reset(new SimpleMatrix(sizeDS, nslawSize2));

        UR2->getRightUnitaryBlockForDS(*itDS, rightUnitaryBlock);
        // centralUnitaryBlock contains a lu-factorized matrix and we solve
        // centralUnitaryBlock * X = rightUnitaryBlock with PLU
        //        printf("right bloc: ie B \n");
        //        rightUnitaryBlock->display();
        centralUnitaryBlocks[*itDS]->PLUForwardBackwardInPlace(*rightUnitaryBlock);
        //        printf("W \n");
        //        centralUnitaryBlocks[*itDS]->display();
        //        printf("inv(W)B \n");
        //        rightUnitaryBlock->display();
        //      integration of r with theta method removed
        //      *currentUnitaryBlock += h *Theta[*itDS]* *leftUnitaryBlock * (*rightUnitaryBlock); //left = C, right = W.B
        //gemm(h,*leftUnitaryBlock,*rightUnitaryBlock,1.0,*currentUnitaryBlock);
        *leftUnitaryBlock *= h;
        //        printf("currentbloc : ie D \n");
        //        currentUnitaryBlock->display();
        //        printf("leftUnitaryBlock : ie C \n");
        //        leftUnitaryBlock->display();

        prod(*leftUnitaryBlock, *rightUnitaryBlock, *currentUnitaryBlock, false);
        //left = C, right = W.B
        //        printf("finalbloc \n");
        //        currentBlock->display();

#ifndef WithSmartPtr
        delete rightUnitaryBlock;
#endif
      }
      else RuntimeException::selfThrow("MLCP::computeBlock not yet implemented for relation of type " + relationType1);

#ifndef WithSmartPtr
      delete leftUnitaryBlock;
#endif
    }

  }
}

void MLCP::computeQBlock(SP::UnitaryRelation UR, unsigned int pos)
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

  SP::SiconosMatrix H;
  SP::SiconosVector workX = UR->getWorkXPtr();
  if (osiType == "Moreau" || osiType == "Lsodar")
  {
    if (relationType == FirstOrder)
    {
      if (relationSubType == Type1R) // || relationSubType =="FirstOrderType2R" || relationType =="FirstOrderType3R")
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
      RuntimeException::selfThrow("MLCP::getExtraUnitaryBlock, not yet implemented for relations of subtype " + relationSubType);

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

#ifndef WithSmartPtr
      e = (boost::static_pointer_cast<NewtonImpactNSL>(mainInteraction->getNonSmoothLawPtr()))->getE();
#else
      e = (boost::static_pointer_cast<NewtonImpactNSL>(mainInteraction->getNonSmoothLawPtr()))->getE();
#endif

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
        RuntimeException::selfThrow("MLCP::computeQBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {

#ifndef WithSmartPtr
      e = (boost::static_pointer_cast<NewtonImpactFrictionNSL>(mainInteraction->getNonSmoothLawPtr()))->getEn();
#else
      e = (boost::static_pointer_cast<NewtonImpactFrictionNSL>(mainInteraction->getNonSmoothLawPtr()))->getEn();
#endif

      // Only the normal part is multiplied by e
      if (simulationType == "TimeStepping")
        (*q)(pos) +=  e * (*UR->getYOldPtr(levelMin))(0);

      else RuntimeException::selfThrow("MLCP::computeQBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);

    }
    else
      RuntimeException::selfThrow("MLCP::computeQBlock not yet implemented for this type of relation and a non smooth law of type " + nslawType);
  }
}

void MLCP::computeQ(double time)
{
  if (q->size() != sizeOutput)
    q->resize(sizeOutput);
  q->zero();

  // === Get index set from Simulation ===
  SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);
  // === Loop through "active" Unitary Relations (ie present in indexSets[level]) ===

  unsigned int pos = 0;
  UnitaryRelationsIterator itCurrent, itLinked;
  string simulationType = simulation->getType();
  for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // *itCurrent is a SP::UnitaryRelation.

    // Compute q, this depends on the type of non smooth problem, on the relation type and on the non smooth law
    pos = M->getPositionOfUnitaryBlock(*itCurrent);
    computeQBlock((*itCurrent), pos); // free output is saved in y

  }
}

void MLCP::preCompute(double time)
{
  // This function is used to prepare data for the MixedLinearComplementarity_Problem
  // - computation of M and q
  // - set sizeOutput
  // - check dim. for z,w

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, sizeOutput have been computed in initialize and are uptodate.

  // Get topology
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
    if (z->size() != sizeOutput)
    {
      z->resize(sizeOutput, false);
      z->zero();
    }

    if (w->size() != sizeOutput)
    {
      w->resize(sizeOutput);
      w->zero();
    }
  }

  // Computes q of MLCP
  computeQ(time);

}
void displayNM(const NumericsMatrix* const m)
{
  if (!m)
  {
    fprintf(stderr, "Numerics, NumericsMatrix display failed, NULL input.\n");
    exit(EXIT_FAILURE);
  }
  int storageType = m->storageType;
  if (storageType == 0)
  {
    printf("\n ========== Numerics Matrix of dim %dX%d\n", m->size0, m->size1);
    printf("[");
    for (int i = 0; i < m->size1 * m->size0; i++)
    {
      printf("%lf ", m->matrix0[i]);
      if ((i + 1) % m->size1 == 0)
        printf("\n");
    }
    printf("]");
    printf("\n (warning: column-major) \n");
  }
  else if (storageType == 1)
    fprintf(stderr, "storageType NumericsdisplayNM.\n");

}
int MLCP::compute(double time)
{
  // --- Prepare data for MLCP computing ---
  preCompute(time);

  int info = 0;
  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (sizeOutput != 0)
  {
    numerics_problem.q = q->getArray();
    int nbSolvers = 1;
    // Call MLCP Driver
    //printf("MLCP display");
    //printf("n %d m %d",n,m);
    //displayNM(numerics_problem.M);
    //      exit(1);
    //mlcpDefaultSolver *pSolver = new mlcpDefaultSolver(m,n);
    //      displayMLCP(&numerics_problem);
    info = mlcp_driver(&numerics_problem, &*z->getArray(), &*w->getArray(), &*solver->getNumericsSolverOptionsPtr(), &*numerics_options);

    // --- Recovering of the desired variables from MLCP output ---
    postCompute();

  }

  return info;
}

void MLCP::postCompute()
{
  // This function is used to set y/lambda values using output from lcp_driver (w,z).
  // Only UnitaryRelations (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::UnitaryRelationsSet indexSet = simulation->getIndexSetPtr(levelMin);

  // y and lambda vectors
  SP::SiconosVector lambda;
  SP::SiconosVector y;

  // === Loop through "active" Unitary Relations (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;

  for (UnitaryRelationsIterator itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // size of the unitaryBlock that corresponds to the current UnitaryRelation
    nsLawSize = (*itCurrent)->getNonSmoothLawSize();
    // Get the relative position of UR-unitaryBlock in the vector w or z
    pos = M->getPositionOfUnitaryBlock(*itCurrent);

    // Get Y and Lambda for the current Unitary Relation
    y = (*itCurrent)->getYPtr(levelMin);
    lambda = (*itCurrent)->getLambdaPtr(levelMin);
    // Copy w/z values, starting from index pos into y/lambda.
    setBlock(w, y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
    setBlock(z, lambda, lambda->size(), pos, 0);
  }
}

void MLCP::display() const
{
  cout << "======= MLCP of size " << sizeOutput << " with: " << endl;
  cout << "======= m " << m << " n " << n << endl;
  cout << "M  ";
  if (M) M->display();
  else cout << "-> NULL" << endl;
  cout << endl << " q : " ;
  if (q) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;
}

void MLCP::saveNSProblemToXML()
{
  OneStepNSProblem::saveNSProblemToXML();
  //   if(onestepnspbxml != NULL)
  //     {
  // //       (static_cast<MLCPXML>(onestepnspbxml))->setM(*M);
  //       (static_cast<MLCPXML>(onestepnspbxml))->setQ(*q);
  //     }
  //   else RuntimeException::selfThrow("MLCP::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  RuntimeException::selfThrow("MLCP::saveNSProblemToXML - Not yet implemented.");
}

MLCP* MLCP::convert(OneStepNSProblem* osnsp)
{
  MLCP* lcp = dynamic_cast<MLCP*>(osnsp);
  return lcp;
}



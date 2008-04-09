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
#include "LCP.h"
#include "Topology.h"
#include "UnitaryRelation.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "Relation.h"
#include "DynamicalSystem.h"
#include "TimeDiscretisation.h"
#include "LinearComplementarity_Problem.h" // Numerics structure

using namespace std;

// xml constructor
LCP::LCP(OneStepNSProblemXML* onestepnspbxml, Simulation* newSimu):
  OneStepNSProblem("LCP", onestepnspbxml, newSimu),

#ifndef WithSmartPtr
  w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false),
#endif

  MStorageType(0)
{}

// Constructor from a set of data
LCP::LCP(Simulation* newSimu, NonSmoothSolver* newSolver, const string& newId):
  OneStepNSProblem("LCP", newSimu, newId, newSolver),

#ifndef WithSmartPtr
  w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false),
#endif

  MStorageType(0)
{}

// destructor
LCP::~LCP()
{

#ifndef WithSmartPtr
  if (isWAllocatedIn) delete w;
  w = NULL;
  if (isZAllocatedIn) delete z;
  z = NULL;
  if (isMAllocatedIn)
    delete M;
  M = NULL;
  if (isQAllocatedIn) delete q;
  q = NULL;
#endif

}

// Setters

void LCP::setW(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("LCP: setW, inconsistent size between given w size and problem size. You should set sizeOutput before");

  if (! w)
  {
#ifndef WithSmartPtr
    w = new SimpleVector(sizeOutput);
    isWAllocatedIn = true;
#else
    w.reset(new SimpleVector(sizeOutput));
#endif
  }
  else if (w->size() != sizeOutput)
    RuntimeException::selfThrow("LCP: setW, w size differs from sizeOutput");

  *w = newValue;
}

void LCP::setWPtr(SiconosVectorSPtr newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("LCP: setWPtr, inconsistent size between given w size and problem size. You should set sizeOutput before");

#ifndef WithSmartPtr
  if (isWAllocatedIn) delete w;
  isWAllocatedIn = false;
#endif

  w = newPtr;

}


void LCP::setZ(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("LCP: setZ, inconsistent size between given z size and problem size. You should set sizeOutput before");

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

void LCP::setZPtr(SiconosVectorSPtr newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("LCP: setZPtr, inconsistent size between given z size and problem size. You should set sizeOutput before");

#ifndef WithSmartPtr
  if (isZAllocatedIn) delete z;
  isZAllocatedIn = false;
#endif

  z = newPtr;

}

void LCP::setM(const OSNSMatrix& newValue)
{
  // Useless ?
  RuntimeException::selfThrow("LCP: setM, forbidden operation. Try setMPtr().");
}

void LCP::setMPtr(OSNSMatrixSPtr newPtr)
{

#ifndef WithSmartPtr
  // Note we do not test if newPtr and M sizes are equal. Not necessary?
  if (isMAllocatedIn)
    delete M;
  isMAllocatedIn = false;
#endif

  M = newPtr;

}

void LCP::setQ(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("LCP: setQ, inconsistent size between given q size and problem size. You should set sizeOutput before");

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

void LCP::setQPtr(SiconosVectorSPtr newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("LCP: setQPtr, inconsistent size between given q size and problem size. You should set sizeOutput before");

#ifndef WithSmartPtr
  if (isQAllocatedIn) delete q;
  isQAllocatedIn = false;
#endif

  q = newPtr;


}

void LCP::initialize()
{
  // - Checks memory allocation for main variables (M,q,w,z)
  // - Formalizes the problem if the topology is time-invariant

  // This function performs all steps that are time-invariant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize();

  // Memory allocation for w, M, z and q.
  // If one of them has already been allocated, nothing is done.
  // We suppose that user has chosen a correct size.
  if (! w)
  {

#ifndef WithSmartPtr
    w = new SimpleVector(maxSize);
    isWAllocatedIn = true;
#else
    w.reset(new SimpleVector(maxSize));
#endif

  }
  else
  {
    if (w->size() != maxSize)
      w->resize(maxSize);
  }
  if (! z)
  {

#ifndef WithSmartPtr
    z = new SimpleVector(maxSize);
    isZAllocatedIn = true;
#else
    z.reset(new SimpleVector(maxSize));
#endif

  }
  else
  {
    if (z->size() != maxSize)
      z->resize(maxSize);
  }

  if (! q)
  {

#ifndef WithSmartPtr
    q = new SimpleVector(maxSize);
    isQAllocatedIn = true;
#else
    q.reset(new SimpleVector(maxSize));
#endif

  }

  // get topology
  Topology * topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // Note that blocks is up to date since updateBlocks has been called during OneStepNSProblem::initialize()

  // If the topology is TimeInvariant ie if M structure does not change during simulation:
  if (topology->isTimeInvariant() &&   !OSNSInteractions->isEmpty())
  {
    // Get index set from Simulation
    UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);
    if (! M)
    {
      // Creates and fills M using UR of indexSet

#ifndef WithSmartPtr
      M = new OSNSMatrix(indexSet, blocks, MStorageType);
      isMAllocatedIn = true;
#else
      M.reset(new OSNSMatrix(indexSet, blocks, MStorageType));
#endif
    }
    else
    {
      M->setStorageType(MStorageType);
      M->fill(indexSet, blocks);
    }
    sizeOutput = M->size();
  }
  else // in that case, M will be updated during preCompute
  {
    // Default size for M = maxSize
    if (! M)
    {
      if (MStorageType == 0)
      {

#ifndef WithSmartPtr
        M = new OSNSMatrix(maxSize, 0);
#else
        M.reset(new OSNSMatrix(maxSize, 0));
#endif
      }

      else   // if(MStorageType == 1) size = number of blocks = number of UR in the largest considered indexSet
      {

#ifndef WithSmartPtr
        M = new OSNSMatrix(simulation->getIndexSetPtr(levelMin)->size(), 1);
        isMAllocatedIn = true;
#else
        M.reset(new OSNSMatrix(maxSize, 0));
#endif
      }
    }
  }
}

void LCP::computeBlock(UnitaryRelation* UR1, UnitaryRelation* UR2)
{

  // Computes matrix blocks[UR1][UR2] (and allocates memory if necessary) if UR1 and UR2 have commond DynamicalSystem.
  // How blocks are computed depends explicitely on the type of Relation of each UR.

  // Get DS common between UR1 and UR2
  DynamicalSystemsSet commonDS;
  intersection(*UR1->getDynamicalSystemsPtr(), *UR2->getDynamicalSystemsPtr(), commonDS);

  if (!commonDS.isEmpty()) // Nothing to be done if there are no common DS between the two UR.
  {
    DSIterator itDS;
    // Warning: we suppose that at this point, all non linear operators (G for lagrangian relation for example) have been computed through plug-in mechanism.

    // Get dimension of the NonSmoothLaw (ie dim of the block)
    unsigned int nslawSize1 = UR1->getNonSmoothLawSize();
    unsigned int nslawSize2 = UR2->getNonSmoothLawSize();
    // Check allocation
    if (blocks[UR1][UR2] == NULL)
      blocks[UR1][UR2] = new SimpleMatrix(nslawSize1, nslawSize2);

    // Get the W and Theta maps of one of the Unitary Relation - Warning: in the current version, if OSI!=Moreau, this fails.
    // If OSI = MOREAU, centralBlocks = W
    // if OSI = LSODAR, centralBlocks = M (mass matrices)
    MapOfMatrices centralBlocks;
    MapOfDouble Theta; // If OSI = LSODAR, Theta remains empty
    getOSIMaps(UR1, centralBlocks, Theta);

    SiconosMatrix* currentBlock = blocks[UR1][UR2];
    SiconosMatrix *leftBlock = NULL, *rightBlock = NULL;
    unsigned int sizeDS;
    string relationType1, relationType2;
    double h = simulation->getTimeDiscretisationPtr()->getH();

    // General form of the block is :   block = a*extraBlock + b * leftBlock * centralBlocks * rightBlock
    // a and b are scalars, centralBlocks a matrix depending on the integrator (and on the DS), the simulation type ...
    // left, right and extra depend on the relation type and the non smooth law.
    relationType1 = UR1->getRelationType();
    relationType2 = UR2->getRelationType();
    // ==== First Order Relations - Specific treatment for diagonal blocks ===
    if (UR1 == UR2)
      UR1->getExtraBlock(currentBlock);
    else
      currentBlock->zero();


    // loop over the common DS
    for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
    {
      sizeDS = (*itDS)->getDim();
      // get blocks corresponding to the current DS
      // These blocks depends on the relation type.
      leftBlock = new SimpleMatrix(nslawSize1, sizeDS);

      UR1->getLeftBlockForDS(*itDS, leftBlock);
      // Computing depends on relation type -> move this in UnitaryRelation method?
      if (relationType1 == "FirstOrder" && relationType2 == "FirstOrder")
      {
        rightBlock = new SimpleMatrix(sizeDS, nslawSize2);
        UR2->getRightBlockForDS(*itDS, rightBlock);

        // centralBlock contains a lu-factorized matrix and we solve
        // centralBlock * X = rightBlock with PLU
        centralBlocks[*itDS]->PLUForwardBackwardInPlace(*rightBlock);
        //      integration of r with theta method removed
        //      *currentBlock += h *Theta[*itDS]* *leftBlock * (*rightBlock); //left = C, right = W.B
        //gemm(h,*leftBlock,*rightBlock,1.0,*currentBlock);
        *leftBlock *= h;
        prod(*leftBlock, *rightBlock, *currentBlock, false);
        //left = C, right = W.B
        delete rightBlock;
      }
      else if (relationType1 == "Lagrangian" || relationType2 == "Lagrangian")
      {
        if (UR1 == UR2)
        {
          SimpleMatrix * work = new SimpleMatrix(*leftBlock);
          work->trans();
          centralBlocks[*itDS]->PLUForwardBackwardInPlace(*work);
          //*currentBlock +=  *leftBlock ** work;
          prod(*leftBlock, *work, *currentBlock, false);
          //      gemm(CblasNoTrans,CblasNoTrans,1.0,*leftBlock,*work,1.0,*currentBlock);
          delete work;
        }
        else
        {
          rightBlock = new SimpleMatrix(nslawSize2, sizeDS);
          UR2->getLeftBlockForDS(*itDS, rightBlock);
          // Warning: we use getLeft for Right block because right = transpose(left) and because of size checking inside the getBlock function,
          // a getRight call will fail.
          rightBlock->trans();
          centralBlocks[*itDS]->PLUForwardBackwardInPlace(*rightBlock);
          //*currentBlock +=  *leftBlock ** work;
          prod(*leftBlock, *rightBlock, *currentBlock, false);
          delete rightBlock;
        }
      }
      else RuntimeException::selfThrow("LCP::computeBlock not yet implemented for relation of type " + relationType1);
      delete leftBlock;
    }
  }
}

void LCP::computeQ(double time)
{
  if (q->size() != sizeOutput)
    q->resize(sizeOutput);
  q->zero();

  // === Get index set from Simulation ===
  UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);
  // === Loop through "active" Unitary Relations (ie present in indexSets[level]) ===

  unsigned int pos = 0;
  UnitaryRelationsIterator itCurrent, itLinked;
  string simulationType = simulation->getType();
  for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // *itCurrent is a UnitaryRelation*.

    // Compute q, this depends on the type of non smooth problem, on the relation type and on the non smooth law
    pos = M->getPositionOfBlock(*itCurrent);

#ifndef WithSmartPtr
    (*itCurrent)->computeEquivalentY(time, levelMin, simulationType, q, pos);
#else
    (*itCurrent)->computeEquivalentY(time, levelMin, simulationType, q.get(), pos);
#endif
  }
}

void LCP::preCompute(double time)
{
  // This function is used to prepare data for the LinearComplementarity_Problem
  // - computation of M and q
  // - set sizeOutput
  // - check dim. for z,w

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  Topology * topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  if (!topology->isTimeInvariant())
  {
    // Computes new blocks if required
    updateBlocks();

    // Updates matrix M
    UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);
    M->fill(indexSet, blocks);
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

  // Computes q of LCP
  computeQ(time);

}

int LCP::compute(double time)
{
  // --- Prepare data for LCP computing ---
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
    // The LCP in Numerics format
    LinearComplementarity_Problem numerics_problem;
    numerics_problem.M = M->getNumericsMatrix();
    numerics_problem.q = q->getArray();
    numerics_problem.size = sizeOutput;
    int nbSolvers = 1;
    // Call LCP Driver
    info = lcp_driver(&numerics_problem, z->getArray() , w->getArray() , solver->getNumericsSolverOptionsPtr(), nbSolvers, numerics_options);

    // --- Recovering of the desired variables from LCP output ---
    postCompute();

  }

  return info;
}

void LCP::postCompute()
{
  // This function is used to set y/lambda values using output from lcp_driver (w,z).
  // Only UnitaryRelations (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);

  // y and lambda vectors
  SiconosVector *lambda, *y;

  // === Loop through "active" Unitary Relations (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;

  for (UnitaryRelationsIterator itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // size of the block that corresponds to the current UnitaryRelation
    nsLawSize = (*itCurrent)->getNonSmoothLawSize();
    // Get the relative position of UR-block in the vector w or z
    pos = M->getPositionOfBlock(*itCurrent);

    // Get Y and Lambda for the current Unitary Relation
    y = (*itCurrent)->getYPtr(levelMin);
    lambda = (*itCurrent)->getLambdaPtr(levelMin);
    // Copy w/z values, starting from index pos into y/lambda.

#ifndef WithSmartPtr
    setBlock(w, y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
    setBlock(z, lambda, lambda->size(), pos, 0);
#else
    setBlock(w.get(), y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
    setBlock(z.get(), lambda, lambda->size(), pos, 0);
#endif

  }
}

void LCP::display() const
{
  cout << "======= LCP of size " << sizeOutput << " with: " << endl;
  cout << "M  ";
  if (M) M->display();
  else cout << "-> NULL" << endl;
  cout << endl << " q : " ;
  if (q) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;
}

void LCP::saveNSProblemToXML()
{
  OneStepNSProblem::saveNSProblemToXML();
  //   if(onestepnspbxml != NULL)
  //     {
  // //       (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
  //       (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  //     }
  //   else RuntimeException::selfThrow("LCP::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  RuntimeException::selfThrow("LCP::saveNSProblemToXML - Not yet implemented.");
}

LCP* LCP::convert(OneStepNSProblem* osnsp)
{
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}



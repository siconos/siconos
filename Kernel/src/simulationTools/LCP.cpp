/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
#include "LCPXML.h"
#include "Topology.h"
#include "UnitaryRelation.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "Relation.h"
#include "DynamicalSystem.h"
#include "TimeDiscretisation.h"

using namespace std;

// xml constructor
LCP::LCP(OneStepNSProblemXML* onestepnspbxml, Simulation* newSimu):
  OneStepNSProblem("LCP", onestepnspbxml, newSimu), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false) ,
  MSparseBlock(NULL), prepareM(NULL)
{
  LCPXML * xmllcp = (static_cast<LCPXML*>(onestepnspbxml));

  // If both q and M are given in xml file, check if sizes are consistent
  if (xmllcp->hasQ() && xmllcp->hasM() && ((xmllcp->getM()).size(0) != (xmllcp->getQ()).size()))
    RuntimeException::selfThrow("LCP: xml constructor, inconsistent sizes between given q and M");

  // The dim of the problem is given by M if given in the xml input file.
  if (xmllcp->hasM())
  {
    sizeOutput = (xmllcp->getM()).size(0);
    M = new SimpleMatrix(xmllcp->getM());
    isMAllocatedIn = true;
  }

  if (xmllcp->hasQ())
  {
    // get sizeOutput if necessary
    if (M == NULL)
      sizeOutput = (xmllcp->getQ()).size();

    q = new SimpleVector(xmllcp->getQ());
    isQAllocatedIn = true;
  }

  // Note that M and q are allocated only if given in the xml file.

}

// Constructor from a set of data (appart from Simulation and id, all arguments are optional)
LCP::LCP(Simulation* newSimu, const string& newId, const string& newSolver, unsigned int MaxIter,
         double Tolerance, unsigned int Verbose,  const string&  NormType, double  SearchDirection, double Rho):
  OneStepNSProblem("LCP", newSimu, newId), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false), MSparseBlock(NULL), prepareM(NULL)
{
  // set solver:
  solver = new Solver(nspbType, newSolver, MaxIter, Tolerance, Verbose, NormType, SearchDirection, Rho);
  isSolverAllocatedIn = true;
}

// Constructor from a set of data
LCP::LCP(Solver* newSolver, Simulation* newSimu, const string& newId):
  OneStepNSProblem("LCP", newSimu, newId, newSolver), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false), MSparseBlock(NULL), prepareM(NULL)
{}

// destructor
LCP::~LCP()
{
  if (isWAllocatedIn) delete w;
  w = NULL;
  if (isZAllocatedIn) delete z;
  z = NULL;
  if (isMAllocatedIn)
  {
    if (solver->useBlocks())
    {
      delete MSparseBlock;
    }
    else delete M;
  }
  M = NULL;
  MSparseBlock = NULL;
  if (isQAllocatedIn) delete q;
  q = NULL;
  prepareM = NULL;
}

// Setters

void LCP::setW(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("LCP: setW, inconsistent size between given w size and problem size. You should set sizeOutput before");

  if (w == NULL)
  {
    w = new SimpleVector(sizeOutput);
    isWAllocatedIn = true;
  }
  else if (w->size() != sizeOutput)
    RuntimeException::selfThrow("LCP: setW, w size differs from sizeOutput");

  *w = newValue;
}

void LCP::setWPtr(SiconosVector* newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("LCP: setWPtr, inconsistent size between given w size and problem size. You should set sizeOutput before");

  if (isWAllocatedIn) delete w;
  w = newPtr;
  isWAllocatedIn = false;
}


void LCP::setZ(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("LCP: setZ, inconsistent size between given z size and problem size. You should set sizeOutput before");

  if (z == NULL)
  {
    z = new SimpleVector(sizeOutput);
    isZAllocatedIn = true;
  }

  *z = newValue;
}

void LCP::setZPtr(SiconosVector* newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("LCP: setZPtr, inconsistent size between given z size and problem size. You should set sizeOutput before");

  if (isZAllocatedIn) delete z;
  z = newPtr;
  isZAllocatedIn = false;
}

void LCP::setM(const SiconosMatrix& newValue)
{
  if (solver->useBlocks())
    RuntimeException::selfThrow("LCP: setM, impossible since M is sparse by blocks.");
  if (sizeOutput != newValue.size(0) || sizeOutput != newValue.size(1))
    RuntimeException::selfThrow("LCP: setM, inconsistent size between given M size and problem size. You should set sizeOutput before");

  if (M == NULL)
  {
    M = new SimpleMatrix(sizeOutput, sizeOutput);
    isMAllocatedIn = true;
  }

  *M = newValue;
}

void LCP::setMPtr(SiconosMatrix* newPtr)
{
  if (solver->useBlocks())
    RuntimeException::selfThrow("LCP: setMPtr, impossible since M is sparse by blocks.");
  if (sizeOutput != newPtr->size(0) || sizeOutput != newPtr->size(1))
    RuntimeException::selfThrow("LCP: setMPtr, inconsistent size between given M size and problem size. You should set sizeOutput before");

  if (isMAllocatedIn) delete M;
  M = newPtr;
  isMAllocatedIn = false;
}


void LCP::setQ(const SiconosVector& newValue)
{
  if (sizeOutput != newValue.size())
    RuntimeException::selfThrow("LCP: setQ, inconsistent size between given q size and problem size. You should set sizeOutput before");

  if (q == NULL)
  {
    q = new SimpleVector(sizeOutput);
    isQAllocatedIn = true;
  }

  *q = newValue;
}

void LCP::setQPtr(SiconosVector* newPtr)
{
  if (sizeOutput != newPtr->size())
    RuntimeException::selfThrow("LCP: setQPtr, inconsistent size between given q size and problem size. You should set sizeOutput before");

  if (isQAllocatedIn) delete q;
  q = newPtr;
  isQAllocatedIn = false;
}

void LCP::setSolverPtr(Solver * newSolv)
{
  OneStepNSProblem::setSolverPtr(newSolv);
  initSolver();
}

void LCP::initSolver()
{
  if (solver->useBlocks())
    prepareM = &LCP::collectMBlocks;
  else
    prepareM = &LCP::assembleM;
}

void LCP::initialize()
{
  // This function performs all steps that are not time dependant

  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize();

  // Memory allocation for w, M, z and q.
  // If one of them has already been allocated, nothing is done.
  // We suppose that user choose a correct size.
  if (w == NULL)
  {
    w = new SimpleVector(maxSize);
    isWAllocatedIn = true;
  }
  if (z == NULL)
  {
    z = new SimpleVector(maxSize);
    isZAllocatedIn = true;
  }
  if (M == NULL && !solver->useBlocks())
  {
    M = new SimpleMatrix(maxSize, maxSize);
    isMAllocatedIn = true;
  }

  if (MSparseBlock == NULL && solver->useBlocks())
  {
    int nc = simulation->getIndexSetPtr(0)->size();
    MSparseBlock = new SparseBlockMatrix(nc, nc);
    isMAllocatedIn = true;
  }

  if (q == NULL)
  {
    q = new SimpleVector(maxSize);
    isQAllocatedIn = true;
  }
  //   w->zero();
  //   z->zero();

  // if all relative degrees are equal to 0 or 1

  initSolver();

  // get topology
  Topology * topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  if (topology->isTimeInvariant() &&   !OSNSInteractions->isEmpty())
  {
    // computeSizeOutput and updateBlocks were already done in OneStepNSProblem::initialize
    if (sizeOutput != 0)
      (*this.*prepareM)();
  }

  // If topology is time-dependant, all the above commands are called during preCompute function.
}

void LCP::preCompute(double time)
{
  // compute M and q operators for LCP problem

  // get topology
  Topology * topology = simulation->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  if (!topology->isTimeInvariant())
    computeSizeOutput();

  if (sizeOutput != 0)
  {
    if (!topology->isTimeInvariant())
    {
      updateBlocks();
      // Assemble or collect M
      (*this.*prepareM)();

      // check z and w sizes and reset if necessary
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
    computeQ(time);
  }
}

void LCP::computeBlock(UnitaryRelation* UR1, UnitaryRelation* UR2)
{

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

void LCP::collectMBlocks()
{
  // === Description ===
  // This routine is used to fill in the SparseBlockStructuredMatrix Mspbl.
  //
  // For each Unitary Relation which is active (ie which is in indexSet[1] of the Simulation), named itRow (it corresponds to a row of block in M),
  // and for each Unitary Relation connected to itRow and active, named itCol, we get the block[itRow][itCol] and add the embedded double** into Mspbl.
  //
  if (!solver->useBlocks())
    RuntimeException::selfThrow("LCP::collectMBlocks failed. The linked solver does not fit to sparse block storage. Maybe you change your solver and forget to call initSolver.");

  // Get index set 1 from Simulation
  UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);

  MSparseBlock->fill(indexSet, blocks);
  MSparseBlock->convert();
}

void LCP::assembleM() //
{
  // === Description ===
  // For each Unitary Relation which is active (ie which is in indexSet[1] of the Simulation), named itRow (it corresponds to a row of block in M),
  // and for each Unitary Relation connected to itRow and active, named itCol, we get the block[itRow][itCol] and copy it at the right place in M matrix of the LCP
  //
  // Note that if a Unitary Relation is connected to another one (ie if they have common DS), the corresponding block must be in blocks map. This is the role of
  // updateBlocks() to ensure this.
  // But the presence of a block in the map does not imply that its corresponding UR are active, since this block may have been computed at a previous time step.
  // That is why we need to check if itCol is in indexSet1.
  //
  // See updateBlocks function for more details.

  // Get index set 1 from Simulation
  UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);
  UnitaryRelationsIterator itRow;
  unsigned int pos = 0, col = 0; // index position used for block copy into M, see below.
  UnitaryMatrixColumnIterator itCol;

  if (solver->useBlocks())
    RuntimeException::selfThrow("LCP::assembleM failed. The linked solver does not fit to non sparse block storage. Maybe you change your solver and forget to call initSolver.");

  // full matrix M assembling
  // === Memory allocation, if required ===
  // Mem. is allocate only if M==NULL or if its size has changed.

  if (M->size(0) != sizeOutput || M->size(1) != sizeOutput)
    M->resize(sizeOutput, sizeOutput);
  M->zero();
  // === Loop through "active" Unitary Relations (ie present in indexSets[level]) ===

  for (itRow = indexSet->begin(); itRow != indexSet->end(); ++itRow)
  {
    for (itCol = blocks[*itRow].begin(); itCol != blocks[*itRow].end(); ++itCol)
    {
      // *itRow and itCol are UnitaryRelation*.

      // Check that itCol is in Index set 1
      if ((indexSet->find((*itCol).first)) != indexSet->end())
      {
        pos = blocksPositions[*itRow];
        col = blocksPositions[(*itCol).first];
        // copy the block into Mlcp - pos/col: position in M (row and column) of first element of the copied block
        static_cast<SimpleMatrix*>(M)->setBlock(pos, col, *(blocks[*itRow][(*itCol).first])); // \todo avoid copy
      }
    }
  }
  // === end of UnitaryRelations loop ===
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
    pos = blocksPositions[*itCurrent];
    (*itCurrent)->computeEquivalentY(time, levelMin, simulationType, q, pos);
  }
}

int LCP::compute(double time)
{
  // --- Prepare data for LCP computing ---
  preCompute(time);

  int info = 0;
  // --- Call Numerics solver ---
  if (sizeOutput != 0)
  {
    int Nlcp = (int)sizeOutput;
    int iter, titer;
    double err;
    clock_t startLCPsolve;

    method solvingMethod = *(solver->getSolvingMethodPtr());

    startLCPsolve = clock();
    if (solver->useBlocks()) // Use solver block
    {
      info = lcp_solver_block(MSparseBlock->getNumericsMatSparse() , q->getArray() , &solvingMethod , z->getArray() , w->getArray() , &iter , &titer , &err);
    }

    else // Use classical solver
      info = lcp_solver(M->getArray(), q->getArray(), &Nlcp, &solvingMethod, z->getArray(), w->getArray());

    CPUtime += (clock() - startLCPsolve);
    nbIter++;
    // Remark : the output result, info, from solver is treated when this function is call from Simulation

    // --- Recover the desired variables from LCP output ---
    postCompute();
  }
  return info;
}

void LCP::postCompute()
{
  // === Get index set from Topology ===
  UnitaryRelationsSet * indexSet = simulation->getIndexSetPtr(levelMin);

  // y and lambda vectors
  SiconosVector *lambda, *y;

  // === Loop through "active" Unitary Relations (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;
  UnitaryRelationsIterator itCurrent, itLinked;

  for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // *itCurrent is a UnitaryRelation*.
    // size if a block that corresponds to the current UnitaryRelation
    nsLawSize = (*itCurrent)->getNonSmoothLawSize();
    // Get the relative position of UR-block in the vector w or z
    pos = blocksPositions[*itCurrent];

    // Get Y and Lambda for the current Unitary Relation
    y = (*itCurrent)->getYPtr(levelMin);
    lambda = (*itCurrent)->getLambdaPtr(levelMin);
    // Copy w/z values, starting from index pos into y/lambda.
    setBlock(w, y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
    setBlock(z, lambda, lambda->size(), pos, 0);
  }
}

void LCP::display() const
{
  cout << "======= LCP of size " << sizeOutput << " with: " << endl;
  cout << "M  ";
  if (solver->useBlocks())
  {
    cout << "a Sparse Block matrix." << endl;
    if (MSparseBlock != NULL)
      MSparseBlock->display();
    else cout << "M is NULL" << endl;
  }
  else
  {
    if (M != NULL) M->display();
    else cout << "-> NULL" << endl;
  }
  cout << endl << " q : " ;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;
}

void LCP::saveNSProblemToXML()
{
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
    (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("LCP::saveNSProblemToXML - OneStepNSProblemXML object not exists");
}

void LCP::saveMToXML()
{
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("LCP::saveMToXML - OneStepNSProblemXML object not exists");
}

void LCP::saveQToXML()
{
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("LCP::saveQToXML - OneStepNSProblemXML object not exists");
}

LCP* LCP::convert(OneStepNSProblem* osnsp)
{
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}



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
#include "MixedComplementarityConditionNSL.h"
#include "Simulation.h"

using namespace std;
using namespace RELATION;

// xml constructor
MLCP::MLCP(SP::OneStepNSProblemXML onestepnspbxml):
  LinearOSNS(onestepnspbxml, "MLCP") {}

// Constructor from a set of data
MLCP::MLCP(SP::NonSmoothSolver newSolver, const string& newId):
  LinearOSNS("MLCP", newSolver, newId) {}

void MLCP::updateM()
{
  // Get index set from Simulation
  SP::UnitaryRelationsGraph indexSet = simulation->getIndexSetPtr(levelMin);

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
  mlcp_driver_reset(&numerics_problem, (solver->getNumericsSolverOptionsPtr()).get());
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
      (unitaryBlocks[UR1][UR2]).reset(new SimpleMatrix(nslawSize1, nslawSize2));

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
    //      printf("h : %f \n",h);

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

      }
      else RuntimeException::selfThrow("MLCP::computeBlock not yet implemented for relation of type " + relationType1);

    }

  }
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
    info = mlcp_driver(&numerics_problem, _z->getArray(), _w->getArray(), (solver->getNumericsSolverOptionsPtr()).get(), &*numerics_options);

    // --- Recovering of the desired variables from MLCP output ---
    postCompute();

  }

  return info;
}

void MLCP::display() const
{
  cout << "======= MLCP of size " << sizeOutput << " with: " << endl;
  cout << "======= m " << m << " n " << n << endl;
  LinearOSNS::display();
}

MLCP* MLCP::convert(OneStepNSProblem* osnsp)
{
  MLCP* lcp = dynamic_cast<MLCP*>(osnsp);
  return lcp;
}



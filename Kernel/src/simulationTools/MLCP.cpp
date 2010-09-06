/* Siconos-Kernel, Copyright INRIA 2005-2010.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */
#include "MLCP.hpp"
#include "MixedComplementarityConditionNSL.hpp"
#include "EqualityConditionNSL.hpp"
#include "Simulation.hpp"

using namespace std;
using namespace RELATION;

// xml constructor
MLCP::MLCP(SP::OneStepNSProblemXML onestepnspbxml):
  LinearOSNS(onestepnspbxml, "MLCP") {}

// Constructor from a set of data
MLCP::MLCP(const int newNumericsSolverId):
  LinearOSNS(newNumericsSolverId, "MLCP", "unamed")
{
  _numerics_solver_options->solverId = newNumericsSolverId;
  mixedLinearComplementarity_setDefaultSolverOptions(NULL, &*_numerics_solver_options);

  _numerics_problem.blocksLine = (int*)malloc(MLCP_NB_BLOCKS * sizeof(int));
  _numerics_problem.blocksIsComp = (int*)malloc(MLCP_NB_BLOCKS * sizeof(int));
  _numerics_problem.blocksLine[0] = 0;
  _curBlock = 0;

}

void MLCP::updateM()
{
  // Get index set from Simulation
  SP::UnitaryRelationsGraph indexSet = simulation()->indexSet(levelMin());

  if (!_M)
  {
    // Creates and fills M using UR of indexSet
    _M.reset(new OSNSMatrix(indexSet, _unitaryBlocks, _MStorageType));
    _numerics_problem.M = &*_M->getNumericsMatrix();
    _numerics_problem.A = 0;
    _numerics_problem.B = 0;
    _numerics_problem.C = 0;
    _numerics_problem.D = 0;
    _numerics_problem.a = 0;
    _numerics_problem.b = 0;
    _numerics_problem.problemType = 0;
    _numerics_problem.n = _n;
    _numerics_problem.m = _m;
  }
  else
  {
    _M->setStorageType(_MStorageType);
    _M->fill(indexSet, _unitaryBlocks);

  }
  _sizeOutput = _M->size();
}

void  MLCP::reset()
{
  if (_numerics_problem.blocksLine)
    free(_numerics_problem.blocksLine);
  _numerics_problem.blocksLine = 0;
  if (_numerics_problem.blocksIsComp)
    free(_numerics_problem.blocksIsComp);
  _numerics_problem.blocksIsComp = 0;
  mlcp_driver_reset(&_numerics_problem, &*_numerics_solver_options);
}


void MLCP::computeOptions(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{
  DynamicalSystemsSet commonDS;
  intersection(*UR1->dynamicalSystems(), *UR2->dynamicalSystems(), commonDS);
  assert(!commonDS.isEmpty()) ;

  // Get dimension of the NonSmoothLaw (ie dim of the unitaryBlock)
  unsigned int nslawSize1 = UR1->getNonSmoothLawSize();
  unsigned int nslawSize2 = UR2->getNonSmoothLawSize();

  unsigned int equalitySize1 =  0;
  unsigned int equalitySize2 =  0;
  if (Type::value(*(UR1->interaction()->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize1 =  MixedComplementarityConditionNSL::convert(UR1->interaction()->nonSmoothLaw())->getEqualitySize();
  else if (Type::value(*(UR1->interaction()->nonSmoothLaw()))
           == Type::EqualityConditionNSL)
    equalitySize1 = nslawSize1;

  if (Type::value(*(UR2->interaction()->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize2 = MixedComplementarityConditionNSL::
                    convert(UR2->interaction()->nonSmoothLaw())->getEqualitySize();
  else if (Type::value(*(UR2->interaction()->nonSmoothLaw()))
           == Type::EqualityConditionNSL)
    equalitySize2 = nslawSize2;


  if (UR1 == UR2)
  {
    //UR1->getExtraUnitaryBlock(currentUnitaryBlock);
    _m += nslawSize1 - equalitySize1;
    _n += equalitySize1;
    if (_curBlock > MLCP_NB_BLOCKS - 2)
      printf("MLCP.cpp : number of block to small, memory crach below!!!\n");
    /*add an equality block.*/
    if (equalitySize1 > 0)
    {
      _numerics_problem.blocksLine[_curBlock + 1] = _numerics_problem.blocksLine[_curBlock] + equalitySize1;
      _numerics_problem.blocksIsComp[_curBlock] = 0;
      _curBlock++;
    }
    /*add a complementarity block.*/
    if (nslawSize1 - equalitySize1 > 0)
    {
      _numerics_problem.blocksLine[_curBlock + 1] = _numerics_problem.blocksLine[_curBlock] + nslawSize1 - equalitySize1;
      _numerics_problem.blocksIsComp[_curBlock] = 1;
      _curBlock++;
    }
  }
}

void MLCP::computeUnitaryBlock(SP::UnitaryRelation UR1, SP::UnitaryRelation UR2)
{
  computeOptions(UR1, UR2);
  LinearOSNS::computeUnitaryBlock(UR1, UR2);
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
  _numerics_problem.n = _n;
  _numerics_problem.m = _m;
  /*If user has not allocted the working memory, do it. */
  int allocated = mlcp_alloc_working_memory(&_numerics_problem, &*_numerics_solver_options);
  int info = 0;
  // --- Call Numerics driver ---
  // Inputs:
  // - the problem (M,q ...)
  // - the unknowns (z,w)
  // - the options for the solver (name, max iteration number ...)
  // - the global options for Numerics (verbose mode ...)

  if (_sizeOutput != 0)
  {
    _numerics_problem.q = _q->getArray();
    int nbSolvers = 1;
    // Call MLCP Driver
    //printf("MLCP display");
    //printf("n %d m %d",n,m);
    //displayNM(_numerics_problem.M);
    //      exit(1);
    //mlcpDefaultSolver *pSolver = new mlcpDefaultSolver(m,n);
    //displayMLCP(&_numerics_problem);
    try
    {
      //  display();
      info = mlcp_driver(&_numerics_problem, _z->getArray(), _w->getArray(),
                         &*_numerics_solver_options, &*_numerics_options);
    }
    catch (...)
    {
      cout << "exception catched" << endl;
      info = 1;
    }

    // --- Recovering of the desired variables from MLCP output ---
    postCompute();

  }
  if (allocated)
    mlcp_free_working_memory(&_numerics_problem, &*_numerics_solver_options);
  return info;
}

void MLCP::display() const
{
  cout << "======= MLCP of size " << _sizeOutput << " with: " << endl;
  cout << "======= m " << _m << " _n " << _n << endl;
  LinearOSNS::display();
}

MLCP* MLCP::convert(OneStepNSProblem* osnsp)
{
  MLCP* mlcp = dynamic_cast<MLCP*>(osnsp);
  return mlcp;
}

void MLCP::initialize(SP::Simulation sim)
{
  // General initialize for LinearOSNS
  LinearOSNS::initialize(sim);



}
void  MLCP::updateUnitaryBlocks()
{
  _curBlock = 0;
  _m = 0;
  _n = 0;
  LinearOSNS::updateUnitaryBlocks();
}
void  MLCP::computeAllUnitaryBlocks()
{
  _curBlock = 0;
  _m = 0;
  _n = 0;
  LinearOSNS::computeAllUnitaryBlocks();
}

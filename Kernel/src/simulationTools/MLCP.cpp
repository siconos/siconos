/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
//#define MLCP_DEBUG
// xml constructor
MLCP::MLCP(SP::OneStepNSProblemXML onestepnspbxml):
  LinearOSNS(onestepnspbxml, "MLCP")
{
  _n = 0;
  _m = 0;
}

// Constructor from a set of data
MLCP::MLCP(const int newNumericsSolverId):
  LinearOSNS(newNumericsSolverId, "MLCP", "unnamed")
{
  _numerics_solver_options->solverId = newNumericsSolverId;
  mixedLinearComplementarity_setDefaultSolverOptions(NULL, &*_numerics_solver_options);
  _n = 0;
  _m = 0;
  _numerics_problem.blocksLine = (int*)malloc(MLCP_NB_BLOCKS * sizeof(int));
  _numerics_problem.blocksIsComp = (int*)malloc(MLCP_NB_BLOCKS * sizeof(int));
  _numerics_problem.blocksLine[0] = 0;
  _curBlock = 0;

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


void MLCP::computeOptions(SP::Interaction inter1, SP::Interaction inter2)
{
  // Get dimension of the NonSmoothLaw (ie dim of the interactionBlock)
  unsigned int nslawSize1 = inter1->getNonSmoothLawSize();
  unsigned int nslawSize2 = inter2->getNonSmoothLawSize();

  unsigned int equalitySize1 =  0;
  unsigned int equalitySize2 =  0;
  if (Type::value(*(inter1->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize1 =  MixedComplementarityConditionNSL::convert(inter1->nonSmoothLaw())->getEqualitySize();
  else if (Type::value(*(inter1->nonSmoothLaw()))
           == Type::EqualityConditionNSL)
    equalitySize1 = nslawSize1;

  if (Type::value(*(inter2->nonSmoothLaw()))
      == Type::MixedComplementarityConditionNSL)
    equalitySize2 = MixedComplementarityConditionNSL::
                    convert(inter2->nonSmoothLaw())->getEqualitySize();
  else if (Type::value(*(inter2->nonSmoothLaw()))
           == Type::EqualityConditionNSL)
    equalitySize2 = nslawSize2;


  if (inter1 == inter2)
  {
    //inter1->getExtraInteractionBlock(currentInteractionBlock);
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

void MLCP::computeInteractionBlock(const InteractionsGraph::EDescriptor& ed)
{

  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());
  SP::Interaction inter1 = indexSet->bundle(indexSet->source(ed));
  SP::Interaction inter2 = indexSet->bundle(indexSet->target(ed));

  assert(inter1 != inter2);
  bool isLinear = simulation()->model()->nonSmoothDynamicalSystem()->isLinear();

  if (!_hasBeenUpdated || !isLinear)
    LinearOSNS::computeInteractionBlock(ed);
}

void MLCP::computeDiagonalInteractionBlock(const InteractionsGraph::VDescriptor& vd)
{

  SP::InteractionsGraph indexSet = simulation()->indexSet(levelMin());
  SP::DynamicalSystem DS1 = indexSet->properties(vd).source;
  SP::DynamicalSystem DS2 = indexSet->properties(vd).target;
  SP::Interaction inter = indexSet->bundle(vd);

  // commonDS here...
  if (!_hasBeenUpdated)
    computeOptions(inter, inter);
  LinearOSNS::computeDiagonalInteractionBlock(vd);
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

void MLCP::preCompute(double time)
{
  LinearOSNS::preCompute(time);
  _numerics_problem.n = _n;
  _numerics_problem.m = _m;
}

int MLCP::compute(double time)
{
  // --- Prepare data for MLCP computing ---
  preCompute(time);
  //  _hasBeenUpdated=true;
#ifdef MLCP_DEBUG
  printf("MLCP::compute m n :%d,%d\n", _n, _m);
#endif
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

    // Call MLCP Driver
    //printf("MLCP display");
    //printf("n %d m %d",n,m);
    //displayNM(_numerics_problem.M);
    //      exit(1);
    //mlcpDefaultSolver *pSolver = new mlcpDefaultSolver(m,n);
    try
    {
      info = mlcp_driver(&_numerics_problem, _z->getArray(), _w->getArray(),
                         &*_numerics_solver_options, &*_numerics_options);
#ifdef MLCP_DEBUG
      display();
#endif
    }
    catch (...)
    {
      cout << "exception catched" << endl;
      info = 1;
    }

    // --- Recovering of the desired variables from MLCP output ---
    if (!info)
      postCompute();
    else
      printf("MLCP solver failed\n");

  }
  else
  {
#ifdef MLCP_DEBUG
    printf("MLCP::compute : sizeoutput is null\n");
#endif
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

  _numerics_problem.M = &*_M->getNumericsMatrix();
  _numerics_problem.A = 0;
  _numerics_problem.B = 0;
  _numerics_problem.C = 0;
  _numerics_problem.D = 0;
  _numerics_problem.a = 0;
  _numerics_problem.b = 0;
  _numerics_problem.problemType = 0;
}
void  MLCP::updateInteractionBlocks()
{
  if (!_hasBeenUpdated)
  {
    _curBlock = 0;
    _m = 0;
    _n = 0;
  }
  LinearOSNS::updateInteractionBlocks();
}
void  MLCP::computeAllInteractionBlocks()
{
  assert(0);
  _curBlock = 0;
  _m = 0;
  _n = 0;
  LinearOSNS::computeAllInteractionBlocks();
}

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
#include "MLCP2.hpp"
#include "Moreau2.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepIntegrator.hpp"
#include "FirstOrderLinearR.hpp"
#include "Topology.hpp"

using namespace std;

// Constructor from a set of data
MLCP2::MLCP2(const string& newNumericsSolverName, const string& newId):
  MLCP(newNumericsSolverName, newId)
{
  mFirstCall = true;
  m = 0;
  n = 0;

}

void MLCP2::initialize(SP::Simulation simulation)
{
  SP::DynamicalSystemsSet  DSSet = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  if (DSSet->begin() == DSSet->end())
    printf("DSSet is empty\n");
  else
    printf("DSSet is not empty\n");
  updateDSBlocks();
  updateDSInteractionBlocks();
  DSSet = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  if (DSSet->begin() == DSSet->end())
    printf("DSSet is empty\n");
  else
    printf("DSSet is not empty\n");
  updateInteractionDSBlocks();
  DSSet = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  if (DSSet->begin() == DSSet->end())
    printf("DSSet is empty\n");
  else
    printf("DSSet is not empty\n");
  MLCP::initialize(simulation);
  _w->resize(n + m);
  _z->resize(n + m);
}
void MLCP2::updateM()
{
  SP::InteractionsGraph IG = simulation->indexSet(levelMin);
  SP::DynamicalSystemsSet DSSet = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  if (!M)
  {
    // Creates and fills M using Interactionof indexSet
    M.reset(new OSNSMatrix(IG, DSSet, interactionBlocks,  DSBlocks, DSInteractionBlocks,  interactionDSBlocks, MStorageType));

    numerics_problem.M = M->getNumericsMatrix().get();
    numerics_problem.A = 0;
    numerics_problem.B = 0;
    numerics_problem.C = 0;
    numerics_problem.D = 0;
    numerics_problem.a = 0;
    numerics_problem.b = 0;
    numerics_problem.problemType = 0;
    numerics_problem.n = n; //size of the dynamical systems
    numerics_problem.m = m; //size of the NSlaw
  }
  else
  {
    M->setStorageType(MStorageType);
    M->fill(IG, interactionBlocks);

  }
  sizeOutput = M->size();
}

// void MLCP2::updateInteractionBlocks()
// {
//   SP::InteractionsSet indexSet;
//   bool isTimeInvariant;
//   InteractionsIterator itinter1, itinter2;
//   DSIterator itDS;

//   // Get index set from Simulation

//   indexSet = simulation->indexSet(levelMin);
//   isTimeInvariant = simulation->model()->nonSmoothDynamicalSystem()->topology()->isTimeInvariant();

//   if (!isTimeInvariant || mFirstCall){
//     for(itinter1 = indexSet->begin(); itinter1!=indexSet->end();++itinter1){
//       computeblock(&(*itinter1));
//     }
//     DynamicalSystemsSet * allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
//     for(itDS=allDS->begin(); itDS!= allDS->end(); ++itDS){
//       computeblock(&(*itDS));
//       for(itinter1 = indexSet->begin(); itinter1!=indexSet->end();++itinter1){
//  computeblock(&(*itDS),&(*itinter1));
//  computeblock(&(*itinter1),&(*itDS));
//       }
//     }
//   }
//   mFirstCall = false;

// }


// // fill BlocksDS : W from the DS
void MLCP2::computeDSBlock(SP::DynamicalSystem DS)
{
  if (!DSBlocks[DS])
  {
    DSBlocks[DS].reset(new SimpleMatrix(DS->getDim(), DS->getDim()));
    n += DS->getDim();
  }
  SP::OneStepIntegrator Osi = simulation->integratorOfDS(DS); // get OneStepIntegrator of current dynamical system
  const OSI::TYPES  osiType = Osi->getType();
  if (osiType == OSI::MOREAU2)
  {
    DSBlocks[DS] = (boost::static_pointer_cast<Moreau2> (Osi))->W(DS); // get its W matrix ( pointer link!)
    (*DSBlocks[DS]) *= -1.0;
    DSBlocks[DS]->display();
  }
}

// // fill interactionBlocks : D
void MLCP2::computeInteractionBlock(SP::Interaction inter1, SP::Interaction inter2)
{
  if (!interactionBlocks[inter1][inter1])
  {
    interactionBlocks[inter1][inter1].reset(new SimpleMatrix(inter1->getNonSmoothLawSize(), inter1->getNonSmoothLawSize()));
    m += inter1->getNonSmoothLawSize();
  }
  inter1->getExtraInteractionBlock(interactionBlocks[inter1][inter1]);
}


// fill interactionDSBlocks : C
void MLCP2::computeInteractionDSBlock(SP::Interaction inter, SP::DynamicalSystem DS)
{
  if (interactionDSBlocks[inter][DS] == 0)
  {
    interactionDSBlocks[inter][DS].reset(new SimpleMatrix(inter->getNonSmoothLawSize(), DS->getDim()));
  }
  inter->getLeftInteractionBlockForDS(DS, interactionDSBlocks[inter][DS]);
}

//fill DSInteractionBlocks with B
void MLCP2::computeDSInteractionBlock(SP::DynamicalSystem DS, SP::Interaction inter)
{
  double h = simulation->timeStep();
  if (!DSInteractionBlocks[DS][inter])
  {
    DSInteractionBlocks[DS][inter].reset(new SimpleMatrix(DS->getDim(), inter->getNonSmoothLawSize()));
  }
  inter->getRightInteractionBlockForDS(DS, DSInteractionBlocks[DS][inter]);
  *(DSInteractionBlocks[DS][inter]) *= h;

}


void MLCP2::computeq(double time)
{
  if (q->size() != sizeOutput)
    q->resize(sizeOutput);
  q->zero();

  // === Get index set from Simulation ===
  SP::InteractionsGraph indexSet = simulation->indexSet(levelMin);
  // === Loop through "active" Interactions (ie present in indexSets[level]) ===

  unsigned int pos = 0;
  InteractionsGraph::VIterator ui, uiend;
  string simulationType = simulation->getType();
  for (itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // *itCurrent is a SP::Interaction.
    // Compute q, this depends on the type of non smooth problem, on the relation type and on the non smooth law
    pos = M->getPositionOfInteractionBlock(*itCurrent);
    //update e(ti+1)
    SP::SiconosVector  e = boost::static_pointer_cast<FirstOrderLinearR>((*itCurrent)->relation())->e();
    boost::static_pointer_cast<SiconosVector>(q)->addBlock(pos, *e);
  }
  SP::DynamicalSystemsSet  allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (DSIterator itDS = allDS->begin(); itDS != allDS->end(); ++itDS)
  {
    pos = M->getPositionOfDSBlock(*itDS);


    SP::OneStepIntegrator Osi = simulation->integratorOfDS(*itDS); // get OneStepIntegrator of current dynamical system
    const OSI::TYPES osiType = Osi->getType();
    if (osiType == OSI::MOREAU2)
    {
      //update with ffree
      SP::DynamicalSystem  DSaux = *itDS;
      SP::SiconosVector  Vaux = Osi->getWorkX(DSaux);
      boost::static_pointer_cast<SiconosVector>(q)->addBlock(pos, *Vaux);
    }

  }
}

void MLCP2::preCompute(double time)
{
  // This function is used to prepare data for the MixedLinearComplementarityProblem
  // - computation of M and q
  // - set sizeOutput
  // - check dim. for z,w

  // If the topology is time-invariant, only q needs to be computed at each time step.
  // M, sizeOutput have been computed in initialize and are uptodate.

  // Get topology
  SP::Topology topology = simulation->model()->nonSmoothDynamicalSystem()->topology();

  if (!topology->isTimeInvariant() || !simulation->model()->nonSmoothDynamicalSystem()->isLinear())
  {
    // Computes new interactionBlocks if required
    updateInteractionBlocks();

    // Updates matrix M
    SP::InteractionsSet interSet = simulation->indexSet(levelMin);
    SP::DynamicalSystemsSet DSSet = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
    //fill M block
    M->fill(interSet, DSSet, interactionBlocks, DSBlocks, DSInteractionBlocks, interactionDSBlocks);
    sizeOutput = M->size();

    // Checks z and w sizes and reset if necessary
    if (_z->size() != sizeOutput)
    {
      _z->resize(sizeOutput, false);
      _z->zero();
    }

    if (_w->size() != sizeOutput)
    {
      _w->resize(sizeOutput);
      _w->zero();
    }
  }

  // Computes q of MLCP2
  computeq(time);

}
void displayNM_(const NumericsMatrix* const m)
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

int MLCP2::compute(double time)
{
  // --- Prepare data for MLCP2 computing ---
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
    // Call MLCP2 Driver
    //printf("MLCP2 display");
    //printf("n %d m %d",n,m);
    //displayNM_(numerics_problem.M);
    //      exit(1);
    //mlcpDefaultSolver *pSolver = new mlcpDefaultSolver(m,n);
    //      displayMLCP2(&numerics_problem);
    info = mlcp_driver(&numerics_problem, _z->getArray(), _w->getArray(), (solver->numericsSolverOptions()).get(), numerics_options.get());

    // --- Recovering of the desired variables from MLCP2 output ---
    postCompute();

  }

  return info;
}

void MLCP2::postCompute()
{
  // This function is used to set y/lambda values using output from lcp_driver (w,z).
  // Only Interactions (ie Interactions) of indexSet(leveMin) are concerned.

  // === Get index set from Topology ===
  SP::InteractionsSet indexSet = simulation->indexSet(levelMin);

  // y and lambda vectors
  SP::SiconosVector lambda, y, x;

  // === Loop through "active" Interactions (ie present in indexSets[1]) ===

  unsigned int pos = 0;
  unsigned int nsLawSize;

  for (InteractionsIterator itCurrent = indexSet->begin(); itCurrent !=  indexSet->end(); ++itCurrent)
  {
    // size of the block that corresponds to the current Interaction
    nsLawSize = (*itCurrent)->getNonSmoothLawSize();
    // Get the relative position of inter-block in the vector w or z
    pos = M->getPositionOfInteractionBlock(*itCurrent);

    // Get Y and Lambda for the current Interaction
    y = (*itCurrent)->y(levelMin);
    lambda = (*itCurrent)->lambda(levelMin);
    // Copy w/z values, starting from index pos into y/lambda.
    setBlock(*(_w.get()), y, y->size(), pos, 0);// Warning: yEquivalent is saved in y !!
    setBlock(*(_z.get()), lambda, lambda->size(), pos, 0);
  }
  SP::DynamicalSystemsSet allDS = simulation->model()->nonSmoothDynamicalSystem()->dynamicalSystems();
  for (DSIterator itDS = allDS->begin(); itDS != allDS->end(); ++itDS)
  {
    pos = M->getPositionOfDSBlock(*itDS);
    x = (*itDS)->x();
    setBlock(*(_z.get()), x, x->size(), pos, 0);
  }

}

MLCP2* MLCP2::convert(OneStepNSProblem* osnsp)
{
  MLCP2* lcp = dynamic_cast<MLCP2*>(osnsp);
  return lcp;
}


void MLCP2::initialize(SP::Simulation sim)
{
  // General initialize for LinearOSNS
  Linear::initialize(sim);


}

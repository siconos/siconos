/* Siconos-Kernel, Copyright INRIA 2005-2012.
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

#include "EventDriven.hpp"
#include "SimulationXML.hpp"
#include "OneStepNSProblemXML.hpp"
#include "SimulationXML.hpp"
#include "Lsodar.hpp"
#include "LCP.hpp"
#include "Model.hpp"
#include "Interaction.hpp"
#include "EventsManager.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Topology.hpp"
#include "DynamicalSystem.hpp"
#include "LagrangianDS.hpp"
#include "EventFactory.hpp"
#include "BlockMatrix.hpp"
#include "NewMarkAlphaOSI.hpp"

#define DEBUG_MESSAGES

#include <debug.h>

using namespace std;
using namespace RELATION;

// --- XML constructor ---
EventDriven::EventDriven(SP::SimulationXML strxml, double t0, double T,
                         SP::DynamicalSystemsSet dsList,
                         SP::InteractionsSet interList):
  Simulation(strxml, t0, T, dsList, interList), _istate(1)
{
  // === One Step NS Problem === We read data in the xml output
  // (mainly Interactions concerned and solver) and assign them to
  // both one step ns problems ("acceleration" and "impact").  At the
  // time, only LCP is available for Event Driven.

  if (_simulationxml->hasOneStepNSProblemXML()) // ie if OSNSList is not empty
  {
    SetOfOSNSPBXML OSNSList = _simulationxml->getOneStepNSProblemsXML();
    SP::OneStepNSProblemXML osnsXML;
    string type;
    // For EventDriven, two OSNSPb are required, "acceleration" and
    // "impact"
    _numberOfOneStepNSproblems = 2;
    (*_allNSProblems).resize(_numberOfOneStepNSproblems);
    if (OSNSList.size() != 2)
      RuntimeException::selfThrow("EventDriven::xml constructor - Wrong number of OSNS problems: 2 are required.");
    int id = SICONOS_OSNSP_ED_IMPACT;
    for (SetOfOSNSPBXMLIt it = OSNSList.begin(); it != OSNSList.end(); ++it)
    {
      osnsXML = *it;
      type = osnsXML->getNSProblemType();
      if (type == LCP_TAG)  // LCP
      {
        (*_allNSProblems)[id].reset(new LCP(osnsXML));
      }
      else
        RuntimeException::selfThrow("EventDriven::xml constructor - wrong type of NSProblem: inexistant or not yet implemented");
      id = SICONOS_OSNSP_ED_SMOOTH_ACC;
    }
  }
  _newtonTolerance = DEFAULT_TOL_ED;
  _newtonMaxIteration = 50;
  _newtonNbSteps = 0;
  _newtonResiduDSMax = 0.0;
  _newtonResiduYMax = 0.0;
  _localizeEventMaxIter = 100;
}
/** defaut constructor
 *  \param a pointer to a timeDiscretisation (linked to the model that owns this simulation)
 */
EventDriven::EventDriven(SP::TimeDiscretisation td): Simulation(td), _istate(1)
{
  _numberOfOneStepNSproblems = 2;
  (*_allNSProblems).resize(_numberOfOneStepNSproblems);
  _newtonTolerance = DEFAULT_TOL_ED;
  _newtonMaxIteration = 50;
  _newtonNbSteps = 0;
  _newtonResiduDSMax = 0.0;
  _newtonResiduYMax = 0.0;
  _localizeEventMaxIter = 100;
}

EventDriven::EventDriven(SP::TimeDiscretisation td, int nb): Simulation(td), _istate(1)
{
  (*_allNSProblems).resize(nb);
  _newtonTolerance = DEFAULT_TOL_ED;
  _newtonMaxIteration = 50;
  _newtonNbSteps = 0;
  _newtonResiduDSMax = 0.0;
  _newtonResiduYMax = 0.0;
  _localizeEventMaxIter = 100;
}

double EventDriven::TOL_ED = DEFAULT_TOL_ED;




void EventDriven::insertIntegrator(SP::OneStepIntegrator osi)
{
  Simulation::insertIntegrator(osi);
  // Determine the number of OneStepNSproblems depending on the OneStepIntegrator type
  OSI::TYPES  osiType = osi->getType();
  if (osiType == OSI::NEWMARKALPHAOSI) // EventDrivent asscociated with NewMarkAlpha
  {
    _numberOfOneStepNSproblems = 3;
    if (_allNSProblems->size() != 3)
    {
      (*_allNSProblems).resize(_numberOfOneStepNSproblems);
    }
  }
}

void EventDriven::updateIndexSet(unsigned int i)
{
  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  assert(i < topo->indexSetsSize() &&
         "EventDriven::updateIndexSet(i), indexSets[i] does not exist.");
  // IndexSets[0] must not be updated in simulation, since it belongs to Topology.
  assert(i > 0  &&
         "EventDriven::updateIndexSet(i=0), indexSets[0] cannot be updated.");

  // For all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i].
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet1 = topo->indexSet(1);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);
  assert(indexSet0);
  assert(indexSet1);
  assert(indexSet2);

  // DEBUG_PRINTF("update indexSets start : indexSet0 size : %ld\n", indexSet0->size());
  // DEBUG_PRINTF("update IndexSets start : indexSet1 size : %ld\n", indexSet1->size());
  // DEBUG_PRINTF("update IndexSets start : indexSet2 size : %ld\n", indexSet2->size());

  InteractionsGraph::VIterator uibegin, uipend, uip;
  std11::tie(uibegin, uipend) = indexSet0->vertices();
  // loop over all vextice of the indexSet[i-1]
  for (uip = uibegin; uip != uipend; ++uip)
  {
    SP::Interaction inter = indexSet0->bundle(*uip);
    if (i == 1) // IndexSet[1]
    {
      // if indexSet[1]=>getYRef(0): output y
      // if indexSet[2]=>getYRef(1): output ydot
      double y = inter->getYRef(0); // output to define the IndexSets at this Interaction
      /*
         if (i == 1)
         {
         cout << "Id of Interaction: " << inter->number() << endl;
         cout << "Output of level 0 at this Interaction: " << y << endl;
         cout << endl;
         }
         */
      if (y < -TOL_ED) // y[0] < 0
      {
        RuntimeException::selfThrow("EventDriven::updateIndexSet, output of level 0 must be positive!!! ");
      }
      // 1 - If the Interaction is not yet in the set
      if (!indexSet1->is_vertex(inter)) // Interaction is not yet in the indexSet[i]
      {
        if (fabs(y) <= TOL_ED)
        {
          // vertex and edges insertions
          indexSet1->copy_vertex(inter, *indexSet0);
        }
      }
      else // if the Interaction was already in the set
      {
        if (fabs(y) > TOL_ED)
        {
          indexSet1->remove_vertex(inter); // remove the Interaction from IndexSet[1]
          inter->lambda(1)->zero(); // reset the lambda[1] to zero
        }
      }
    }
    else if (i == 2) // IndexSet[2]
    {
      if (indexSet1->is_vertex(inter)) // Interaction is in the indexSet[1]
      {
        double y = inter->getYRef(1); // output of level 1 at this Interaction
        if (!indexSet2->is_vertex(inter)) // Interaction is not yet in the indexSet[2]
        {
          if (fabs(y) <= TOL_ED)
          {
            // vertex and edges insertions
            indexSet2->copy_vertex(inter, *indexSet0);
          }
        }
        else // if the Interaction was already in the set
        {
          if (fabs(y) > TOL_ED)
          {
            indexSet2->remove_vertex(inter); // remove the Interaction from IndexSet[1]
            inter->lambda(2)->zero(); // reset the lambda[i] to zero
          }
        }
      }
      else // Interaction is not in the indexSet[1]
      {
        if (indexSet2->is_vertex(inter)) // Interaction is in the indexSet[2]
        {
          indexSet2->remove_vertex(inter); // remove the Interaction from IndexSet[2]
          inter->lambda(2)->zero(); // reset the lambda[i] to zero
        }
      }
    }
    else
    {
      RuntimeException::selfThrow("EventDriven::updateIndexSet, IndexSet[i > 2] doesn't exist");
    }
  }

  // DEBUG_PRINTF("update indexSets end : indexSet0 size : %ld\n", indexSet0->size());
  // DEBUG_PRINTF("update IndexSets end : indexSet1 size : %ld\n", indexSet1->size());
  // DEBUG_PRINTF("update IndexSets end : indexSet2 size : %ld\n", indexSet2->size());
}

void EventDriven::updateIndexSetsWithDoubleCondition()
{

  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());

  // for all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]

  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  SP::InteractionsGraph indexSet2 = topo->indexSet(2);

  InteractionsGraph::VIterator ui, uiend, vnext;
  std11::tie(ui, uiend) = indexSet2->vertices();

  for (vnext = ui; ui != uiend; ui = vnext)
  {
    ++vnext;

    SP::Interaction inter = indexSet2->bundle(*ui);
    double gamma = inter->getYRef(2);
    double F     = inter->getLambdaRef(2);
    if (fabs(F) < TOL_ED)
      indexSet2->remove_vertex(inter);
    else if ((gamma < -TOL_ED) || (F < -TOL_ED))
      RuntimeException::selfThrow("EventDriven::updateIndexSetsWithDoubleCondition(), output[2] and lambda[2] for Interactionof indexSet[2] must be higher or equal to zero.");
    else if (((fabs(gamma) > TOL_ED) && (fabs(F) > TOL_ED)))
      RuntimeException::selfThrow("EventDriven::updateIndexSetsWithDoubleCondition(), something is wrong for the LCP resolution.");
  }
}

void EventDriven::initializeInteraction(SP::Interaction inter)
{

  RELATION::TYPES pbType = inter->getRelationType();
  if (pbType == Lagrangian)
  {
    //    inter->setDataXFromVelocity();
  }
  else
    RuntimeException::selfThrow("EventDriven::initializeInteractions(SP::interaction inter) - not implemented for Relation of type " + pbType);

}

void EventDriven::initOSNS()
{
  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());
  // === initialization for OneStepIntegrators ===
  OSI::TYPES  osiType = (*_allOSI->begin())->getType();
  for (OSIIterator itosi = _allOSI->begin();  itosi != _allOSI->end(); ++itosi)
  {
    //Check whether OSIs used are of the same type
    if ((*itosi)->getType() != osiType)
      RuntimeException::selfThrow("OSIs used must be of the same type");
    for (DSIterator itds = (*itosi)->dynamicalSystems()->begin();
         itds != (*itosi)->dynamicalSystems()->end(); ++itds)
    {
      // Initialize right-hand side
      (*itds)->initRhs(model()->t0());
    }
  }
  // for all Interactions in indexSet[i-1], compute y[i-1] and
  // update the indexSet[i]
  InteractionsGraph::VIterator ui, uiend;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();

  SP::InteractionsGraph indexSet0 = topo->indexSet(0);

  // For each Interaction in I0 ...
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    // indexSet0->bundle(*ui)->initialize("EventDriven");
    initializeInteraction(indexSet0->bundle(*ui));
  }

  // === update all index sets ===
  updateIndexSets();

  if (!_allNSProblems->empty()) // ie if some Interactions have been
    // declared and a Non smooth problem built.
  {
    if (osiType == OSI::LSODAR) //EventDriven associated with Lsodar OSI
    {
      // === OneStepNSProblem initialization. === First check that
      // there are 2 osns: one "impact" and one "acceleration"
      // if(_allNSProblems->size()!=2)
      //   RuntimeException::selfThrow
      //     (" EventDriven::initialize, \n an EventDriven simulation associated with Lsodar must have two non smooth problem.\n Here, there are "
      //      +_allNSProblems->size());
    }
    else if (osiType == OSI::NEWMARKALPHAOSI) // EventDrivent asscociated with NewMarkAlpha
    {
      if (_allNSProblems->size() != 3)
        RuntimeException::selfThrow
        (" EventDriven::initialize, \n an EventDriven simulation associated with NewMarkAlphaOSI must have three non smooth problem.\n Here, there are "
         + _allNSProblems->size());
      // Initialize OSNSP at position level
      (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_POS]->setLevels(2, 2);
      (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_POS]->initialize(shared_from_this());
    }
    else
    {
      RuntimeException::selfThrow(" EventDriven::initialize, this OneStepIntegrator has not implemented yet.");
    }

    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_IMPACT])) // ie if the impact problem does not
      // exist
      RuntimeException::selfThrow
      ("EventDriven::initialize, an EventDriven simulation must have an 'impact' non smooth problem.");
    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC])) // ie if the acceleration-level problem
      // does not exist
      RuntimeException::selfThrow
      ("EventDriven::initialize, an EventDriven simulation must have an 'acceleration' non smooth problem.");
    // Initialize OSNSP for impact problem and at the acceleration level
    // WARNING: only for Lagrangian systems - To be reviewed for other ones.
    (*_allNSProblems)[SICONOS_OSNSP_ED_IMPACT]->setLevels(1, 1);
    (*_allNSProblems)[SICONOS_OSNSP_ED_IMPACT]->initialize(shared_from_this());
    (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->setLevels(2, 2);
    (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->initialize(shared_from_this());
    //
    // Detect NonSmoothEvent at the beginning of the simulation
    SP::InteractionsGraph indexSet1 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    if (indexSet1->size() != 0) // There is one non-smooth event to be added
    {
      _eventsManager->scheduleNonSmoothEvent(_eventsManager->startingTime(), false);
    };
  }
}

void EventDriven::initOSIs()
{
  for (OSIIterator itosi = _allOSI->begin();  itosi != _allOSI->end(); ++itosi)
  {
    // Initialize the acceleration like for NewMarkAlphaScheme
    if ((*itosi)->getType() == OSI::NEWMARKALPHAOSI)
    {
      SP::NewMarkAlphaOSI osi_NewMark =  std11::static_pointer_cast<NewMarkAlphaOSI>(*itosi);
      for (DSIterator itds = (*itosi)->dynamicalSystems()->begin();
           itds != (*itosi)->dynamicalSystems()->end(); ++itds)
      {
        if ((Type::value(**itds) == Type::LagrangianDS) || (Type::value(**itds) == Type::LagrangianLinearTIDS))
        {
          SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS>(*itds);
          *(d->getWorkVector(DynamicalSystem::acce_like)) = *(d->acceleration()); // set a0 = ddotq0
          // Allocate the memory to stock coefficients of the polynomial for the dense output
          d->allocateWorkMatrix(LagrangianDS::coeffs_denseoutput, (*itds)->getDim(), (osi_NewMark->getOrderDenseOutput() + 1));
        }
      }
    }
  }
}


void EventDriven::initialize(SP::Model m, bool withOSI)
{
  // Initialization for Simulation
  Simulation::initialize(m, withOSI);
  // Initialization for all OneStepIntegrators
  initOSIs();
}



// void EventDriven::initLevelMin()
// {
//   // At the time, we consider that for all systems, levelMin is
//   // equal to the minimum value of the relative degree
//   _levelMin = model()->nonSmoothDynamicalSystem()
//               ->topology()->minRelativeDegree();
//   if(_levelMin==0)
//     _levelMin++;
// }



// void EventDriven::initLevelMax()
// {
//   _levelMax = model()->nonSmoothDynamicalSystem()->topology()->maxRelativeDegree();
//   // Interactions initialization (here, since level depends on the
//   // type of simulation) level corresponds to the number of Y and
//   // Lambda derivatives computed.
//   if(_levelMax==0)
//     _levelMax++;
// }

void EventDriven::computef(SP::OneStepIntegrator osi, integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{
  //std::cout << "EventDriven::computef -------------------------> start" <<std::endl;

  // computeF is supposed to fill xdot in, using the definition of the
  // dynamical systems belonging to the osi

  // Check osi type: only lsodar is allowed.
  if (osi->getType() != OSI::LSODAR)
    RuntimeException::selfThrow("EventDriven::computef(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  Lsodar& lsodar = static_cast<Lsodar&>(*osi);
  // fill in xWork vector (ie all the x of the ds of this osi) with x
  lsodar.fillXWork(sizeOfX, x);

  double t = *time;
  model()->setCurrentTime(t);
  // solve a LCP at "acceleration" level if required
  if (!_allNSProblems->empty())
  {
    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->interactions())->isEmpty())
    {
      // Update the state of the DS
      (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(t);
      updateInput(2); // Necessary to compute DS state below
    }
    // Compute the right-hand side ( xdot = f + r in DS) for all the
    //ds, with the new value of input.  lsodar->computeRhs(t);
  }
  // update the DS of the OSI.
  lsodar.computeRhs(t);
  //  for the DS state, ie the ones computed by lsodar (x above)
  // Update Index sets? No !!

  // Get the required value, ie xdot for output.
  DSIterator it;
  unsigned int i = 0;
  for (it = lsodar.dynamicalSystemsBegin(); it != lsodar.dynamicalSystemsEnd(); ++it)
  {
    if (Type::value(**it) == Type::LagrangianDS ||
        Type::value(**it) == Type::LagrangianLinearTIDS)
    {
      LagrangianDS& LDS = *std11::static_pointer_cast<LagrangianDS>(*it);
      SiconosVector& qDotTmp = *LDS.velocity();
      SiconosVector& qDotDotTmp = *LDS.acceleration();
      for (unsigned int j = 0 ; j < (*it)->getDim() ; ++j)
        xdot[i++] = qDotTmp(j);
      for (unsigned int j = 0 ; j < (*it)->getDim() ; ++j)
        xdot[i++] = qDotDotTmp(j);
    }
    else
    {
      SiconosVector& xtmp2 = *(*it)->rhs(); // Pointer link !
      for (unsigned int j = 0 ; j < (*it)->getN() ; ++j) // Warning: getN, not getDim !!!!
        xdot[i++] = xtmp2(j);
    }
  }

  //std::cout << "EventDriven::computef -------------------------> stop" <<std::endl;

}

void EventDriven::computeJacobianfx(SP::OneStepIntegrator osi,
                                    integer *sizeOfX,
                                    doublereal *time,
                                    doublereal *x,
                                    doublereal *jacob)
{
  if (osi->getType() != OSI::LSODAR)
    RuntimeException::selfThrow("EventDriven::computeJacobianfx(osi, ...), not yet implemented for a one step integrator of type " + osi->getType());

  SP::Lsodar lsodar = std11::static_pointer_cast<Lsodar>(osi);

  // Remark A: according to DLSODAR doc, each call to jacobian is
  // preceded by a call to f with the same arguments NEQ, T, and Y.
  // Thus to gain some efficiency, intermediate quantities shared by
  // both calculations may be saved in class members?  fill in xWork
  // vector (ie all the x of the ds of this osi) with x fillXWork(x);
  // -> copy
  // Maybe this step is not necessary?  because of
  // remark A above

  // Compute the jacobian of the vector field according to x for the
  // current ds
  double t = *time;
  model()->setCurrentTime(t);
  lsodar->computeJacobianRhs(t);

  // Save jacobianX values from dynamical system into current jacob
  // (in-out parameter)

  DSIterator it;
  unsigned int i = 0;
  for (it = lsodar->dynamicalSystemsBegin(); it != lsodar->dynamicalSystemsEnd(); ++it)
  {
    if (Type::value(**it) == Type::LagrangianDS ||
        Type::value(**it) == Type::LagrangianLinearTIDS)
    {
      LagrangianDS& lds = *std11::static_pointer_cast<LagrangianDS>(*it);
      BlockMatrix& jacotmp = *lds.jacobianRhsx();
      for (unsigned int j = 0; j < (*it)->getN(); ++j)
      {
        for (unsigned int k = 0; k < (*it)->getDim(); ++k)
          jacob[i++] = jacotmp(k, j);
      }
    }
    else
    {
      SiconosMatrix& jacotmp = *(*it)->jacobianRhsx(); // Pointer link !
      for (unsigned int j = 0; j < (*it)->getN(); ++j)
      {
        for (unsigned int k = 0; k < (*it)->getDim(); ++k)
          jacob[i++] = jacotmp(k, j);
      }
    }
  }
}

unsigned int EventDriven::computeSizeOfg()
{
  SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  return (indexSet0->size());
}


void EventDriven::computeg(SP::OneStepIntegrator osi,
                           integer * sizeOfX, doublereal* time,
                           doublereal* x, integer * ng,
                           doublereal * gOut)
{
  assert(!_model.expired());
  assert(model()->nonSmoothDynamicalSystem());
  assert(model()->nonSmoothDynamicalSystem()->topology());
  InteractionsGraph::VIterator ui, uiend;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);
  unsigned int nsLawSize, k = 0 ;
  SP::SiconosVector y, ydot, lambda;
  SP::Lsodar lsodar = std11::static_pointer_cast<Lsodar>(osi);
  // Solve LCP at acceleration level to calculate the lambda[2] at Interaction of indexSet[2]
  lsodar->fillXWork(sizeOfX, x);
  //
  double t = *time;
  model()->setCurrentTime(t);
  if (!_allNSProblems->empty())
  {
    if (!((*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->interactions())->isEmpty())
    {
      (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(t);
    }
  };
  /*
     double * xdottmp = (double *)malloc(*sizeOfX*sizeof(double));
     computef(osi, sizeOfX,time,x,xdottmp);
     free(xdottmp);
     */
  // Update the output from level 0 to level 1
  updateOutput(0);
  updateOutput(1);
  //
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui);
    nsLawSize = inter->getNonSmoothLawSize();
    y = inter->y(0);   // output y at this Interaction
    ydot = inter->y(1); // output of level 1 at this Interaction
    lambda = inter->lambda(2); // input of level 2 at this Interaction
    if (!(indexSet2->is_vertex(inter))) // if Interaction is not in the indexSet[2]
    {
      for (unsigned int i = 0; i < nsLawSize; ++i)
      {
        if ((*y)(i) >= TOL_ED) // y[0] > 0
        {
          gOut[k] = (*y)(i);
        }
        else // y[0] = 0
        {
          if ((*ydot)(i) >= 0) // if y[1] >= 0;
          {
            gOut[k] = TOL_ED; // g = TOL_ED temporarily
          }
          else  // if y[1] < 0
          {
            gOut[k] = (*y)(i); // g = y[0]
          }
        }
        k++;
      }
    }
    else // If Interaction is in the indexSet[2]
    {
      for (unsigned int i = 0; i < nsLawSize; ++i)
      {
        gOut[k] = (*lambda)(i); // g = lambda[2]
        k++;
      }
    }
  }
}
void EventDriven::updateImpactState()
{
  OSIIterator itOSI;
  // Compute input = R(lambda[1])
  updateInput(1);

  // Compute post-impact velocity
  for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    (*itOSI)->updateState(1);
}

void EventDriven::updateSmoothState()
{
  // Update input of level 2
  updateInput(2);
  OSIIterator itOSI;
  // Compute acceleration
  for (itOSI = _allOSI->begin(); itOSI != _allOSI->end() ; ++itOSI)
    (*itOSI)->updateState(2);
}


void EventDriven::update(unsigned int levelInput)
{
  assert(levelInput <= 2);
  if (levelInput == 1)
  {
    updateImpactState();
  }
  else
  {
    updateSmoothState();
  }
  // Update output (y)
  updateOutput(levelInput);
  // Warning: index sets are not updated in this function !!
}

void EventDriven::advanceToEvent()
{
  _tinit = _eventsManager->startingTime();
  _tend =  _eventsManager->nextTime();
  _tout = _tend;
  bool isNewEventOccur = false;  // set to true if a new event occur during integration
  OSI::TYPES  osiType = (*_allOSI->begin())->getType(); // Type of OSIs
  double _maxConstraint = 0.0;
  if (osiType == OSI::NEWMARKALPHAOSI)
  {
    newtonSolve(_newtonTolerance, _newtonMaxIteration);
    // Update after Newton iteration
    // Update input of level 2 >>> has already been done in newtonSolve
    // Update state of all Dynamicall Systems >>>  has already been done in newtonSolve
    // Update outputs of levels 0, 1, 2
    updateOutput(0);
    updateOutput(1);
    updateOutput(2);
    // Detect whether or not some events occur during the integration step
    _maxConstraint = detectEvents();
    LocalizeFirstEvent();
    //
#ifdef DEBUG_MESSAGES
    cout << "========== EventDriven::advanceToEvent =============" << endl;
    cout << "Istate: " << _istate << endl;
    cout << "Maximum value of constraint functions: " << _maxConstraint << endl;
#endif
    //
    if (_istate != 2) //some events occur
    {
      cout << "In EventDriven::advanceToEvent, some events are detected!!!" << endl;
      if (_maxConstraint < TOL_ED) // events occur at the end of the integration step
      {
        isNewEventOccur = true;
      }
      else // events need to be localized
      {
        isNewEventOccur = true;
        // LocalizeFirstEvent();
      }
    }
  }
  else if (osiType == OSI::LSODAR)
  {
    // WARNING: this is supposed to work for only one OSI, including all
    // the DS.  To be reviewed for multiple OSI case (if it has sense?).

    // ---> Step 1: integrate the smooth dynamics from current event to
    // next event; Starting event = last accessed event.  Next event =
    // next time step or first root of the 'g' function found by
    // integrator (Lsodar)

    // if _istate == 1 => first call. It this case we suppose that _tinit
    // and _tend have been initialized before
    // if(_istate == 2 || _istate == 3)
    //  {
    //    _tinit = _eventsManager->startingTime();
    //    _tend =  _eventsManager->nextTime();
    //  }
    // _tout = _tend;
    // call integrate method for each OSI, between _tinit and _tend.
    OSIIterator it;
    for (it = _allOSI->begin(); it != _allOSI->end(); ++it)
    {

      (*it)->resetNonSmoothPart();
      //====================================================================================
      //    std::cout << " Start of Lsodar integration" << std::endl;
      (*it)->integrate(_tinit, _tend, _tout, _istate); // integrate must
      //    std::cout << " End of Lsodar integration" << std::endl;
      // return a flag (_istate) telling if _tend has been  reached or not.
      //====================================================================================

      if (_printStat)
      {
        statOut << " =================> Results after advanceToEvent <================= " << endl;
        statOut << " Starting time: " << _tinit << endl;
        statOut << " _istate " << _istate << endl;
      }
      if (_istate == 3) // ie if _tout is not equal to _tend: one or more roots have been found.
      {
        isNewEventOccur = true;
        // Add an event into the events manager list
        _eventsManager->scheduleNonSmoothEvent(_tout);
        if (_printStat)
          statOut << " -----------> New non-smooth event at time " << _tout << endl;
      }
      if (_printStat)
      {
        SP::Lsodar lsodar = std11::static_pointer_cast<Lsodar>(*it);
        statOut << "Results at time " << _tout << ":" << endl;
        SA::integer iwork = lsodar->getIwork();
        SA::doublereal Rwork = lsodar->getRwork();
        statOut << "Number of steps: " << iwork[10] << ", number of f evaluations: " << iwork[11] << ", number of jacobianF eval.: " << iwork[12] << "." << endl;
      }
    }
    // Set model time to _tout
    model()->setCurrentTime(_tout);
    //update output[0], output[1]
    updateOutput(0);
    updateOutput(1);
    // Update all the index sets ...
    updateIndexSets();
    //update lambda[2], input[2] and indexSet[2] with double consitions for the case there is no new event added during time integration, otherwise, this
    // update is done when the new event is processed
    if (!isNewEventOccur)
    {
      if (!_allNSProblems->empty())
      {
        // Solve LCP at acceleration level
        if (!((*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->interactions())->isEmpty())
        {
          (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(_tout);
          updateInput(2); //
          //  updateInput(1); // this is done to reset the nonsmoothinput at the level of impact
        }
        // update indexSet[2] with double condition
        updateIndexSetsWithDoubleCondition();
      }
    }
  }
  else
  {
    RuntimeException::selfThrow("In EventDriven::advanceToEvent, this type of OneStepIntegrator does not exist for Event-Driven scheme!!!");
  }
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double EventDriven::computeResiduGaps()
{
  // Make sure that the state of all Dynamical Systems was updated
  double t = nextTime(); // time at the end of the step
  SP::InteractionsGraph indexSet2 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(2);
  SP::SiconosVector _y;
  // Loop over all interactions of indexSet2
  InteractionsGraph::VIterator ui, uiend;
  double _maxResiduGap = 0.0;
  for (std11::tie(ui, uiend) = indexSet2->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet2->bundle(*ui);
    inter->computeOutput(t, 0); // compute y[0] for the interaction at the end time
    _y = inter->y(0);
    if (_maxResiduGap < abs((*_y)(0))) // (*_y)[0] gives gap at this interaction
    {
      _maxResiduGap = abs((*_y)(0));
    }
    //
#ifdef DEBUG_MESSAGES
    cout << "gap at contact: ";
    _y->display();
#endif
    //
  }
  //
#ifdef DEBUG_MESSAGES
  cout << "Maximum gap residu: " << _maxResiduGap << endl;
#endif
  //
  return _maxResiduGap;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EventDriven::prepareNewtonIteration()
{
  // At this stage, we do
  // (1) compute iteration matrix W for all DSs belonging to all OSIs
  // (2) compute free residu for all DSs belonging to all OSIs and get maximum residu
  // (3) compute free state for all DSs belonging to all OSIs
  // (4) compute maximum gap residu over all interactions of indexSet 2
  _newtonResiduDSMax = 0.0;
  _newtonResiduYMax = 0.0;
  double _maxResidu;
  // Update input of level 2
  updateInput(2);
  // Loop over all OSIs
  OSI::TYPES  osiType;
  for (OSIIterator itosi = _allOSI->begin(); itosi != _allOSI->end(); ++itosi)
  {
    osiType = (*itosi)->getType();
    if (osiType != OSI::NEWMARKALPHAOSI)
    {
      RuntimeException::selfThrow("In EventDriven::prepareNewtonIteration, the current OSI is not NewMarkAlpha scheme!!!");
    }
    // Compute iteration matrix W
    (*itosi)->prepareNewtonIteration(nextTime()); // preparation for each OSI
    // Compute free residus, maximum residu
    _maxResidu = (*itosi)->computeResidu();
    if (_newtonResiduDSMax < _maxResidu)
    {
      _newtonResiduDSMax = _maxResidu;
    }
    // Compute free state
    (*itosi)->computeFreeState();
  }
  // Compute maximum gap residu
  _newtonResiduYMax = computeResiduGaps();
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
bool EventDriven::newtonCheckConvergence(double _tol)
{
  bool checkConvergence = true;
  if ((_newtonResiduDSMax > _tol) || (_newtonResiduYMax > _tol))
  {
    checkConvergence = false;
  }
  return checkConvergence;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EventDriven::predictionNewtonIteration()
{
  // Prediction of the state for all Dynamical Systems before Newton iteration
  for (OSIIterator itosi = _allOSI->begin(); itosi != _allOSI->end(); ++itosi)
  {
    if ((*itosi)->getType() == OSI::NEWMARKALPHAOSI)
    {
      SP::NewMarkAlphaOSI osi_NewMark = std11::static_pointer_cast<NewMarkAlphaOSI>(*itosi);
      osi_NewMark->prediction();
    }
    else
    {
      RuntimeException::selfThrow("In EventDriven::predictionNewtonIteration, the current OSI must be NewMarkAlpha scheme!!!");
    }
  }
  // Prediction of the output and lambda for all Interactions before Newton iteration
  double t = nextTime(); // time at the end of the step
  SP::InteractionsGraph indexSet0 = model()->nonSmoothDynamicalSystem()->topology()->indexSet(0);
  // Loop over all interactions
  InteractionsGraph::VIterator ui, uiend;
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui);
    inter->computeOutput(t, 0); // compute y[0] for the interaction at the end time with the state predicted for Dynamical Systems
    inter->lambda(2)->zero(); // reset lambda[2] to zero
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EventDriven::correctionNewtonIteration()
{
  //Update the input of level 2 for all Dynamical Systems after each iteration
  updateInput(2);
  // Correction
  for (OSIIterator itosi = _allOSI->begin(); itosi != _allOSI->end(); ++itosi)
  {
    if ((*itosi)->getType() == OSI::NEWMARKALPHAOSI)
    {
      SP::NewMarkAlphaOSI osi_NewMark = std11::static_pointer_cast<NewMarkAlphaOSI>(*itosi);
      osi_NewMark->correction();
    }
    else
    {
      RuntimeException::selfThrow("In EventDriven::correctionNewtonIteration, the current OSI must be NewMarkAlpha scheme!!!");
    }
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EventDriven::newtonSolve(double criterion, unsigned int maxStep)
{
  _isNewtonConverge = false;
  _newtonNbSteps = 0; // number of Newton iterations
  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
  int info = 0;
  _istate = 1; // beginning of time integration
  // Prediction
  predictionNewtonIteration();
  while (1 != 0)
  {
    _newtonNbSteps++;
    // Prepare for iteration
    prepareNewtonIteration();
    // Check convergence
    _isNewtonConverge = newtonCheckConvergence(_newtonTolerance);
    //
#ifdef DEBUG_MESSAGES
    cout << "Iteration: " << _newtonNbSteps << endl;
    cout << "Convergence: " << _isNewtonConverge << endl;
#endif
    //
    if (_isNewtonConverge)
    {
      break;
    }
    if (_newtonNbSteps >  maxStep)
    {
      cout << "Warning!!!In EventDriven::newtonSolve: Number of iterations is greater than the maximum value " << maxStep << endl;
    }
    // If no convergence, solve LCP
    if (!_allNSProblems->empty() && !allInteractions->isEmpty())
    {
      info = computeOneStepNSProblem(SICONOS_OSNSP_ED_SMOOTH_POS);
      if (info != 0)
      {
        cout << "Warning!!!In EventDriven::newtonSolve: LCP solver may fail" << endl;
      }
    }
    // Correction of the state of all Dynamical Systems
    correctionNewtonIteration();
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
double EventDriven::detectEvents()
{
  double _maxResiduOutput = 0.0; // maximum of abs(g_i) with i running over all activated or deactivated contacts
  // Loop over all interactions to detect whether some constraints are activated or deactivated
  bool _IsContactClosed = false;
  bool _IsContactOpened = false;
  InteractionsGraph::VIterator ui, uiend;
  SP::SiconosVector y, ydot, lambda;
  SP::Topology topo = model()->nonSmoothDynamicalSystem()->topology();
  SP::InteractionsGraph indexSet0 = topo->indexSet(0);
  SP::InteractionsGraph indexSet2 = topo->indexSet(2);
  //
#ifdef DEBUG_MESSAGES
  cout << "======== In EventDriven::detectEvents =========" << endl;
#endif
  for (std11::tie(ui, uiend) = indexSet0->vertices(); ui != uiend; ++ui)
  {
    SP::Interaction inter = indexSet0->bundle(*ui);
    double nsLawSize = inter->getNonSmoothLawSize();
    if (nsLawSize != 1)
    {
      RuntimeException::selfThrow("In EventDriven::detectEvents, the interaction size > 1 has not been implemented yet!!!");
    }
    y = inter->y(0);   // output y at this Interaction
    ydot = inter->y(1); // output of level 1 at this Interaction
    lambda = inter->lambda(2); // input of level 2 at this Interaction
    if (!(indexSet2->is_vertex(inter))) // if Interaction is not in the indexSet[2]
    {
      if ((*y)(0) < TOL_ED) // gap at the current interaction <= 0
      {
        _IsContactClosed = true;
        if (_maxResiduOutput < std::abs((*y)(0)))
        {
          _maxResiduOutput = std::abs((*y)(0));
        }
      }
    }
    else // If interaction is in the indexSet[2]
    {
      if ((*lambda)(0) < TOL_ED) // normal force at the current interaction <= 0
      {
        _IsContactOpened = true;
        if (_maxResiduOutput < std::abs((*lambda)(0)))
        {
          _maxResiduOutput = std::abs((*lambda)(0));
        }
      }
    }
    //
#ifdef DEBUG_MESSAGES
    cout << "Contact number: " << inter->number() << endl;
    cout << "Contact gap: " << (*y)(0) << endl;
    cout << "Contact force: " << (*lambda)(0) << endl;
    cout << "Is contact is closed: " << _IsContactClosed << endl;
    cout << "Is contact is opened: " << _IsContactOpened << endl;
#endif
    //
  }
  //
  if ((!_IsContactClosed) && (!_IsContactOpened))
  {
    _istate = 2; //no event is detected
  }
  else if ((_IsContactClosed) && (!_IsContactOpened))
  {
    _istate = 3; // Only some contacts are closed
  }
  else if ((!_IsContactClosed) && (_IsContactOpened))
  {
    _istate = 4; // Only some contacts are opened
  }
  else
  {
    _istate = 5; // Some contacts are closed AND some contacts are opened
  }
  //
  return  _maxResiduOutput;
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
void EventDriven::LocalizeFirstEvent()
{
  // We localize the first event occuring during the integration step when the flag _istate = 3 or 4
  // Compute the coefficients of the dense output polynomial for all DSs
  for (OSIIterator itosi = _allOSI->begin(); itosi != _allOSI->end(); ++itosi)
  {
    if ((*itosi)->getType() == OSI::NEWMARKALPHAOSI)
    {
      SP::NewMarkAlphaOSI osi_NewMark = std11::static_pointer_cast<NewMarkAlphaOSI>(*itosi);
      osi_NewMark->prepareEventLocalization();
    }
    else
    {
      RuntimeException::selfThrow("In EventDriven::LocalizeFirstEvent, the current OSI must be NewMarkAlpha scheme!!!");
    }
  }
  //
  /*
  double t_a = startingTime();
  double t_b = nextTime();
  double _maxConstraint = 0.0;
  bool found = false;
  unsigned int _numIter = 0;
  while (!found)
    {
      _numIter++;
      double t_i = (t_b - t_a)/2.0; // mid-time of the current interval
      // set t_i as the current time
      model()->setCurrentTime(t_i);
      // Generate dense output for all DSs at the time t_i
      for(OSIIterator itosi = _allOSI->begin(); itosi != _allOSI->end(); ++itosi)
  {
    if ((*itosi)->getType() == OSI::NEWMARKALPHAOSI)
      {
        SP::NewMarkAlphaOSI osi_NewMark = std11::static_pointer_cast<NewMarkAlphaOSI>(*itosi);
        osi_NewMark->DenseOutputallDSs(t_i);
      }
  }
      // If _istate = 3 or 5, i.e. some contacts are closed, we need to compute y[0] for all interactions
      if ((_istate == 3)||(_istate == 5)) // some contacts are closed
  {
    updateOutput(0);
  }
      // If _istate = 4 or 5, i.e. some contacts are detached, we need to solve LCP at the acceleration level to compute contact forces
      if ((_istate == 4)||(_istate == 5)) // some contacts are opened
  {
    if(!_allNSProblems->empty())
      {
        (*_allNSProblems)[SICONOS_OSNSP_ED_SMOOTH_ACC]->compute(t_i);
      }
  }
      // Check whether or not some events occur in the interval [t_a, t_i]
       _maxConstraint = detectEvents();
       if ((_istate != 2)&&(_maxConstraint < TOL_ED)) // first event is found
   {
     _tout = t_i;
     found = true;
   }
      // if some events are detected in the interval [t_a, t_i] (if _istate != 2), set t_b = t_i
      if (_istate != 2)
  {
    t_b = t_i;
  }
      else // if no event is detected in [t_a, t_i], then we have to detect events in the interval [t_i, t_b]
  {
    t_a = t_i;
  }
      //
      if (_numIter > _localizeEventMaxIter)
  {
    RuntimeException::selfThrow("In EventDriven::LocalizeFirstEvent, the numbner of iterations performed is too large!!!");
  }
    }
  */
}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
EventDriven* EventDriven::convert(Simulation *str)
{
  EventDriven* ed = dynamic_cast<EventDriven*>(str);
  return ed;
}

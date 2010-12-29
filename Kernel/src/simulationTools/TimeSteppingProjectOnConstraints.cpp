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

#include "TimeStepping.hpp"

#include "TimeSteppingProjectOnConstraints.hpp"
#include "NewtonEulerDS.hpp"

using namespace std;


TimeSteppingProjectOnConstraints::TimeSteppingProjectOnConstraints(SP::TimeDiscretisation td,
    SP::OneStepIntegrator osi,
    SP::OneStepNSProblem osnspb_velo,
    SP::OneStepNSProblem osnspb_pos)
  : TimeStepping(td, osi, osnspb_velo)
{
  (*_allNSProblems).resize(SICONOS_NB_OSNSP_TSP);
  insertNonSmoothProblem(osnspb_pos, SICONOS_OSNSP_TS_POS);
}

// --- Destructor ---
TimeSteppingProjectOnConstraints::~TimeSteppingProjectOnConstraints()
{
}


void TimeSteppingProjectOnConstraints::newtonSolve(double criterion, unsigned int maxStep)
{
  TimeStepping::newtonSolve(criterion, maxStep);
  int info = 0;
  //cout<<"begin projection:\n";
  SP::InteractionsSet allInteractions = model()->nonSmoothDynamicalSystem()->interactions();
  for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  {
    (*it)->relation()->computeJach(getTkp1());
    (*it)->relation()->computeh(getTkp1());
  }

  info = computeOneStepNSProblem(SICONOS_OSNSP_TS_POS);
  //cout<<"||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||  Z:"<<endl;
  //(*_allNSProblems)[SICONOS_OSNSP_TS_POS]->display();
  //(boost::static_pointer_cast<LinearOSNS>((*_allNSProblems)[SICONOS_OSNSP_TS_POS]))->z()->display();
  if (info)
    cout << "TimeSteppingProjectOnConstraints project on constraints failed." << endl ;
  //cout<<"during projection before normalizing of q:\n";
  //for (InteractionsIterator it = allInteractions->begin(); it != allInteractions->end(); it++)
  //{
  //  (*it)->relation()->computeh(getTkp1());
  //}
  SP::DynamicalSystemsGraph dsGraph = model()->nonSmoothDynamicalSystem()->dynamicalSystems();

  for (DynamicalSystemsGraph::VIterator vi = dsGraph->begin(); vi != dsGraph->end(); ++vi)
  {
    SP::DynamicalSystem ds = dsGraph->bundle(*vi);
    Type::Siconos dsType = Type::value(*ds);
    if (dsType != Type::NewtonEulerDS)
      RuntimeException::selfThrow("TS:: - ds is not from NewtonEulerDS.");
    SP::NewtonEulerDS neds = boost::static_pointer_cast<NewtonEulerDS>(ds);
    neds->normalizeq();
    // SP::SiconosVector dotq = neds->dotq();
    // SP::SiconosVector q = neds->q();
    // SP::SiconosVector qold = neds->qMemory()->getSiconosVector(0);
    // double h = timeStep();
    // dotq->setValue(3,(q->getValue(3)-qold->getValue(3))/h);
    // dotq->setValue(4,(q->getValue(4)-qold->getValue(4))/h);
    // dotq->setValue(5,(q->getValue(5)-qold->getValue(5))/h);
    // dotq->setValue(6,(q->getValue(6)-qold->getValue(6))/h);
    neds->updateT();
  }

}
bool TimeSteppingProjectOnConstraints::predictorDeactivate(SP::UnitaryRelation ur, unsigned int i)
{
  double h = timeStep();
  double y = ur->getYRef(i - 1);
  double yDot = ur->getYRef(1);
  y += 0.5 * h * yDot;
  assert(!isnan(y));
  return (y > 10e-4);
}


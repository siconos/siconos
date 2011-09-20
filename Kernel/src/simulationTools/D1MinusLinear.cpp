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

#include "D1MinusLinear.hpp"
#include "Simulation.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"

using namespace std;
using namespace RELATION;

struct D1MinusLinear::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  OneStepNSProblem *osnsp;
  SP::UnitaryRelation UR;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::UnitaryRelation UR) : osnsp(p), UR(UR) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = UR->getNonSmoothLawSize();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    subscal(e, *(UR->y_k(osnsp->levelMin())), *(UR->yp()), subCoord, false);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e = nslaw.en();
    (*(UR->yp()))(0) +=  e * (*(UR->y_k(osnsp->levelMin())))(0);
  }

  void visit(const EqualityConditionNSL& nslaw) {}

  void visit(const MixedComplementarityConditionNSL& nslaw) {}
};

// --- constructor from a ds ---
D1MinusLinear::D1MinusLinear(SP::DynamicalSystem newDS) :
  OneStepIntegrator(OSI::D1MINUSLINEAR)
{
  OSIDynamicalSystems->insert(newDS);
}

// --- constructor from a list of ds ---
D1MinusLinear::D1MinusLinear(DynamicalSystemsSet& newDS): OneStepIntegrator(OSI::D1MINUSLINEAR, newDS) {}

const SimpleMatrix D1MinusLinear::getW(SP::DynamicalSystem ds)
{
  assert(ds && "D1MinusLinear::getW(ds) - ds == NULL.");
  assert(WMap[ds] && "D1MinusLinear::getW(ds) - W[ds] == NULL.");
  return *(WMap[ds]); // copy
}

SP::SimpleMatrix D1MinusLinear::W(SP::DynamicalSystem ds)
{
  assert(ds && "D1MinusLinear::W(ds) - ds == NULL.");
  return WMap[ds];
}

void D1MinusLinear::setW(const SiconosMatrix& newValue, SP::DynamicalSystem ds)
{
  // check if ds is in the OSI
  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("D1MinusLinear::setW(newVal,ds) - ds does not belong to this Integrator...");

  // check dimensions consistency
  unsigned int line = newValue.size(0);
  unsigned int col  = newValue.size(1);

  if (line != col)
    RuntimeException::selfThrow("D1MinusLinear::setW(newVal,ds) - newVal is not square.");

  if (!ds)
    RuntimeException::selfThrow("D1MinusLinear::setW(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim();
  if (line != sizeW)
    RuntimeException::selfThrow("D1MinusLinear::setW(newVal,ds) - inconsistent dimension between newVal and dynamical system to be integrated.");

  // memory allocation for W, if required
  if (!WMap[ds])
    WMap[ds].reset(new SimpleMatrix(newValue));
  else if (line == WMap[ds]->size(0) && col == WMap[ds]->size(1))
    *(WMap[ds]) = newValue;
  else
    RuntimeException::selfThrow("D1MinusLinear::setW(newVal,ds) - inconsistent dimensions with problem size for given input matrix W.");
}

void D1MinusLinear::setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds)
{
  unsigned int line = newPtr->size(0);
  unsigned int col  = newPtr->size(1);
  if (line != col)
    RuntimeException::selfThrow("D1MinusLinear::setWPtr(newVal,ds) - newVal is not square.");

  if (!ds)
    RuntimeException::selfThrow("D1MinusLinear::setWPtr(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim();
  if (line != sizeW)
    RuntimeException::selfThrow("D1MinusLinear::setW(newVal,ds) - inconsistent dimension between newVal and dynamical system to be integrated.");

  WMap[ds] = newPtr;
}

void D1MinusLinear::initialize()
{
  OneStepIntegrator::initialize();

  double t0 = simulationLink->model()->t0();

  ConstDSIterator itDS;

  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    initW(t0, *itDS);
    (*itDS)->allocateWorkVector(DynamicalSystem::local_buffer, WMap[*itDS]->size(0));
  }
}

void D1MinusLinear::initW(double t, SP::DynamicalSystem ds)
{
  if (!ds)
    RuntimeException::selfThrow("D1MinusLinear::initW(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("D1MinusLinear::initW(t,ds) - ds does not belong to the OSI.");

  if (WMap.find(ds) != WMap.end())
    RuntimeException::selfThrow("D1MinusLinear::initW(t,ds) - W(ds) is already in the map and has been initialized.");

  Type::Siconos dsType = Type::value(*ds);

  // Lagrangian Systems
  if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    WMap[ds].reset(new SimpleMatrix(*(d->mass())));
    SP::SiconosMatrix W = WMap[ds];
  }
  // Newton Euler Systems
  else if (dsType == Type::NewtonEulerDS) {}
  else RuntimeException::selfThrow("D1MinusLinear::initW(t,ds) - not yet implemented for Dynamical system type: " + dsType);
}

void D1MinusLinear::computeW(double t, SP::DynamicalSystem ds)
{
  assert(ds && "D1MinusLinear::computeW(t,ds) - ds == NULL.");
  assert((WMap.find(ds) != WMap.end()) && "D1MinusLinear::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

  Type::Siconos dsType = Type::value(*ds);
  SP::SiconosMatrix W = WMap[ds];

  // Lagrangian Nonlinear Systems
  if (dsType == Type::LagrangianDS)
  {
    SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    d->computeMass();
    *W = *(d->mass());
  }
  // Lagrangian Linear Systems
  else if (dsType == Type::LagrangianLinearTIDS) {}
  // Newton Euler Systems
  else if (dsType == Type::NewtonEulerDS)
  {
    SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);
    *(d->luW()) = *(d->massMatrix());
    d->luW()->PLUFactorizationInPlace();
  }
  else RuntimeException::selfThrow("D1MinusLinear::computeW(t,ds) - not yet implemented for Dynamical system type: " + dsType);
}

double D1MinusLinear::computeResidu()
{
  double t = simulationLink->nextTime(); // end of the time step
  double told = simulationLink->startingTime(); // beginning of the time step
  double h = simulationLink->timeStep(); // time step length

  DSIterator it; // iteration through the set of DS
  SP::DynamicalSystem ds; // current DS
  Type::Siconos dsType; // type of the current DS

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it;
    dsType = Type::value(*ds);
    SP::SiconosVector residuFree = ds->residuFree(); // pointer constructor

    // Lagrangian Nonlinear Systems
    if (dsType == Type::LagrangianDS)
    {
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      SP::SiconosMatrix Mold = d->mass(); // TODO OLD mass matrix

      // initialize ds->residuFree and predicted right velocity (left limit)
      residuFree->zero();
      SP::SiconosVector vpred(new SimpleVector(*vold));

      // -- LEFT SIDE --
      if (d->fL())
      {
        d->computeFL(told, qold, vold); // left force vector
        scal(-0.5 * h, *(d->fL()), *residuFree, false);

        Mold->PLUForwardBackwardInPlace(*(d->fL())); // predicted right velocity (left limit)
        scal(h, *(d->fL()), *vpred, false);
      }

      // get current information for right side
      SP::SiconosVector q = d->q(); // pointer constructor
      *q = *qold;
      scal(0.5 * h, *vold + *vpred, *q, false);
      d->computeMass();
      SP::SiconosMatrix M = d->mass();

      // -- RIGHT SIDE --
      if (d->fL())
      {
        d->computeFL(t, q, vpred);
        scal(-0.5 * h, *(d->fL()), *residuFree, false);
      }

      *(d->workFree()) = *residuFree; // copy residuFree in workFree
      *(d->workFree()) -= *(d->p(1)); // subtract nonsmooth reaction on impulse level
    }
    // Lagrangian Linear Systems
    else if (dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      SP::SiconosMatrix M = d->mass(); // constant mass
      SP::SiconosMatrix C = d->C(); // constant dissipation
      SP::SiconosMatrix K = d->K(); // constant stiffness

      // current time dependent force
      SP::SiconosVector Fext = d->fExt();

      // initialize ds->residuFree and predicted right velocity (left limit)
      residuFree->zero();
      SP::SiconosVector vpred(new SimpleVector(vold->size()));

      // introduce calculation of contact forces
      SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP

      // -- LEFT SIDE --
      if (Fext)
      {
        d->computeFExt(told);
        scal(-0.5 * h, *Fext, *residuFree, false);
        scal(h, *Fext, *vpred, false);
      }

      if (K)
      {
        SP::SiconosVector dummy(new SimpleVector(qold->size()));
        prod(*K, *qold, *dummy);
        scal(0.5 * h, *dummy, *residuFree, false);
        scal(-h, *dummy, *vpred, false);
      }

      if (C)
      {
        SP::SiconosVector dummy(new SimpleVector(vold->size()));
        prod(*C, *vold, *dummy);
        scal(0.5 * h, *dummy, *residuFree, false);
        scal(-h, *dummy, *vpred, false);
      }

      M->PLUForwardBackwardInPlace(*vpred);
      *vpred += *vold;

      // get current information for right side
      SP::SiconosVector q = d->q(); // pointer constructor
      *q = *qold;
      scal(0.5 * h, *vold + *vpred, *q, false);

      // -- RIGHT SIDE --
      if (Fext)
      {
        d->computeFExt(t);
        scal(-0.5 * h, *Fext, *residuFree, false);
      }

      if (C)
      {
        prod(0.5 * h, *C, *vpred, *residuFree, false);
      }

      if (K)
      {
        prod(0.5 * h, *K, *q, *residuFree, false);
      }

      *(d->workFree()) = *residuFree; // copy residuFree in workFree
      *(d->workFree()) -= *(d->p(1)); // subtract nonsmooth reaction on impulse level
    }
    // Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

      // get left state from memory
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit
      SP::SiconosMatrix Mold = d->massMatrix(); // constant mass
      SP::SiconosMatrix T = d->T(); // explicit usage of T TODO

      // initialize ds->residuFree and predicted right velocity (left limit)
      residuFree->zero();
      SP::SiconosVector vpred(new SimpleVector(*vold));

      // -- LEFT SIDE --
      if (d->fL())
      {
        d->computeFL(told, qold, vold); // left force vector
        scal(-0.5 * h, *(d->fL()), *residuFree, false);

        Mold->PLUForwardBackwardInPlace(*(d->fL())); // predicted right velocity (left limit)
        scal(h, *(d->fL()), *vpred, false);
      }

      // get current information for right side
      SP::SiconosVector q = d->q(); // pointer constructor
      *q = *qold;
      prod(0.5 * h, *T, *vold + *vpred, *q, false);

      // -- RIGHT SIDE --
      if (d->fL())
      {
        d->computeFL(t, q, vpred);
        scal(-0.5 * h, *(d->fL()), *residuFree, false);
      }

      d->updateT();

      *(d->workFree()) = *residuFree; // copy residuFree in workFree
      *(d->workFree()) -= *(d->p(1)); // subtract nonsmooth reaction on impulse level
    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeResidu() - not yet implemented for Dynamical system type: " + dsType);
  }
  return 0.; // there is no Newton iteration and the residuum is assumed to vanish
}

void D1MinusLinear::computeFreeState()
{
  DSIterator it; // iterator through the set of DS
  SP::DynamicalSystem ds; // current DS
  SP::SiconosMatrix W; // W iteration matrix of the current DS
  Type::Siconos dsType ; // type of the current DS

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it;
    dsType = Type::value(*ds);
    W = WMap[ds];

    // Lagrangian Systems
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      // get left state from memory
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // right limit

      // get right information
      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector vfree = d->workFree(); // pointer constructor
      (*vfree) = *(d->residuFree());

      M->PLUForwardBackwardInPlace(*vfree);
      *vfree *= -1.;
      *vfree += *vold;
    }
    // Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

      // get left state from memory
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      SP::SiconosVector vfree = d->workFree(); // pointer constructor
      (*vfree) = *(d->residuFree());

      d->luW()->PLUForwardBackwardInPlace(*vfree);
      *vfree *= -1.;
      *vfree += *vold;
    }
    else
      RuntimeException::selfThrow("D1MinusLinear::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }
}

void D1MinusLinear::computeFreeOutput(SP::UnitaryRelation UR, OneStepNSProblem* osnsp)
{
  SP::OneStepNSProblems allOSNS  = simulationLink->oneStepNSProblems(); // all OSNSP

  // get relation and non smooth law information
  RELATION::TYPES relationType = UR->getRelationType(); // relation
  RELATION::SUBTYPES relationSubType = UR->getRelationSubType();
  unsigned int relativePosition = UR->getRelativePosition();
  unsigned int sizeY = UR->getNonSmoothLawSize(); // related NSL
  SP::DynamicalSystem ds = *(UR->interaction()->dynamicalSystemsBegin()); // related DS

  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  SP::SiconosMatrix C;
  SP::SiconosVector Xq = UR->xq();
  SP::SiconosVector Yp = UR->yp();
  SP::SiconosVector Xfree = UR->workFree();

  assert(Xfree);

  SP::Interaction mainInteraction = UR->interaction();
  assert(mainInteraction);
  assert(mainInteraction->relation());

  // Newton Euler Systems
  if (relationType == NewtonEuler)
  {
    SP::SiconosMatrix CT =  boost::static_pointer_cast<NewtonEulerR>(mainInteraction->relation())->jachqT();

    if (CT)
    {
      assert(Xfree);
      assert(Yp);

      coord[3] = CT->size(1);
      coord[5] = CT->size(1);
      subprod(*CT, *Xfree, *Yp, coord, true);
    }
  }
  // Lagrangian Systems
  else
  {
    C = mainInteraction->relation()->C();

    if (C)
    {
      assert(Xfree);
      assert(Yp);
      assert(Xq);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, *Yp, coord, true);
    }
    if (relationType == Lagrangian)
    {
      C = (boost::static_pointer_cast<LagrangianR>(mainInteraction->relation()))->jachq();
      SP::SiconosMatrix C2 = mainInteraction->relation()->C();

      if (C)
      {
        assert(Xfree);
        assert(Yp);
        assert(Xq);

        coord[3] = C->size(1);
        coord[5] = C->size(1);
        subprod(*C, *Xfree, *Yp, coord, true);
      }

      SP::SiconosMatrix ID(new SimpleMatrix(sizeY, sizeY));
      ID->eye();

      Index xcoord(8);
      xcoord[0] = 0;
      xcoord[1] = sizeY;
      xcoord[2] = 0;
      xcoord[3] = sizeY;
      xcoord[4] = 0;
      xcoord[5] = sizeY;
      xcoord[6] = 0;
      xcoord[7] = sizeY;

      if (relationSubType == RheonomousR)
      {
        if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() == osnsp)
        {
          boost::static_pointer_cast<LagrangianRheonomousR>(UR->interaction()->relation())->computehDot(simulation()->getTkp1());
          subprod(*ID, *(boost::static_pointer_cast<LagrangianRheonomousR>(UR->interaction()->relation())->hDot()), *Yp, xcoord, false);
        }
        else
          RuntimeException::selfThrow("D1MinusLinear::computeFreeOutput(ur,osnsp) - not yet implemented for SICONOS_OSNSP.");
      }
      if (relationSubType == ScleronomousR) {}
    }
  }

  // Lagrangian Systems AND Newton Euler Systems
  if (UR->getRelationType() == Lagrangian || UR->getRelationType() == NewtonEuler)
  {
    SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, UR));
    UR->interaction()->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
  }
}

void D1MinusLinear::updateState(unsigned int level)
{
  DSIterator it; // iterator through the set of DS
  SP::SiconosMatrix W; // W iteration matrix of the current DS

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds];
    Type::Siconos dsType = Type::value(*ds);

    // Lagrangian Systems
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      SP::SiconosVector v = d->velocity(); // pointer constructor
      *v = *(d->p(level)); // value = nonsmooth impulse
      W->PLUForwardBackwardInPlace(*v); // solution for its velocity equivalent
      *v += *(ds->workFree()); // add free velocity
    }
    // Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

      SP::SiconosVector v = d->velocity(); // pointer constructor
      *v = *(d->p(level)); // value = nonsmooth impulse
      d->luW()->PLUForwardBackwardInPlace(*v); // solution for its velocity equivalent
      *v += *(ds->workFree()); // add free velocity
    }
    else RuntimeException::selfThrow("D1MinusLinear::updateState(level) - not yet implemented for Dynamical system type: " + dsType);
  }
}

void D1MinusLinear::display()
{
  OneStepIntegrator::display();

  cout << "====== D1MinusLinear OSI display ======" << endl;
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    cout << "--------------------------------" << endl;
    cout << "--> W of dynamical system number " << (*it)->number() << ": " << endl;
    if (WMap[*it]) WMap[*it]->display();
    else cout << "-> NULL" << endl;
  }
  cout << "================================" << endl;
}

void D1MinusLinear::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
}

D1MinusLinear* D1MinusLinear::convert(OneStepIntegrator* osi)
{
  return dynamic_cast<D1MinusLinear*>(osi);
}

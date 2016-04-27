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
#include "SchatzmanPaoliOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianLinearTIR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "OneStepNSProblem.hpp"

using namespace RELATION;

// --- constructor from a set of data ---
SchatzmanPaoliOSI::SchatzmanPaoliOSI(double theta):
  OneStepIntegrator(OSI::SCHATZMANPAOLIOSI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _theta = theta;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

// --- constructor from a set of data ---
SchatzmanPaoliOSI::SchatzmanPaoliOSI(double theta, double gamma):
  OneStepIntegrator(OSI::SCHATZMANPAOLIOSI), _useGammaForRelation(false)
{
  _theta = theta;
  _gamma = gamma;
  _useGamma = true;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

const SimpleMatrix SchatzmanPaoliOSI::getW(SP::DynamicalSystem ds)
{
  assert(ds &&
         "SchatzmanPaoliOSI::getW(ds): ds == NULL.");
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "SchatzmanPaoliOSI::getW(ds): W[ds] == NULL.");
  return *_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W; // Copy !!
}

SP::SimpleMatrix SchatzmanPaoliOSI::W(SP::DynamicalSystem ds)
{
  assert(ds && "SchatzmanPaoliOSI::W(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W;
;
}

const SimpleMatrix SchatzmanPaoliOSI::getWBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds &&
         "SchatzmanPaoliOSI::getWBoundaryConditions(ds): ds == NULL.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).W &&
         "SchatzmanPaoliOSI::getWBoundaryConditions(ds): WBoundaryConditions[ds] == NULL.");
  return *(_dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions); // Copy !!
}

SP::SiconosMatrix SchatzmanPaoliOSI::WBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds && "SchatzmanPaoliOSI::WBoundaryConditions(ds): ds == NULL.");
  return _dynamicalSystemsGraph->properties(_dynamicalSystemsGraph->descriptor(ds)).WBoundaryConditions;
}


void SchatzmanPaoliOSI::initialize(Model& m)
{
  OneStepIntegrator::initialize(m);
  // Get initial time
  double t0 = _simulation->startingTime();
  // Compute W(t0) for all ds
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    Type::Siconos dsType = Type::value(*ds);
    if (dsType == Type::LagrangianLinearTIDS)
    {
      // Computation of the first step for starting
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);

      SP::SiconosVector q0  = d->q0();
      SP::SiconosVector q  = d->q();
      SP::SiconosVector v0  = d->velocity0();
      SP::SiconosVector velocity  = d->velocity();

      //  std::cout << " q0 = " << std::endl;
      // q0->display();
      //  std::cout << " v0 = " << std::endl;
      // v0->display();
      // We first swap the initial value contained in q and v after initialization.

      d->qMemory()->swap(*q);
      d->velocityMemory()->swap(*velocity);

      // we compute the new state values
      double h = _simulation->timeStep();
      *q = *q0 + h* * v0;
      //*velocity=*velocity; we do nothing for the velocity

      // This value will swapped when OneStepIntegrator::saveInMemory will be called
      // by the rest of  Simulation::initialize (_eventsManager->preUpdate();)

      // SP::SiconosVector qprev = d->qMemory()->getSiconosVector(0);
      // SP::SiconosVector qprev2 = d->qMemory()->getSiconosVector(1);
      // SP::SiconosVector vprev = d->velocityMemory()->getSiconosVector(0);
      //  std::cout << " qprev = " << std::endl;
      // qprev->display();
      //  std::cout << " qprev2 = " << std::endl;
      // qprev2->display();
      //  std::cout << " vprev = " << std::endl;
      // vprev->display();



    }
    // Memory allocation for workX. workX[ds*] corresponds to xfree (or vfree in lagrangian case).
    // workX[*itDS].reset(new SiconosVector((*itDS)->getDim()));

    // W initialization
    initW(t0, ds, *dsi);

    //      if ((*itDS)->getType() == Type::LagrangianDS || (*itDS)->getType() == Type::FirstOrderNonLinearDS)
    ds->allocateWorkVector(DynamicalSystem::local_buffer,_dynamicalSystemsGraph->properties(*dsi).W->size(0));
  }
}
void SchatzmanPaoliOSI::initW(double t, SP::DynamicalSystem ds, DynamicalSystemsGraph::VDescriptor& dsv)
{
  // This function:
  // - allocate memory for a matrix W

  if (!ds)
    RuntimeException::selfThrow("SchatzmanPaoliOSI::initW(t,ds) - ds == NULL");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("SchatzmanPaoliOSI::initW(t,ds) - ds does not belong to the OSI.");

  if (_dynamicalSystemsGraph->properties(dsv).W)
    RuntimeException::selfThrow("SchatzmanPaoliOSI::initW(t,ds) - W(ds) is already in the map and has been initialized.");


  //unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  // Memory allocation for W

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);


  // 1 - Lagrangian non linear systems
  if (dsType == Type::LagrangianDS)
  {

    RuntimeException::selfThrow("SchatzmanPaoliOSI::initW - not yet implemented for Dynamical system type :" + dsType);

  }
  // 4 - Lagrangian linear systems
  else if (dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);
    SP::SiconosMatrix K = d->K();
    SP::SiconosMatrix C = d->C();
    _dynamicalSystemsGraph->properties(dsv).W.reset(new SimpleMatrix(*d->mass())); //*W = *d->mass();
    SP::SiconosMatrix W = _dynamicalSystemsGraph->properties(dsv).W;

    if (C)
      scal(1 / 2.0 * h * _theta, *C, *W, false); // W += 1/2.0*h*_theta *C

    if (K)
      scal(h * h * _theta * _theta, *K, *W, false); // W = h*h*_theta*_theta*K

    // WBoundaryConditions initialization
    if (d->boundaryConditions())
      initWBoundaryConditions(d,dsv);


  }

  // === ===
  else if (dsType == Type::NewtonEulerDS)
  {
    RuntimeException::selfThrow("SchatzmanPaoliOSI::initW - not yet implemented for Dynamical system type :" + dsType);
  }
  else RuntimeException::selfThrow("SchatzmanPaoliOSI::initW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.

}


void SchatzmanPaoliOSI::initWBoundaryConditions(SP::DynamicalSystem ds, DynamicalSystemsGraph::VDescriptor& dsv)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("SchatzmanPaoliOSI::initWBoundaryConditions(t,ds) - ds == NULL");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("SchatzmanPaoliOSI::initWBoundaryConditions(t,ds) - ds does not belong to the OSI.");

  Type::Siconos dsType = Type::value(*ds);

  RuntimeException::selfThrow("SchatzmanPaoliOSI::initWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void SchatzmanPaoliOSI::computeWBoundaryConditions(SP::DynamicalSystem ds, SiconosMatrix& WBoundaryConditions)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initWBoundaryConditions.

  assert(ds &&
         "SchatzmanPaoliOSI::computeWBoundaryConditions(t,ds) - ds == NULL");

  Type::Siconos dsType = Type::value(*ds);

  RuntimeException::selfThrow("SchatzmanPaoliOSI::computeWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void SchatzmanPaoliOSI::computeW(double t, SP::DynamicalSystem ds, SiconosMatrix& W)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  assert(ds &&
         "SchatzmanPaoliOSI::computeW(t,ds) - ds == NULL");

  //double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  // 1 - Lagrangian non linear systems
  if (dsType == Type::LagrangianDS)
  {

    RuntimeException::selfThrow("SchatzmanPaoliOSI::computeW - not yet implemented for Dynamical system type :" + dsType);

  }
  // 4 - Lagrangian linear systems
  else if (dsType == Type::LagrangianLinearTIDS)
  {
    // Nothing: W does not depend on time.
  }

  // === ===
  else if (dsType == Type::NewtonEulerDS)
  {
    RuntimeException::selfThrow("SchatzmanPaoliOSI::computeW - not yet implemented for Dynamical system type :" + dsType);
  }
  else RuntimeException::selfThrow("SchatzmanPaoliOSI::computeW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}




double SchatzmanPaoliOSI::computeResidu()
{

  // This function is used to compute the residu for each "SchatzmanPaoliOSI-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds); // Its type
    SP::SiconosVector residuFree = ds->workspace(DynamicalSystem::freeresidu);

    // 1 - Lagrangian Non Linear Systems
    if (dsType == Type::LagrangianDS)
    {
      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
    // 2 - Lagrangian Linear Systems
    else if (dsType == Type::LagrangianLinearTIDS)
    {
      // ResiduFree =  M(-q_{k}+q_{k-1})  + h^2 (K q_k)+  h^2 C (\theta \Frac{q_k-q_{k-1}}{2h}+ (1-\theta) v_k))  (1)
      // This formulae is only valid for the first computation of the residual for q = q_k
      // otherwise the complete formulae must be applied, that is
      // ResiduFree   M(q-2q_{k}+q_{k-1})  + h^2 (K(\theta q+ (1-\theta) q_k)))+  h^2 C (\theta \Frac{q-q_{k-1}}{2h}+ (1-\theta) v_k))  (2)
      // for q != q_k, the formulae (1) is wrong.
      // in the sequel, only the equation (1) is implemented

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector q_k = d->qMemory()->getSiconosVector(0); // q_k
      SP::SiconosVector q_k_1 = d->qMemory()->getSiconosVector(1); // q_{k-1}
      SP::SiconosVector v_k = d->velocityMemory()->getSiconosVector(0); //v_k
      //  std::cout << "SchatzmanPaoliOSI::computeResidu - q_k_1 =" <<std::endl;
      // q_k_1->display();
      //  std::cout << "SchatzmanPaoliOSI::computeResidu - q_k =" <<std::endl;
      // q_k->display();
      //  std::cout << "SchatzmanPaoliOSI::computeResidu - v_k =" <<std::endl;
      // v_k->display();

      // --- ResiduFree computation Equation (1) ---
      residuFree->zero();
      double coeff;
      // -- No need to update W --

      //SP::SiconosVector v = d->velocity(); // v = v_k,i+1

      SP::SiconosMatrix M = d->mass();
      prod(*M, (*q_k_1 - *q_k), *residuFree); // residuFree = M(-q_{k}+q_{k-1})

      SP::SiconosMatrix K = d->K();
      if (K)
      {
        prod(h * h, *K, *q_k, *residuFree, false); // residuFree += h^2*K*qi
      }

      SP::SiconosMatrix C = d->C();
      if (C)
        prod(h * h, *C, (1.0 / (2.0 * h)*_theta * (*q_k - *q_k_1) + (1.0 - _theta)* *v_k)  , *residuFree, false);
      // residufree += h^2 C (\theta \Frac{q-q_{k-1}}{2h}+ (1-\theta) v_k))


      SP::SiconosVector Fext = d->fExt();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = -h * h * (1 - _theta);
        scal(coeff, *Fext, *residuFree, false); // residufree -= h^2*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = -h * h * _theta;
        scal(coeff, *Fext, *residuFree, false); // residufree -= h^2*_theta * fext(ti+1)
      }



      //  std::cout << "SchatzmanPaoliOSI::ComputeResidu LagrangianLinearTIDS residufree :"  << std::endl;
      // residuFree->display();


      (* d->workspace(DynamicalSystem::free)) = *residuFree; // copy residuFree in Workfree
      if (d->p(0))
        *(d->workspace(DynamicalSystem::free)) -= *d->p(0); // Compute Residu in Workfree Notation !!

      //  std::cout << "SchatzmanPaoliOSI::ComputeResidu LagrangianLinearTIDS p(0) :"  << std::endl;
      //  if (d->p(0))
      //    d->p(0)->display();
      //  else
      //     std::cout << " p(0) :"  << std::endl;
      //  std::cout << "SchatzmanPaoliOSI::ComputeResidu LagrangianLinearTIDS residu :"  << std::endl;
      // d->workspace(DynamicalSystem::free)->display();



      //     normResidu = d->workspace(DynamicalSystem::free)->norm2();
      normResidu = 0.0; // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();

    }
    else if (dsType == Type::NewtonEulerDS)
    {
      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
    else
      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    if (normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void SchatzmanPaoliOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  //double t = _simulation->nextTime(); // End of the time step
  //double told = _simulation->startingTime(); // Beginning of the time step
  //double h = t-told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  SP::DynamicalSystem ds; // Current Dynamical System.
  SP::SiconosMatrix W; // W SchatzmanPaoliOSI matrix of the current DS.
  Type::Siconos dsType ; // Type of the current DS.

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;

    ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds); // Its type
    W =  _dynamicalSystemsGraph->properties(*dsi).W; // Its W SchatzmanPaoliOSI matrix of iteration.

    //1 - Lagrangian Non Linear Systems
    if (dsType == Type::LagrangianDS)
    {

      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
    }
    // 2 - Lagrangian Linear Systems
    else if (dsType == Type::LagrangianLinearTIDS)
    {
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // vFree = v_i + W^{-1} ResiduFree    // with
      // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 + h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0); // q_k
      //   SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); //v_k

      // --- ResiduFree computation ---

      // vFree pointer is used to compute and save ResiduFree in this first step.


      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SP::SiconosVector qfree = d->workspace(DynamicalSystem::free);//workX[d];
      (*qfree) = *(d->workspace(DynamicalSystem::freeresidu));

      W->PLUForwardBackwardInPlace(*qfree);
      *qfree *= -1.0;
      *qfree += *qold;

    }
    // 3 - Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
    }
    else
      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }

}

void SchatzmanPaoliOSI::prepareNewtonIteration(double time)
{
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    computeW(time, ds, *_dynamicalSystemsGraph->properties(*dsi).W );
  }
}

struct SchatzmanPaoliOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem* _osnsp;
  SP::Interaction _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter) :
    _osnsp(p), _inter(inter) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = _inter->nonSmoothLaw()->size();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    // Only the normal part is multiplied by e
    SP::SiconosVector y_k_1 ;
    y_k_1 = _inter->yMemory(_osnsp->inputOutputLevel())->getSiconosVector(1);


    //  std::cout << "y_k_1 " << std::endl;
    // y_k_1->display();
    subscal(e, *y_k_1, *(_inter->yForNSsolver()), subCoord, false);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    SP::SiconosVector y_k_1 ;
    y_k_1 = _inter->yMemory(_osnsp->inputOutputLevel())->getSiconosVector(1);
    (*_inter->yForNSsolver())(0) +=  e * (*y_k_1)(0);

  }
  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
};


void SchatzmanPaoliOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */

  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);
  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
  // Get relation and non smooth law types
  RELATION::TYPES relationType = inter->relation()->getType();
  RELATION::SUBTYPES relationSubType = inter->relation()->getSubType();
  unsigned int sizeY = inter->nonSmoothLaw()->size();

  unsigned int relativePosition = 0;



  Index coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = 0;
  coord[7] = sizeY;
  SP::SiconosMatrix  C;
  SP::SiconosMatrix  D;
  SP::SiconosMatrix  F;
  SP::BlockVector deltax;
  SiconosVector& yForNSsolver = *inter->yForNSsolver();
  SP::SiconosVector e;
  SP::BlockVector Xfree;

  if (relationType == NewtonEuler)
  {
    Xfree = DSlink[NewtonEulerR::xfree];
  }
  else if (relationType == Lagrangian)
  {
    Xfree = DSlink[LagrangianR::xfree];
  }

  assert(Xfree);

  assert(Xfree);


  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if (relationSubType == LinearTIR)
  {

    if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() != osnsp)
      RuntimeException::selfThrow("SchatzmanPaoliOSI::computeFreeOutput not yet implemented for SICONOS_OSNSP ");

    C = mainInteraction->relation()->C();

    if (C)
    {

      assert(Xfree);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interactionfor each ds of the
      // current Interaction.

      if (_useGammaForRelation)
      {
        assert(deltax);
        subprod(*C, *deltax, yForNSsolver, coord, true);
      }
      else
      {
        subprod(*C, *Xfree, yForNSsolver, coord, true);
        //        subprod(*C,*(*(mainInteraction->dynamicalSystemsBegin()))->workspace(DynamicalSystem::free),*Yp,coord,true);
        //        if (mainInteraction->dynamicalSystems()->size() == 2)
        //        {
        //          subprod(*C,*(*++(mainInteraction->dynamicalSystemsBegin()))->workspace(DynamicalSystem::free),*Yp,coord,false);
        //        }
      }

    }
    SP::LagrangianLinearTIR ltir = std11::static_pointer_cast<LagrangianLinearTIR> (mainInteraction->relation());
    e = ltir->e();
    if (e)
    {
      yForNSsolver += *e;
    }

  }
  else
    RuntimeException::selfThrow("SchatzmanPaoliOSI::ComputeFreeOutput not yet implemented  for relation of Type : " + relationType);



  if (inter->relation()->getSubType() == LinearTIR)
  {
    SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, inter));
    inter->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
  }


}
void SchatzmanPaoliOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  RuntimeException::selfThrow("SchatzmanPaoliOSI::integrate - not yet implemented :");
}

void SchatzmanPaoliOSI::updateState(const unsigned int level)
{

  double h = _simulation->timeStep();

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  SP::SiconosMatrix W;
  SP::DynamicalSystem ds;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);

    W = _dynamicalSystemsGraph->properties(*dsi).W;
    // Get the DS type

    Type::Siconos dsType = Type::value(*ds);

    // 1 - Lagrangian Systems
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      // get dynamical system
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);

      //    SiconosVector *vfree = d->velocityFree();
      SP::SiconosVector q = d->q();
      bool baux = dsType == Type::LagrangianDS && useRCC && _simulation->relativeConvergenceCriterionHeld();
      if (level != LEVELMAX)
      {
        // To compute q, we solve W(q - qfree) = p
        if (d->p(level))
        {
          *q = *d->p(level); // q = p
          W->PLUForwardBackwardInPlace(*q);
        }

        // if (d->boundaryConditions())
        //   for (vector<unsigned int>::iterator
        //        itindex = d->boundaryConditions()->velocityIndices()->begin() ;
        //        itindex != d->boundaryConditions()->velocityIndices()->end();
        //        ++itindex)
        //     v->setValue(*itindex, 0.0);
        *q +=  * ds->workspace(DynamicalSystem::free);

      }
      else
        *q =  * ds->workspace(DynamicalSystem::free);



      // Computation of the velocity

      SP::SiconosVector v = d->velocity();
      SP::SiconosVector q_k_1 = d->qMemory()->getSiconosVector(1); // q_{k-1}

      //  std::cout << "SchatzmanPaoliOSI::updateState - q_k_1 =" <<std::endl;
      // q_k_1->display();
      //  std::cout << "SchatzmanPaoliOSI::updateState - q =" <<std::endl;
      // q->display();

      *v = 1.0 / (2.0 * h) * (*q - *q_k_1);
      //  std::cout << "SchatzmanPaoliOSI::updateState - v =" <<std::endl;
      // v->display();

      // int bc=0;
      // SP::SiconosVector columntmp(new SiconosVector(ds->getDim()));

      // if (d->boundaryConditions())
      // {
      //   for (vector<unsigned int>::iterator  itindex = d->boundaryConditions()->velocityIndices()->begin() ;
      //        itindex != d->boundaryConditions()->velocityIndices()->end();
      //        ++itindex)
      //   {
      //     _WBoundaryConditionsMap[ds]->getCol(bc,*columntmp);
      //     /*\warning we assume that W is symmetric in the Lagrangian case*/
      //     double value = - inner_prod(*columntmp, *v);
      //     value += (d->p(level))->getValue(*itindex);
      //     /* \warning the computation of reactionToBoundaryConditions take into
      //        account the contact impulse but not the external and internal forces.
      //        A complete computation of the residue should be better */
      //     d->reactionToBoundaryConditions()->setValue(bc,value) ;
      //     bc++;
      //   }

      if (baux)
      {
        ds->subWorkVector(q, DynamicalSystem::local_buffer);
        double aux = ((ds->workspace(DynamicalSystem::local_buffer))->norm2()) / (ds->normRef());
        if (aux > RelativeTol)
          _simulation->setRelativeConvergenceCriterionHeld(false);
      }

    }
    //2 - Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      //  // get dynamical system
      //       SP::NewtonEulerDS d = std11::static_pointer_cast<NewtonEulerDS> (ds);
      //       SP::SiconosVector v = d->velocity();
      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"SchatzmanPaoliOSI::updatestate prev v"<<endl;
      //       v->display();
      // #endif

      //       /*d->p has been fill by the Relation->computeInput, it contains
      //            B \lambda _{k+1}*/
      //       *v = *d->p(level); // v = p
      //       d->luW()->PLUForwardBackwardInPlace(*v);

      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"SchatzmanPaoliOSI::updatestate hWB lambda"<<endl;
      //       v->display();
      // #endif

      //       *v +=  * ds->workspace(DynamicalSystem::free);

      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"SchatzmanPaoliOSI::updatestate work free"<<endl;
      //       ds->workspace(DynamicalSystem::free)->display();
      //       std::cout<<"SchatzmanPaoliOSI::updatestate new v"<<endl;
      //       v->display();
      // #endif
      //       //compute q
      //       //first step consists in computing  \dot q.
      //       //second step consists in updating q.
      //       //
      //       SP::SiconosMatrix T = d->T();
      //       SP::SiconosVector dotq = d->dotq();
      //       prod(*T,*v,*dotq,true);
      //       // std::cout<<"SchatzmanPaoliOSI::updateState v"<<endl;
      //       // v->display();
      //       // std::cout<<"SchatzmanPaoliOSI::updateState dotq"<<endl;
      //       // dotq->display();




      //       SP::SiconosVector q = d->q();

      //       //  -> get previous time step state
      //       SP::SiconosVector dotqold = d->dotqMemory()->getSiconosVector(0);
      //       SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      //       // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
      //       double coeff = h*_theta;
      //       scal(coeff, *dotq, *q) ; // q = h*theta*v
      //       coeff = h*(1-_theta);
      //       scal(coeff,*dotqold,*q,false); // q += h(1-theta)*vold
      //       *q += *qold;
      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       std::cout<<"new q before normalizing"<<endl;
      //       q->display();
      // #endif

      //       //q[3:6] must be normalized
      //       d->normalizeq();
      //       dotq->setValue(3,(q->getValue(3)-qold->getValue(3))/h);
      //       dotq->setValue(4,(q->getValue(4)-qold->getValue(4))/h);
      //       dotq->setValue(5,(q->getValue(5)-qold->getValue(5))/h);
      //       dotq->setValue(6,(q->getValue(6)-qold->getValue(6))/h);
      //       d->updateT();
      RuntimeException::selfThrow("SchatzmanPaoliOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
    }
    else RuntimeException::selfThrow("SchatzmanPaoliOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}


void SchatzmanPaoliOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== SchatzmanPaoliOSI OSI display ======" <<std::endl;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--------------------------------" <<std::endl;
    std::cout << "--> W of dynamical system number " << ds->number() << ": " <<std::endl;
    if (_dynamicalSystemsGraph->properties(*dsi).W)  _dynamicalSystemsGraph->properties(*dsi).W->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "--> and corresponding theta is: " << _theta <<std::endl;
  }
  std::cout << "================================" <<std::endl;
}

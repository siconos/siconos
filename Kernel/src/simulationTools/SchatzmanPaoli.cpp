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
#include "SchatzmanPaoli.hpp"
//#include "SchatzmanPaoliXML.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "NewtonEulerR.hpp"
#include "LagrangianRheonomousR.hpp"
#include "LagrangianLinearTIR.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderLinearR.hpp"
#include "NewtonImpactNSL.hpp"
#include "MultipleImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"

using namespace std;
using namespace RELATION;

// --- constructor from a minimum set of data ---
SchatzmanPaoli::SchatzmanPaoli(SP::DynamicalSystem newDS, double newTheta) :
  OneStepIntegrator(OSI::SCHATZMANPAOLI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  OSIDynamicalSystems->insert(newDS);
  _theta = newTheta;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

// --- constructor from a set of data ---
SchatzmanPaoli::SchatzmanPaoli(DynamicalSystemsSet& newDS, double newTheta):
  OneStepIntegrator(OSI::SCHATZMANPAOLI, newDS), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _theta = newTheta;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

// --- constructor from a set of data ---
SchatzmanPaoli::SchatzmanPaoli(double newTheta):
  OneStepIntegrator(OSI::SCHATZMANPAOLI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _theta = newTheta;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

// --- constructor from a minimum set of data ---
SchatzmanPaoli::SchatzmanPaoli(SP::DynamicalSystem newDS, double newTheta, double newGamma) :
  OneStepIntegrator(OSI::SCHATZMANPAOLI), _useGammaForRelation(false)
{
  OSIDynamicalSystems->insert(newDS);
  _theta = newTheta;
  _gamma = newGamma;
  _useGamma = true;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

// --- constructor from a set of data ---
SchatzmanPaoli::SchatzmanPaoli(DynamicalSystemsSet& newDS, double newTheta, double newGamma):
  OneStepIntegrator(OSI::SCHATZMANPAOLI, newDS), _useGammaForRelation(false)
{
  _theta = newTheta;
  _gamma = newGamma;
  _useGamma = true;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}

// --- constructor from a set of data ---
SchatzmanPaoli::SchatzmanPaoli(double newTheta, double newGamma):
  OneStepIntegrator(OSI::SCHATZMANPAOLI), _useGammaForRelation(false)
{
  _theta = newTheta;
  _gamma = newGamma;
  _useGamma = true;
  _sizeMem = SCHATZMANPAOLISTEPSINMEMORY ;
}


// Note: OSIDynamicalSystems and thetaMap must disappear
void SchatzmanPaoli::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
}
const SimpleMatrix SchatzmanPaoli::getW(SP::DynamicalSystem ds)
{
  assert(ds &&
         "SchatzmanPaoli::getW(ds): ds == NULL.");
  //    return *(WMap[0]);
  assert(WMap[ds] &&
         "SchatzmanPaoli::getW(ds): W[ds] == NULL.");
  return *(WMap[ds]); // Copy !!
}

SP::SimpleMatrix SchatzmanPaoli::W(SP::DynamicalSystem ds)
{
  assert(ds && "SchatzmanPaoli::W(ds): ds == NULL.");
  //  return WMap[0];
  //  if(WMap[ds]==NULL)
  //    RuntimeException::selfThrow("SchatzmanPaoli::W(ds): W[ds] == NULL.");
  return WMap[ds];
}

void SchatzmanPaoli::setW(const SiconosMatrix& newValue, SP::DynamicalSystem ds)
{
  // Check if ds is in the OSI
  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("SchatzmanPaoli::setW(newVal,ds) - ds does not belong to this Integrator ...");

  // Check dimensions consistency
  unsigned int line = newValue.size(0);
  unsigned int col  = newValue.size(1);

  if (line != col) // Check that newValue is square
    RuntimeException::selfThrow("SchatzmanPaoli::setW(newVal,ds) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("SchatzmanPaoli::setW(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("SchatzmanPaoli::setW(newVal,ds) - unconsistent dimension between newVal and dynamical system to be integrated ");

  // Memory allocation for W, if required
  if (!WMap[ds]) // allocate a new W if required
  {
    WMap[ds].reset(new SimpleMatrix(newValue));
  }
  else  // or fill-in an existing one if dimensions are consistent.
  {
    if (line == WMap[ds]->size(0) && col == WMap[ds]->size(1))
      *(WMap[ds]) = newValue;
    else
      RuntimeException::selfThrow("SchatzmanPaoli - setW: inconsistent dimensions with problem size for given input matrix W");
  }
}

void SchatzmanPaoli::setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds)
{
  unsigned int line = newPtr->size(0);
  unsigned int col  = newPtr->size(1);
  if (line != col) // Check that newPtr is square
    RuntimeException::selfThrow("SchatzmanPaoli::setWPtr(newVal) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("SchatzmanPaoli::setWPtr(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("SchatzmanPaoli::setW(newVal) - unconsistent dimension between newVal and dynamical system to be integrated ");

  WMap[ds] = newPtr;                  // link with new pointer
}



const SimpleMatrix SchatzmanPaoli::getWBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds &&
         "SchatzmanPaoli::getWBoundaryConditions(ds): ds == NULL.");
  //    return *(WBoundaryConditionsMap[0]);
  assert(_WBoundaryConditionsMap[ds] &&
         "SchatzmanPaoli::getWBoundaryConditions(ds): WBoundaryConditions[ds] == NULL.");
  return *(_WBoundaryConditionsMap[ds]); // Copy !!
}

SP::SiconosMatrix SchatzmanPaoli::WBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds && "SchatzmanPaoli::WBoundaryConditions(ds): ds == NULL.");
  //  return WBoundaryConditionsMap[0];
  //  if(WBoundaryConditionsMap[ds]==NULL)
  //    RuntimeException::selfThrow("SchatzmanPaoli::WBoundaryConditions(ds): W[ds] == NULL.");
  return _WBoundaryConditionsMap[ds];
}


void SchatzmanPaoli::initialize()
{
  OneStepIntegrator::initialize();
  // Get initial time
  double t0 = simulationLink->model()->t0();
  // Compute W(t0) for all ds
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    Type::Siconos dsType = Type::value(*(*itDS));
    if (dsType == Type::LagrangianLinearTIDS)
    {
      // Computation of the first step for starting
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (*itDS);

      SP::SiconosVector q0  = d->q0();
      SP::SiconosVector q  = d->q();
      SP::SiconosVector v0  = d->velocity0();
      SP::SiconosVector velocity  = d->velocity();

      // std::cout << " q0 = " << std::endl;
      // q0->display();
      // std::cout << " v0 = " << std::endl;
      // v0->display();
      // We first swap the initial value contained in q and v after initialization.

      d->qMemory()->swap(q);
      d->velocityMemory()->swap(velocity);

      // we compute the new state values
      double h = simulationLink->timeStep();
      *q = *q0 + h* * v0;
      //*velocity=*velocity; we do nothing for the velocity

      // This value will swapped when OneStepIntegrator::saveInMemory will be called
      // by the rest of  Simulation::initialize (_eventsManager->preUpdate();)

      // SP::SiconosVector qprev = d->qMemory()->getSiconosVector(0);
      // SP::SiconosVector qprev2 = d->qMemory()->getSiconosVector(1);
      // SP::SiconosVector vprev = d->velocityMemory()->getSiconosVector(0);
      // std::cout << " qprev = " << std::endl;
      // qprev->display();
      // std::cout << " qprev2 = " << std::endl;
      // qprev2->display();
      // std::cout << " vprev = " << std::endl;
      // vprev->display();



    }
    // Memory allocation for workX. workX[ds*] corresponds to xfree (or vfree in lagrangian case).
    // workX[*itDS].reset(new SimpleVector((*itDS)->getDim()));

    // W initialization
    initW(t0, *itDS);

    //      if ((*itDS)->getType() == Type::LagrangianDS || (*itDS)->getType() == Type::FirstOrderNonLinearDS)
    (*itDS)->allocateWorkVector(DynamicalSystem::local_buffer, WMap[*itDS]->size(0));
  }
}
void SchatzmanPaoli::initW(double t, SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix W
  // - insert this matrix into WMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("SchatzmanPaoli::initW(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("SchatzmanPaoli::initW(t,ds) - ds does not belong to the OSI.");

  if (WMap.find(ds) != WMap.end())
    RuntimeException::selfThrow("SchatzmanPaoli::initW(t,ds) - W(ds) is already in the map and has been initialized.");


  //unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  // Memory allocation for W
  //  WMap[ds].reset(new SimpleMatrix(sizeW,sizeW));
  //   SP::SiconosMatrix W = WMap[ds];

  double h = simulationLink->timeStep();
  Type::Siconos dsType = Type::value(*ds);


  // 1 - Lagrangian non linear systems
  if (dsType == Type::LagrangianDS)
  {
    // SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    // SP::SiconosMatrix K = d->jacobianqFL(); // jacobian according to q
    // SP::SiconosMatrix C = d->jacobianqDotFL(); // jacobian according to velocity
    // WMap[ds].reset(new SimpleMatrix(*d->mass())); //*W = *d->mass();

    // SP::SiconosMatrix W = WMap[ds];

    // if (C)
    //   scal(-h*_theta, *C,*W,false); // W -= h*_theta*C

    // if (K)
    //   scal(-h*h*_theta*_theta,*K,*W,false); //*W -= h*h*_theta*_theta**K;

    // // WBoundaryConditions initialization
    // if (d->boundaryConditions())
    //   initWBoundaryConditions(d);
    RuntimeException::selfThrow("SchatzmanPaoli::initW - not yet implemented for Dynamical system type :" + dsType);

  }
  // 4 - Lagrangian linear systems
  else if (dsType == Type::LagrangianLinearTIDS)
  {
    SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);
    SP::SiconosMatrix K = d->K();
    SP::SiconosMatrix C = d->C();
    WMap[ds].reset(new SimpleMatrix(*d->mass())); //*W = *d->mass();
    SP::SiconosMatrix W = WMap[ds];

    if (C)
      scal(1 / 2.0 * h * _theta, *C, *W, false); // W += 1/2.0*h*_theta *C

    if (K)
      scal(h * h * _theta * _theta, *K, *W, false); // W = h*h*_theta*_theta*K

    // WBoundaryConditions initialization
    if (d->boundaryConditions())
      initWBoundaryConditions(d);


  }

  // === ===
  else if (dsType == Type::NewtonEulerDS)
  {
    //WMap[ds].reset(new SimpleMatrix(3,3));
    RuntimeException::selfThrow("SchatzmanPaoli::initW - not yet implemented for Dynamical system type :" + dsType);
  }
  else RuntimeException::selfThrow("SchatzmanPaoli::initW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.

}


void SchatzmanPaoli::initWBoundaryConditions(SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("SchatzmanPaoli::initWBoundaryConditions(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("SchatzmanPaoli::initWBoundaryConditions(t,ds) - ds does not belong to the OSI.");

  Type::Siconos dsType = Type::value(*ds);


  // if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
  // {



  //   if (_WBoundaryConditionsMap.find(ds) != _WBoundaryConditionsMap.end())
  //     RuntimeException::selfThrow("SchatzmanPaoli::initWBoundaryConditions(t,ds) - WBoundaryConditions(ds) is already in the map and has been initialized.");

  //   // Memory allocation for WBoundaryConditions
  //   unsigned int sizeWBoundaryConditions = ds->getDim(); // n for first order systems, ndof for lagrangian.

  //   SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

  //   unsigned int numberBoundaryConditions = d->boundaryConditions()->velocityIndices()->size();
  //   _WBoundaryConditionsMap[ds].reset(new SimpleMatrix(sizeWBoundaryConditions,numberBoundaryConditions));
  //   computeWBoundaryConditions(ds);
  // }
  // else
  RuntimeException::selfThrow("SchatzmanPaoli::initWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void SchatzmanPaoli::computeWBoundaryConditions(SP::DynamicalSystem ds)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initWBoundaryConditions.

  assert(ds &&
         "SchatzmanPaoli::computeWBoundaryConditions(t,ds) - ds == NULL");

  Type::Siconos dsType = Type::value(*ds);

  // if (dsType == Type::LagrangianLinearTIDS || dsType == Type::LagrangianDS)
  // {

  //   assert((_WBoundaryConditionsMap.find(ds) != _WBoundaryConditionsMap.end()) &&
  //          "SchatzmanPaoli::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

  //   SP::SimpleMatrix WBoundaryConditions = _WBoundaryConditionsMap[ds];

  //   SP::SiconosVector columntmp(new SimpleVector(ds->getDim()));

  //   int columnindex =0;

  //   vector<unsigned int>::iterator itindex;

  //   SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

  //   for (itindex = d->boundaryConditions()->velocityIndices()->begin() ;
  //        itindex != d->boundaryConditions()->velocityIndices()->end();
  //        ++itindex)
  //   {

  //     WMap[ds]->getCol(*itindex, *columntmp);
  //     /*\warning we assume that W is symmetric in the Lagrangian case
  //       we store only the column and not the row */

  //     WBoundaryConditions->setCol(columnindex, *columntmp);
  //     double diag = (*columntmp)(*itindex);
  //     columntmp->zero();
  //     (*columntmp)(*itindex)=diag;

  //     WMap[ds]->setCol(*itindex, *columntmp);
  //     WMap[ds]->setRow(*itindex, *columntmp);

  //     columnindex ++;
  //   }
  // }
  // else
  RuntimeException::selfThrow("SchatzmanPaoli::computeWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void SchatzmanPaoli::computeW(double t, SP::DynamicalSystem ds)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, WMap[ds] is supposed to exist and not to be null
  // Memory allocation has been done during initW.

  assert(ds &&
         "SchatzmanPaoli::computeW(t,ds) - ds == NULL");

  assert((WMap.find(ds) != WMap.end()) &&
         "SchatzmanPaoli::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

  //double h = simulationLink->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  SP::SiconosMatrix W = WMap[ds];

  // 1 - Lagrangian non linear systems
  if (dsType == Type::LagrangianDS)
  {
    // SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    // SP::SiconosMatrix K = d->jacobianqFL(); // jacobian according to q
    // SP::SiconosMatrix C = d->jacobianqDotFL(); // jacobian according to velocity

    // d->computeMass();
    // *W = *d->mass();

    // if (C)
    // {
    //   d->computeJacobianqDotFL(t);
    //   scal(-h*_theta, *C,*W,false); // W -= h*_theta*C
    // }

    // if (K)
    // {
    //   d->computeJacobianqFL(t);
    //   scal(-h*h*_theta*_theta,*K,*W,false); //*W -= h*h*_theta*_theta**K;
    // }
    RuntimeException::selfThrow("SchatzmanPaoli::computeW - not yet implemented for Dynamical system type :" + dsType);

  }
  // 4 - Lagrangian linear systems
  else if (dsType == Type::LagrangianLinearTIDS)
  {
    // Nothing: W does not depend on time.
  }

  // === ===
  else if (dsType == Type::NewtonEulerDS)
  {
    // SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);
    // d->computeJacobianvFL(t);
    // double thetaFL=_theta;
    // *(d->luW())=*(d->jacobianvFL());
    // scal(h*thetaFL,*(d->jacobianvFL()),*(d->luW()),true);
    // *(d->luW())+=*(d->massMatrix());

    // //cout<<"SchatzmanPaoli::computeW luW before LUFact\n";
    // //d->luW()->display();

    // d->luW()->PLUFactorizationInPlace();
    RuntimeException::selfThrow("SchatzmanPaoli::computeW - not yet implemented for Dynamical system type :" + dsType);
  }
  else RuntimeException::selfThrow("SchatzmanPaoli::computeW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}




double SchatzmanPaoli::computeResidu()
{

  // This function is used to compute the residu for each "SchatzmanPaoli-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double t = simulationLink->nextTime(); // End of the time step
  double told = simulationLink->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it;
  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    dsType = Type::value(*ds); // Its type
    SP::SiconosVector residuFree = ds->residuFree();

    // 1 - Lagrangian Non Linear Systems
    if (dsType == Type::LagrangianDS)
    {
      // // residu = M(q*)(v_k,i+1 - v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi) - pi+1

      //       // -- Convert the DS into a Lagrangian one.
      //       SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      //       // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //       SP::SiconosVector qold =d->qMemory()->getSiconosVector(0);
      //       SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      //       SP::SiconosVector q =d->q();


      //       d->computeMass();
      //       SP::SiconosMatrix M = d->mass();
      //       SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      //       //residuFree->zero();


      //       //   std::cout << "(*v-*vold)->norm2()" << (*v-*vold).norm2() << std::endl;

      //       prod(*M, (*v-*vold), *residuFree); // residuFree = M(v - vold)


      //       if (d->fL())  // if fL exists
      //       {
      //         // computes fL(ti,vi,qi)
      //         d->computeFL(told,qold,vold);
      //         double coef = -h*(1-_theta);
      //         // residuFree += coef * fL_i
      //         scal(coef, *d->fL(), *residuFree, false);

      //         // computes fL(ti+1, v_k,i+1, q_k,i+1) = fL(t,v,q)
      //         //d->computeFL(t);
      //         // or  fL(ti+1, v_k,i+1, q(v_k,i+1))
      //         //or
      //         SP::SiconosVector qbasedonv(new SimpleVector(*qold));
      //         *qbasedonv +=  h*( (1-_theta)* *vold + _theta * *v );
      //         d->computeFL(t,qbasedonv,v);
      //         coef = -h*_theta;
      //         // residuFree += coef * fL_k,i+1
      //         scal(coef, *d->fL(), *residuFree, false);
      //       }

      //       if (d->boundaryConditions())
      //       {

      //         d->boundaryConditions()->computePrescribedVelocity(t);

      //         unsigned int columnindex=0;
      //         SP::SimpleMatrix WBoundaryConditions = _WBoundaryConditionsMap[ds];
      //         SP::SiconosVector columntmp(new SimpleVector(ds->getDim()));

      //         for (vector<unsigned int>::iterator  itindex = d->boundaryConditions()->velocityIndices()->begin() ;
      //              itindex != d->boundaryConditions()->velocityIndices()->end();
      //              ++itindex)
      //         {

      //           double DeltaPrescribedVelocity =
      //             d->boundaryConditions()->prescribedVelocity()->getValue(columnindex)
      //             - vold->getValue(columnindex);

      //           WBoundaryConditions->getCol(columnindex,*columntmp);
      //           *residuFree -= *columntmp * (DeltaPrescribedVelocity);

      //           residuFree->setValue(*itindex, columntmp->getValue(*itindex)   *(DeltaPrescribedVelocity));

      //           columnindex ++;

      //         }
      //       }

      //       *(d->workFree())=*residuFree; // copy residuFree in Workfree
      // //      std::cout << "SchatzmanPaoli::ComputeResidu LagrangianDS residufree :"  << std::endl;
      // //      residuFree->display();
      //       if (d->p(1))
      //         *(d->workFree()) -= *d->p(1); // Compute Residu in Workfree Notation !!
      // //      std::cout << "SchatzmanPaoli::ComputeResidu LagrangianDS residu :"  << std::endl;
      // //      d->workFree()->display();
      //         normResidu = d->workFree()->norm2();
      RuntimeException::selfThrow("SchatzmanPaoli::computeResidu - not yet implemented for Dynamical system type: " + dsType);
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
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector q_k = d->qMemory()->getSiconosVector(0); // q_k
      SP::SiconosVector q_k_1 = d->qMemory()->getSiconosVector(1); // q_{k-1}
      SP::SiconosVector v_k = d->velocityMemory()->getSiconosVector(0); //v_k
      // std::cout << "SchatzmanPaoli::computeResidu - q_k_1 =" <<std::endl;
      // q_k_1->display();
      // std::cout << "SchatzmanPaoli::computeResidu - q_k =" <<std::endl;
      // q_k->display();
      // std::cout << "SchatzmanPaoli::computeResidu - v_k =" <<std::endl;
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


      // if (d->boundaryConditions())
      // {
      //   d->boundaryConditions()->computePrescribedVelocity(t);

      //   unsigned int columnindex=0;
      //   SP::SimpleMatrix WBoundaryConditions = _WBoundaryConditionsMap[ds];
      //   SP::SiconosVector columntmp(new SimpleVector(ds->getDim()));

      //   for (vector<unsigned int>::iterator  itindex = d->boundaryConditions()->velocityIndices()->begin() ;
      //        itindex != d->boundaryConditions()->velocityIndices()->end();
      //        ++itindex)
      //   {

      //     double DeltaPrescribedVelocity =
      //       d->boundaryConditions()->prescribedVelocity()->getValue(columnindex)
      //       -vold->getValue(columnindex);

      //     WBoundaryConditions->getCol(columnindex,*columntmp);
      //     *residuFree += *columntmp * (DeltaPrescribedVelocity);

      //     residuFree->setValue(*itindex, - columntmp->getValue(*itindex)   *(DeltaPrescribedVelocity));

      //     columnindex ++;

      //   }
      // }


      // std::cout << "SchatzmanPaoli::ComputeResidu LagrangianLinearTIDS residufree :"  << std::endl;
      // residuFree->display();


      (* d->workFree()) = *residuFree; // copy residuFree in Workfree
      if (d->p(0))
        *(d->workFree()) -= *d->p(0); // Compute Residu in Workfree Notation !!

      // std::cout << "SchatzmanPaoli::ComputeResidu LagrangianLinearTIDS p(0) :"  << std::endl;
      //  if (d->p(0))
      //    d->p(0)->display();
      //  else
      //    std::cout << " p(0) :"  << std::endl;
      // std::cout << "SchatzmanPaoli::ComputeResidu LagrangianLinearTIDS residu :"  << std::endl;
      // d->workFree()->display();



      //     normResidu = d->workFree()->norm2();
      normResidu = 0.0; // we assume that v = vfree + W^(-1) p
      //     normResidu = realresiduFree->norm2();

    }
    else if (dsType == Type::NewtonEulerDS)
    {
      // // residu = M(q*)(v_k,i+1 - v_i) - h*_theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-_theta)*fL(ti,vi,qi) - pi+1

      //     // -- Convert the DS into a Lagrangian one.
      //     SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

      //     // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //     SP::SiconosVector qold =d->qMemory()->getSiconosVector(0);
      //     SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      //     SP::SiconosVector q =d->q();


      //     SP::SiconosMatrix massMatrix = d->massMatrix();
      //     SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      //     prod(*massMatrix, (*v-*vold), *residuFree); // residuFree = M(v - vold)
      //     if (d->fL())  // if fL exists
      //     {
      //       // computes fL(ti,vi,qi)
      //       SP::SiconosVector fLold=d->fLMemory()->getSiconosVector(0);
      //       double _thetaFL=0.5;
      //       double coef = -h*(1-_thetaFL);
      //       // residuFree += coef * fL_i
      //       scal(coef, *fLold, *residuFree, false);
      //       d->computeFL(t);
      // //        printf("cpmputeFreeState d->FL():\n");
      // //   d->fL()->display();
      //       coef = -h*_thetaFL;
      //       scal(coef, *d->fL(), *residuFree, false);
      //     }
      //     *(d->workFree())=*residuFree;
      //     //cout<<"SchatzmanPaoli::computeResidu :\n";
      //     // residuFree->display();
      //     if ( d->p(1) )
      //     *(d->workFree()) -= *d->p(1);
      //     normResidu = d->workFree()->norm2();
      RuntimeException::selfThrow("SchatzmanPaoli::computeResidu - not yet implemented for Dynamical system type: " + dsType);
    }
    else
      RuntimeException::selfThrow("SchatzmanPaoli::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    if (normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void SchatzmanPaoli::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  //double t = simulationLink->nextTime(); // End of the time step
  //double told = simulationLink->startingTime(); // Beginning of the time step
  //double h = t-told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SimpleVector *rold = static_cast<SimpleVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it; // Iterator through the set of DS.

  SP::DynamicalSystem ds; // Current Dynamical System.
  SP::SiconosMatrix W; // W SchatzmanPaoli matrix of the current DS.
  Type::Siconos dsType ; // Type of the current DS.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    dsType = Type::value(*ds); // Its type
    W = WMap[ds]; // Its W SchatzmanPaoli matrix of iteration.

    //1 - Lagrangian Non Linear Systems
    if (dsType == Type::LagrangianDS)
    {
      // // IN to be updated at current time: W, M, q, v, fL
      //       // IN at told: qi,vi, fLi

      //       // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      //       // Index k stands for Newton iteration and thus corresponds to the last computed
      //       // value, ie the one saved in the DynamicalSystem.
      //       // "i" values are saved in memory vectors.

      //       // vFree = v_k,i+1 - W^{-1} ResiduFree
      //       // with
      //       // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi)

      //       // -- Convert the DS into a Lagrangian one.
      //       SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      //       // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //       SP::SiconosVector qold =d->qMemory()->getSiconosVector(0);
      //       SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      //       // --- ResiduFree computation ---
      //       // ResFree = M(v-vold) - h*[theta*fL(t) + (1-theta)*fL(told)]
      //       //
      //       // vFree pointer is used to compute and save ResiduFree in this first step.
      //       SP::SiconosVector vfree = d->workFree();//workX[d];
      //       (*vfree)=*(d->residuFree());

      //       // -- Update W --
      //       // Note: during computeW, mass and jacobians of fL will be computed/
      //       computeW(t,d);
      //       SP::SiconosVector v = d->velocity(); // v = v_k,i+1

      // // -- vfree =  v - W^{-1} ResiduFree --
      //       // At this point vfree = residuFree
      //       // -> Solve WX = vfree and set vfree = X
      //       W->PLUForwardBackwardInPlace(*vfree);
      //       // -> compute real vfree
      //       *vfree *= -1.0;
      //       *vfree += *v;
      RuntimeException::selfThrow("SchatzmanPaoli::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
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
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0); // q_k
      //   SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); //v_k

      // --- ResiduFree computation ---

      // vFree pointer is used to compute and save ResiduFree in this first step.


      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SP::SiconosVector qfree = d->workFree();//workX[d];
      (*qfree) = *(d->residuFree());

      W->PLUForwardBackwardInPlace(*qfree);
      *qfree *= -1.0;
      *qfree += *qold;

    }
    // 3 - Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      // // IN to be updated at current time: W, M, q, v, fL
      // // IN at told: qi,vi, fLi

      // // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // // Index k stands for Newton iteration and thus corresponds to the last computed
      // // value, ie the one saved in the DynamicalSystem.
      // // "i" values are saved in memory vectors.

      // // vFree = v_k,i+1 - W^{-1} ResiduFree
      // // with
      // // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi)

      // // -- Convert the DS into a Lagrangian one.
      // SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);
      // computeW(t,d);
      // // Get state i (previous time step) from Memories -> var. indexed with "Old"
      // SP::SiconosVector qold =d->qMemory()->getSiconosVector(0);
      // SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      // // --- ResiduFree computation ---
      // // ResFree = M(v-vold) - h*[theta*fL(t) + (1-theta)*fL(told)]
      // //
      // // vFree pointer is used to compute and save ResiduFree in this first step.
      // SP::SiconosVector vfree = d->workFree();//workX[d];
      // (*vfree)=*(d->residuFree());
      // //*(d->vPredictor())=*(d->residuFree());

      // // -- Update W --
      // // Note: during computeW, mass and jacobians of fL will be computed/
      // SP::SiconosVector v = d->velocity(); // v = v_k,i+1

      // // -- vfree =  v - W^{-1} ResiduFree --
      // // At this point vfree = residuFree
      // // -> Solve WX = vfree and set vfree = X
      // //   cout<<"SchatzmanPaoli::computeFreeState residu free"<<endl;
      // //   vfree->display();
      // d->luW()->PLUForwardBackwardInPlace(*vfree);
      // //   cout<<"SchatzmanPaoli::computeFreeState -WRfree"<<endl;
      // //   vfree->display();
      // //   scal(h,*vfree,*vfree);
      // // -> compute real vfree
      // *vfree *= -1.0;
      // *vfree += *v;
      RuntimeException::selfThrow("SchatzmanPaoli::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
    }
    else
      RuntimeException::selfThrow("SchatzmanPaoli::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }

}

void SchatzmanPaoli::prepareNewtonIteration(double time)
{
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    computeW(time, *itDS);
  }
}

struct SchatzmanPaoli::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  OneStepNSProblem *osnsp;
  SP::UnitaryRelation UR;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::UnitaryRelation UR) :
    osnsp(p), UR(UR) {};

  void visit(const NewtonImpactNSL& nslaw)
  {
    double e;
    e = nslaw.e();
    Index subCoord(4);
    subCoord[0] = 0;
    subCoord[1] = UR->getNonSmoothLawSize();
    subCoord[2] = 0;
    subCoord[3] = subCoord[1];
    // Only the normal part is multiplied by e

    SP::SiconosVector y_k_1 ;
    /** \warning V.A. 21/09/2011
     * The goal of the following hack is to circumvent a problem
     * when we initialize the interaction. Is it for the moment very hard to get
     * the right pile in the SiconosMemeory for the interaction
     */

    if (UR->yMemory(osnsp->levelMin())->nbVectorsInMemory() > 1)
    {
      y_k_1 = UR->yMemory(osnsp->levelMin(), 1);
    }
    else
    {
      y_k_1 = UR->yMemory(osnsp->levelMin(), 0);
    }
    // std::cout << "y_k_1 " << std::endl;
    // y_k_1->display();
    subscal(e, *y_k_1, *(UR->yp()), subCoord, false);
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    double e;
    e = nslaw.en();
    // Only the normal part is multiplied by e
    SP::SiconosVector y_k_1 ;
    if (UR->yMemory(osnsp->levelMin())->nbVectorsInMemory() > 1)
    {
      y_k_1 = UR->yMemory(osnsp->levelMin(), 1);
    }
    else
    {
      y_k_1 = UR->yMemory(osnsp->levelMin(), 0);
    }
    // std::cout << "y_k_1 " << std::endl;
    // y_k_1->display();

    (*UR->yp())(0) +=  e * (*y_k_1)(0);

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


void SchatzmanPaoli::computeFreeOutput(SP::UnitaryRelation UR, OneStepNSProblem * osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */


  SP::OneStepNSProblems  allOSNS  = simulationLink->oneStepNSProblems();

  // Get relation and non smooth law types
  RELATION::TYPES relationType = UR->getRelationType();
  RELATION::SUBTYPES relationSubType = UR->getRelationSubType();

  SP::DynamicalSystem ds = *(UR->interaction()->dynamicalSystemsBegin());

  unsigned int sizeY = UR->getNonSmoothLawSize();

  unsigned int relativePosition = UR->getRelativePosition();



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
  SP::SiconosVector Xq;
  SP::SiconosVector Yp;
  SP::SiconosVector e;
  SP::SiconosVector Xfree;

  SP::SiconosVector H_alpha;


  // All of these values should be stored in the node corrseponding to the UR when a SchatzmanPaoli scheme is used.
  Xq = UR->xq();
  Yp = UR->yp();

  Xfree = UR->workFree();

  assert(Xfree);


  SP::Interaction mainInteraction = UR->interaction();
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if (relationSubType == LinearTIR)
  {

    if (((*allOSNS)[SICONOS_OSNSP_TS_VELOCITY]).get() != osnsp)
      RuntimeException::selfThrow("SchatzmanPaoli::computeFreeOutput not yet implemented for SICONOS_OSNSP ");

    C = mainInteraction->relation()->C();

    if (C)
    {

      assert(Xfree);
      assert(Yp);
      assert(Xq);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      if (_useGammaForRelation)
      {
        subprod(*C, *Xq, *Yp, coord, true);
      }
      else
      {
        subprod(*C, *Xfree, *Yp, coord, true);
      }
    }
    SP::LagrangianLinearTIR ltir = boost::static_pointer_cast<LagrangianLinearTIR> (mainInteraction->relation());
    e = ltir->e();
    if (e)
    {
      *Yp += *e;
    }

  }
  else
    RuntimeException::selfThrow("SchatzmanPaoli::ComputeFreeOutput not yet implemented  for relation of Type : " + relationType);



  if (UR->getRelationSubType() == LinearTIR)
  {
    SP::SiconosVisitor nslEffectOnFreeOutput(new _NSLEffectOnFreeOutput(osnsp, UR));
    UR->interaction()->nonSmoothLaw()->accept(*nslEffectOnFreeOutput);
  }


}
void SchatzmanPaoli::integrate(double& tinit, double& tend, double& tout, int&)
{
  RuntimeException::selfThrow("SchatzmanPaoli::integrate - not yet implemented :");
}

void SchatzmanPaoli::updateState(unsigned int level)
{

  double h = simulationLink->timeStep();

  const double& RelativeTol = simulationLink->relativeConvergenceTol();
  bool useRCC = simulationLink->useRelativeConvergenceCriteron();
  if (useRCC)
    simulationLink->setRelativeConvergenceCriterionHeld(true);

  DSIterator it;
  SP::SiconosMatrix W;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds];
    // Get the DS type

    Type::Siconos dsType = Type::value(*ds);

    // 1 - Lagrangian Systems
    if (dsType == Type::LagrangianDS || dsType == Type::LagrangianLinearTIDS)
    {
      // get dynamical system
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      //    SiconosVector *vfree = d->velocityFree();
      SP::SiconosVector q = d->q();
      bool baux = dsType == Type::LagrangianDS && useRCC && simulationLink->relativeConvergenceCriterionHeld();
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
        *q +=  * ds->workFree();

      }
      else
        *q =  * ds->workFree();



      // Computation of the velocity

      SP::SiconosVector v = d->velocity();
      SP::SiconosVector q_k_1 = d->qMemory()->getSiconosVector(1); // q_{k-1}

      // std::cout << "SchatzmanPaoli::updateState - q_k_1 =" <<std::endl;
      // q_k_1->display();
      // std::cout << "SchatzmanPaoli::updateState - q =" <<std::endl;
      // q->display();

      *v = 1.0 / (2.0 * h) * (*q - *q_k_1);
      // std::cout << "SchatzmanPaoli::updateState - v =" <<std::endl;
      // v->display();

      // int bc=0;
      // SP::SimpleVector columntmp(new SimpleVector(ds->getDim()));

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
        double aux = ((ds->getWorkVector(DynamicalSystem::local_buffer))->norm2()) / (ds->normRef());
        if (aux > RelativeTol)
          simulationLink->setRelativeConvergenceCriterionHeld(false);
      }

    }
    //2 - Newton Euler Systems
    else if (dsType == Type::NewtonEulerDS)
    {
      //  // get dynamical system
      //       SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);
      //       SP::SiconosVector v = d->velocity();
      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       cout<<"SchatzmanPaoli::updatestate prev v"<<endl;
      //       v->display();
      // #endif

      //       /*d->p has been fill by the Relation->computeInput, it contains
      //            B \lambda _{k+1}*/
      //       *v = *d->p(level); // v = p
      //       d->luW()->PLUForwardBackwardInPlace(*v);

      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       cout<<"SchatzmanPaoli::updatestate hWB lambda"<<endl;
      //       v->display();
      // #endif

      //       *v +=  * ds->workFree();

      // #ifdef SCHATZMANPAOLI_NE_DEBUG
      //       cout<<"SchatzmanPaoli::updatestate work free"<<endl;
      //       ds->workFree()->display();
      //       cout<<"SchatzmanPaoli::updatestate new v"<<endl;
      //       v->display();
      // #endif
      //       //compute q
      //       //first step consists in computing  \dot q.
      //       //second step consists in updating q.
      //       //
      //       SP::SiconosMatrix T = d->T();
      //       SP::SiconosVector dotq = d->dotq();
      //       prod(*T,*v,*dotq,true);
      //       // cout<<"SchatzmanPaoli::updateState v"<<endl;
      //       // v->display();
      //       // cout<<"SchatzmanPaoli::updateState dotq"<<endl;
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
      //       cout<<"new q before normalizing"<<endl;
      //       q->display();
      // #endif

      //       //q[3:6] must be normalized
      //       d->normalizeq();
      //       dotq->setValue(3,(q->getValue(3)-qold->getValue(3))/h);
      //       dotq->setValue(4,(q->getValue(4)-qold->getValue(4))/h);
      //       dotq->setValue(5,(q->getValue(5)-qold->getValue(5))/h);
      //       dotq->setValue(6,(q->getValue(6)-qold->getValue(6))/h);
      //       d->updateT();
      RuntimeException::selfThrow("SchatzmanPaoli::updateState - not yet implemented for Dynamical system type: " + dsType);
    }
    else RuntimeException::selfThrow("SchatzmanPaoli::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}


void SchatzmanPaoli::display()
{
  OneStepIntegrator::display();

  cout << "====== SchatzmanPaoli OSI display ======" << endl;
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    cout << "--------------------------------" << endl;
    cout << "--> W of dynamical system number " << (*it)->number() << ": " << endl;
    if (WMap[*it]) WMap[*it]->display();
    else cout << "-> NULL" << endl;
    cout << "--> and corresponding theta is: " << _theta << endl;
  }
  cout << "================================" << endl;
}


SchatzmanPaoli* SchatzmanPaoli::convert(OneStepIntegrator* osi)
{
  SchatzmanPaoli* moreau = dynamic_cast<SchatzmanPaoli*>(osi);
  return moreau;
}

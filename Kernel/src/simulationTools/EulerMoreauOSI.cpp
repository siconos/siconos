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
#include "EulerMoreauOSI.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "FirstOrderLinearTIDS.hpp"
#include "FirstOrderLinearTIR.hpp"
#include "FirstOrderLinearR.hpp"
#include "FirstOrderType2R.hpp"
#include "FirstOrderType1R.hpp"
#include "NonSmoothLaw.hpp"
#include "CxxStd.hpp"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

using namespace RELATION;

// --- constructor from a minimum set of data ---
EulerMoreauOSI::EulerMoreauOSI(SP::DynamicalSystem newDS, double newTheta) :
  OneStepIntegrator(OSI::EULERMOREAUOSI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  OSIDynamicalSystems->insert(newDS);
  _theta = newTheta;
}

// --- constructor with theta parameter value  ---
EulerMoreauOSI::EulerMoreauOSI(double newTheta):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _theta = newTheta;
}

// --- constructor from a minimum set of data ---
EulerMoreauOSI::EulerMoreauOSI(SP::DynamicalSystem newDS, double newTheta, double newGamma) :
  OneStepIntegrator(OSI::EULERMOREAUOSI), _useGammaForRelation(false)
{
  OSIDynamicalSystems->insert(newDS);
  _theta = newTheta;
  _gamma = newGamma;
  _useGamma = true;
}

// --- constructor from a set of data ---
EulerMoreauOSI::EulerMoreauOSI(double newTheta, double newGamma):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _useGammaForRelation(false)
{
  _theta = newTheta;
  _gamma = newGamma;
  _useGamma = true;
}


// Note: OSIDynamicalSystems and thetaMap must disappear
void EulerMoreauOSI::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
}
const SimpleMatrix EulerMoreauOSI::getW(SP::DynamicalSystem ds)
{
  int dsN = ds->number();
  assert(ds &&
         "EulerMoreauOSI::getW(ds): ds == NULL.");
  //    return *(WMap[0]);
  assert(WMap[dsN] &&
         "EulerMoreauOSI::getW(ds): W[ds] == NULL.");
  return *(WMap[dsN]); // Copy !!
}

SP::SimpleMatrix EulerMoreauOSI::W(SP::DynamicalSystem ds)
{
  assert(ds && "EulerMoreauOSI::W(ds): ds == NULL.");
  //  return WMap[0];
  //  if(WMap[ds]==NULL)
  //    RuntimeException::selfThrow("EulerMoreauOSI::W(ds): W[ds] == NULL.");
  return WMap[ds->number()];
}

void EulerMoreauOSI::setW(const SiconosMatrix& newValue, SP::DynamicalSystem ds)
{
  // Check if ds is in the OSI
  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("EulerMoreauOSI::setW(newVal,ds) - ds does not belong to this Integrator ...");

  // Check dimensions consistency
  unsigned int line = newValue.size(0);
  unsigned int col  = newValue.size(1);

  if (line != col) // Check that newValue is square
    RuntimeException::selfThrow("EulerMoreauOSI::setW(newVal,ds) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::setW(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  unsigned int dsN = ds->number();
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("EulerMoreauOSI::setW(newVal,ds) - unconsistent dimension between newVal and dynamical system to be integrated ");

  // Memory allocation for W, if required
  if (!WMap[dsN]) // allocate a new W if required
  {
    WMap[dsN].reset(new SimpleMatrix(newValue));
  }
  else  // or fill-in an existing one if dimensions are consistent.
  {
    if (line == WMap[dsN]->size(0) && col == WMap[dsN]->size(1))
      *(WMap[dsN]) = newValue;
    else
      RuntimeException::selfThrow("EulerMoreauOSI - setW: inconsistent dimensions with problem size for given input matrix W");
  }
}

void EulerMoreauOSI::setWPtr(SP::SimpleMatrix newPtr, SP::DynamicalSystem ds)
{
  unsigned int line = newPtr->size(0);
  unsigned int col  = newPtr->size(1);
  if (line != col) // Check that newPtr is square
    RuntimeException::selfThrow("EulerMoreauOSI::setWPtr(newVal) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::setWPtr(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("EulerMoreauOSI::setW(newVal) - unconsistent dimension between newVal and dynamical system to be integrated ");

  WMap[ds->number()] = newPtr;                  // link with new pointer
}



const SimpleMatrix EulerMoreauOSI::getWBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds &&
         "EulerMoreauOSI::getWBoundaryConditions(ds): ds == NULL.");
  //    return *(WBoundaryConditionsMap[0]);
  unsigned int dsN = ds->number();
  assert(_WBoundaryConditionsMap[dsN] &&
         "EulerMoreauOSI::getWBoundaryConditions(ds): WBoundaryConditions[ds] == NULL.");
  return *(_WBoundaryConditionsMap[dsN]); // Copy !!
}

SP::SiconosMatrix EulerMoreauOSI::WBoundaryConditions(SP::DynamicalSystem ds)
{
  assert(ds && "EulerMoreauOSI::WBoundaryConditions(ds): ds == NULL.");
  //  return WBoundaryConditionsMap[0];
  //  if(WBoundaryConditionsMap[ds]==NULL)
  //    RuntimeException::selfThrow("EulerMoreauOSI::WBoundaryConditions(ds): W[ds] == NULL.");
  return _WBoundaryConditionsMap[ds->number()];
}


void EulerMoreauOSI::initialize()
{
  OneStepIntegrator::initialize();
  // Get initial time
  double t0 = simulationLink->model()->t0();
  // Compute W(t0) for all ds
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    // Memory allocation for workX. workX[ds*] corresponds to xfree (or vfree in lagrangian case).
    // workX[*itDS].reset(new SiconosVector((*itDS)->getDim()));

    // W initialization
    initW(t0, *itDS);

    //      if ((*itDS)->getType() == Type::LagrangianDS || (*itDS)->getType() == Type::FirstOrderNonLinearDS)
    (*itDS)->allocateWorkVector(DynamicalSystem::local_buffer, WMap[(*itDS)->number()]->size(0));
  }
}
void EulerMoreauOSI::initW(double t, SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix W
  // - insert this matrix into WMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - ds does not belong to the OSI.");
  unsigned int dsN = ds->number();
  if (WMap.find(dsN) != WMap.end())
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - W(ds) is already in the map and has been initialized.");


  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  // Memory allocation for W
  //  WMap[ds].reset(new SimpleMatrix(sizeW,sizeW));
  //   SP::SiconosMatrix W = WMap[ds];

  double h = simulationLink->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  // 1 - First order non linear systems
  if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    //    // Memory allocation for W
    //     WMap[ds].reset(new SimpleMatrix(sizeW,sizeW));
    //     SP::SiconosMatrix W = WMap[ds];

    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = std11::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Copy M or I if M is Null into W


    //    SP::SiconosMatrix W = WMap[ds];

    if (d->M())
      //      *W = *d->M();
      WMap[dsN].reset(new SimpleMatrix(*d->M()));

    else
    {
      //W->eye();
      // WMap[ds].reset(new SimpleMatrix(sizeW,sizeW,Siconos::IDENTITY));
      WMap[dsN].reset(new SimpleMatrix(sizeW, sizeW)); // Warning if the Jacobian is a sparse matrix
      WMap[dsN]->eye();
    }
    SP::SiconosMatrix W = WMap[dsN];


    // d->computeJacobianfx(t); // Computation of JacxF is not required here
    // since it must have been done in OSI->initialize, before a call to this function.

    // Add -h*_theta*jacobian_XF to W
    scal(-h * _theta, *d->jacobianfx(), *W, false);
  }
  else RuntimeException::selfThrow("EulerMoreauOSI::initW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.




}


void EulerMoreauOSI::initWBoundaryConditions(SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix WBoundaryConditions
  // - insert this matrix into WBoundaryConditionsMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initWBoundaryConditions(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("EulerMoreauOSI::initWBoundaryConditions(t,ds) - ds does not belong to the OSI.");

  Type::Siconos dsType = Type::value(*ds);


  RuntimeException::selfThrow("EulerMoreauOSI::initWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void EulerMoreauOSI::computeWBoundaryConditions(SP::DynamicalSystem ds)
{
  // Compute WBoundaryConditions matrix of the Dynamical System ds, at
  // time t and for the current ds state.

  // When this function is called, WBoundaryConditionsMap[ds] is
  // supposed to exist and not to be null Memory allocation has been
  // done during initWBoundaryConditions.

  assert(ds &&
         "EulerMoreauOSI::computeWBoundaryConditions(t,ds) - ds == NULL");

  Type::Siconos dsType = Type::value(*ds);
  //unsigned int dsN = ds->number();
  RuntimeException::selfThrow("EulerMoreauOSI::computeWBoundaryConditions - not yet implemented for Dynamical system type :" + dsType);
}


void EulerMoreauOSI::computeW(double t, SP::DynamicalSystem ds)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, WMap[ds] is supposed to exist and not to be null
  // Memory allocation has been done during initW.

  assert(ds &&
         "EulerMoreauOSI::computeW(t,ds) - ds == NULL");
  unsigned int dsN = ds->number();
  assert((WMap.find(dsN) != WMap.end()) &&
         "EulerMoreauOSI::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

  double h = simulationLink->timeStep();
  Type::Siconos dsType = Type::value(*ds);

  SP::SiconosMatrix W = WMap[dsN];

  // 1 - First order non linear systems
  if (dsType == Type::FirstOrderNonLinearDS)
  {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = std11::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Copy M or I if M is Null into W
    if (d->M())
      *W = *d->M();
    else
      W->eye();

    d->computeJacobianfx(t);
    // Add -h*_theta*jacobian_XF to W
    scal(-h * _theta, *d->jacobianfx(), *W, false);
  }
  // 2 - First order linear systems
  else if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    SP::FirstOrderLinearDS d = std11::static_pointer_cast<FirstOrderLinearDS> (ds);
    if (dsType == Type::FirstOrderLinearDS)
      d->computeA(t);

    if (d->M())
      *W = *d->M();
    else
      W->eye();
    scal(-h * _theta, *d->A(), *W, false);
  }
  else RuntimeException::selfThrow("EulerMoreauOSI::computeW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}



double EulerMoreauOSI::computeResidu()
{
  DEBUG_PRINT("EulerMoreauOSI::computeResidu(), start\n");
  // This function is used to compute the residu for each "EulerMoreauOSI-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.
  //  $\mathcal R(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
  //  $\mathcal R_{free}(x,r) = x - x_{k} -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

  double t = simulationLink->nextTime(); // End of the time step
  double told = simulationLink->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);


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
    SP::SiconosVector residuFree = ds->workspace(DynamicalSystem::freeresidu);
    // 1 - First Order Non Linear Systems
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS)
    {
      // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
      // Residu = Residu - h*r^k_i+1
      //  $\mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
      //  $\mathcal R_{free}(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      SP::FirstOrderNonLinearDS d = std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

      SP::SiconosVector x = d->x(); // last saved value for x
      SP::SiconosMatrix M = d->M();

      *residuFree = *x;
      *residuFree -= *xold;
      //       std::cout<<"EulerMoreauOSI: x"<<endl;
      //       (x)->display();
      //       std::cout<<"EulerMoreauOSI: xold"<<endl;
      //       (xold)->display();



      if (M)
        prod(*M, *residuFree, *residuFree, true);

      if (d->f())
      {

        double coef = -h * (1 - _theta);
        if (dsType == Type::FirstOrderLinearDS)
        {
          // computes f(ti,xi)
          //This computation is done since fold not  is up to date.
          d->computef(told, xold);
          // residuFree += coef * f_i
          scal(coef, *d->f(), *residuFree, false);
        }
        else
        {
          // residuFree += coef * f_i
          scal(coef, *d->fold(), *residuFree, false);
        }
        //          std::cout<<"EulerMoreauOSI: fold"<<endl;
        //          (*d->fold()).display();
        // computes f(ti+1, x_k,i+1) = f(t,x)
        d->computef(t);
        coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        //          std::cout<<"EulerMoreauOSI: f"<<endl;
        //          (*d->f()).display();
        scal(coef, *d->f(), *residuFree, false);
      }
      //      std::cout<<"EulerMoreauOSI: residu free"<<endl;
      //      (*residuFree).display();
      (*(d->workspace(DynamicalSystem::free))) = *residuFree;
      scal(-h, *d->r(), (*d->workspace(DynamicalSystem::free)), false); // residu = residu - h*r
      normResidu = d->workspace(DynamicalSystem::free)->norm2();
      //    std::cout<<"EulerMoreauOSI: residu "<<endl;
      //    (workX[d])->display();
      //    std::cout<<"EulerMoreauOSI: norm residu :"<<normResidu<<endl;


      //(*d->residur())=(*d->r()) -(*d->gAlpha());

      //      std::cout<<"EulerMoreauOSI Type::FirstOrderNonLinearDS: residu r"<<endl;
      //      (*d->residur()).display();
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == Type::FirstOrderLinearTIDS)
    {
      SP::FirstOrderLinearTIDS d = std11::static_pointer_cast<FirstOrderLinearTIDS>(ds);
      //Don't use W because it is LU factorized
      //Residu : R_{free} = M(x^{\alpha}_{k+1} - x_{k}) -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b_{k+1})
      // because x_k+1=x_k:
      //Residu : R_{free} = -hAx_k -hb_{k+1}
      SP::SiconosVector b = d->b();
      if (b)
        *residuFree = *b;
      else
        residuFree->zero();

      // x value at told
      SP::SiconosVector xBuffer = d->workspace(DynamicalSystem::local_buffer);
      *xBuffer = *(d->xMemory()->getSiconosVector(0));
      //    std::cout<<"EulerMoreauOSI TIDS::computeResidu: x_k"<<endl;
      //    xBuffer->display();

      SP::SiconosMatrix A = d->A();
      if (A)
        prod(*A, *xBuffer, *residuFree, false); // residuFree -= -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b_{k+1}


      *residuFree *= -h;

      DEBUG_EXPR(residuFree->display(););

    }
    else
      RuntimeException::selfThrow("EulerMoreauOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    if (normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void EulerMoreauOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.
  DEBUG_PRINT("EulerMoreauOSI::computeFreeState()");

  double t = simulationLink->nextTime(); // End of the time step
  double told = simulationLink->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it; // Iterator through the set of DS.

  SP::DynamicalSystem ds; // Current Dynamical System.
  SP::SiconosMatrix W; // W EulerMoreauOSI matrix of the current DS.
  Type::Siconos dsType ; // Type of the current DS.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    dsType = Type::value(*ds); // Its type
    W = WMap[ds->number()]; // Its W EulerMoreauOSI matrix of iteration.

    // 1 - First Order Non Linear Systems
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      // xFree = x_k,i+1  - [W_k,i+1]^{-1} * ResiduFree_k,i+1
      // with ResiduFree_k,i+1 = = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      // IN to be updated at current time: W, f
      // IN at told: f
      // IN, not time dependant: M
      SP::FirstOrderNonLinearDS d = std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //    SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

      // --- ResiduFree computation ---
      // ResiduFree = M(x-xold) - h*[theta*f(t) + (1-theta)*f(told)]
      //
      // xFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector xfree = d->workspace(DynamicalSystem::free);//workX[d];
      *xfree = *(d->workspace(DynamicalSystem::freeresidu));

      if (_useGamma)
      {
        SP::SiconosVector rold = d->rMemory()->getSiconosVector(0);
        double coeff = -h * (1 - _gamma);
        scal(coeff, *rold, *xfree, false); //  residuFree += h(1-gamma)*rold
      }

      SP::SiconosVector x = d->x(); // last saved value for x

      // -- xfree =  x - W^{-1} ResiduFree --
      // At this point xfree = residuFree
      // -> Solve WX = xfree and set xfree = X
      // -- Update W --
      if (dsType != Type::FirstOrderLinearTIDS)
        computeW(t, d);

      W->PLUForwardBackwardInPlace(*xfree);

      // -> compute real xfree
      *xfree *= -1.0;
      *xfree += *x;
      //    std::cout<<" moreau::computefreestate xfree"<<endl;
      //    xfree->display();

      //       if (!simulationLink->model()->nonSmoothDynamicalSystem()->isLinear())
      //       {
      SP::SiconosVector xp = d->xp();
      //      std::cout<<"before moreau::computefreestate xp"<<endl;
      //      xp->display();
      W->PLUForwardBackwardInPlace(*xp);
      scal(h, *xp, *xp);
      *xp += *xfree;
      //      std::cout<<"after moreau::computefreestate xp"<<endl;
      //      xp->display();
      SP::SiconosVector xq = d->xq();
      *xq = *xp;
      *xq -= *x;

      //          std::cout <<boolalpha << _useGamma << std::endl;
      //          std::cout <<boolalpha << _useGammaForRelation << std::endl;
      //          std::cout <<_gamma << std::endl;

      if (_useGammaForRelation)
      {


        // This part is ony valid for the linear case. See comments in Dev Notes.

        if (dsType != Type::FirstOrderLinearDS && dsType != Type::FirstOrderLinearTIDS)
        {
          RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState with _useGammaForRelation  - not yet implemented for Dynamical system type: " + dsType);
        }
          *xq = *xfree;
        //            std::cout << "xq before" << std::endl;
        //           xq->display();

        scal(_gamma, *xq, *xq);
        SP::SiconosVector xold = d->xMemory()->getSiconosVector(0);
        //            std::cout << "xold" << std::endl;
        //           xold->display();

        scal(1.0 - _gamma, *xold, *xq, false);
        //           std::cout << "xq after" << std::endl;
        //           xq->display();

      }
      DEBUG_EXPR(xfree->display(););
      DEBUG_EXPR(xp->display(););
      DEBUG_EXPR(xq->display(););


      //      }

    }
    else
      RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }

}

void EulerMoreauOSI::prepareNewtonIteration(double time)
{
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    computeW(time, *itDS);
  }
}


struct EulerMoreauOSI::_NSLEffectOnFreeOutput : public SiconosVisitor
{
  using SiconosVisitor::visit;

  OneStepNSProblem * _osnsp;
  SP::Interaction _inter;

  _NSLEffectOnFreeOutput(OneStepNSProblem *p, SP::Interaction inter) :
    _osnsp(p), _inter(inter) {};

  void visit(const EqualityConditionNSL& nslaw)
  {
    ;
  }
  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    ;
  }
};


void EulerMoreauOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */

  SP::OneStepNSProblems  allOSNS  = simulationLink->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

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
  SP::BlockVector Xq;
  SP::SiconosVector Yp;
  SP::BlockVector Xfree;

  SP::SiconosVector H_alpha;


  /** \todo VA. All of these values should be stored in a node in the interactionGraph
   * corrseponding to the Interaction
   * when a EulerMoreauOSI scheme is used.
   */

  Xq = inter->dataXq();
  Yp = inter->yp();

  if (relationType == FirstOrder)
  {
    Xfree = inter->data(FirstOrderR::free);
  }

  assert(Xfree);


  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if (relationType == FirstOrder && relationSubType == Type2R)
  {


    SP::SiconosVector lambda;
    lambda = inter->lambda(0);
    FirstOrderType2R& rel = *std11::static_pointer_cast<FirstOrderType2R>(mainInteraction->relation());
    C = rel.C();
    D = rel.D();
    assert(lambda);

    if (D)
    {
      coord[3] = D->size(1);
      coord[5] = D->size(1);
      subprod(*D, *lambda, *Yp, coord, true);

      *Yp *= -1.0;
    }
    if (C)
    {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xq, *Yp, coord, false);

    }

    if (_useGammaForRelation)
    {
      RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
    }
    H_alpha = inter->Halpha();
    assert(H_alpha);
    *Yp += *H_alpha;
  }

  else if (relationType == FirstOrder && relationSubType == Type1R)
  {
    FirstOrderType1R& rel = *std11::static_pointer_cast<FirstOrderType1R>(mainInteraction->relation());
    C = rel.C();
    F = rel.F();
    assert(Xfree);
    assert(Xq);

    if (F)
    {
      coord[3] = F->size(1);
      coord[5] = F->size(1);
      subprod(*F, *inter->dataZ(), *Yp, coord, true);

    }
    if (C)
    {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, *Yp, coord, false);

    }

    if (_useGammaForRelation)
    {
      RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
    }
    H_alpha = inter->Halpha();
    assert(H_alpha);
    *Yp += *H_alpha;
  }
  else
  {
    C = mainInteraction->relation()->C();

    if (C)
    {

      assert(Xfree);
      assert(Xq);

      coord[3] = C->size(1);
      coord[5] = C->size(1);
      // creates a POINTER link between workX[ds] (xfree) and the
      // corresponding interactionBlock in each Interactionfor each ds of the
      // current Interaction.
      if (_useGammaForRelation)
      {
        subprod(*C, *Xq, *Yp, coord, true);
      }
      else
      {
        subprod(*C, *Xfree, *Yp, coord, true);
      }
    }

    if (relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
    {
      // In the first order linear case it may be required to add e + FZ to q.
      // q = HXfree + e + FZ
      SP::SiconosVector e;
      if (relationSubType == LinearTIR)
      {
        e = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->e();
        F = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->F();
      }
      else
      {
        e = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->e();
        F = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->F();
      }

      if (e)
        *Yp += *e;

      if (F)
      {
        coord[3] = F->size(1);
        coord[5] = F->size(1);
        subprod(*F, *inter->dataZ(), *Yp, coord, false);
      }
    }

  }

 
}
void EulerMoreauOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  // Last parameter is not used (required for LsodarOSI but not for EulerMoreauOSI).

  //double h = tend - tinit;
  tout = tend;

  DSIterator it;
  SP::SiconosMatrix W;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds->number()];
    Type::Siconos dsType = Type::value(*ds);
    RuntimeException::selfThrow("EulerMoreauOSI::integrate - not yet implemented for Dynamical system type :" + dsType);
  }
}

void EulerMoreauOSI::updateState(const unsigned int level)
{

  DEBUG_PRINT("EulerMoreauOSI::updateState(const unsigned int level)\n");

  double h = simulationLink->timeStep();

  double RelativeTol = simulationLink->relativeConvergenceTol();
  bool useRCC = simulationLink->useRelativeConvergenceCriteron();
  if (useRCC)
    simulationLink->setRelativeConvergenceCriterionHeld(true);

  DSIterator it;
  SP::SiconosMatrix W;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds->number()];
    // Get the DS type

    Type::Siconos dsType = Type::value(*ds);

    // 1 - First Order Systems
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      SP::FirstOrderNonLinearDS fonlds = std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      SP::SiconosVector x = ds->x();
      bool baux = (useRCC && dsType == Type::FirstOrderNonLinearDS && simulationLink->relativeConvergenceCriterionHeld());
      if (level != LEVELMAX)
      {

        //    SP::SiconosVector xFree = fonlds->xFree();

        // Save value of q in local_buffer for relative convergence computation
        if (baux)
          ds->addWorkVector(x, DynamicalSystem::local_buffer);

        //        std::cout <<boolalpha << _useGamma << std::endl;
        //        std::cout <<_gamma << std::endl;
        if (_useGamma)
        {
          //SP::SiconosVector rold =d->rMemory()->getSiconosVector(0);
          // Solve W(x-xfree) = hr
          scal(_gamma * h, *fonlds->r(), *x); // x = gamma*h*r
          // scal((1.0-_gamma)*h,*rold,*x,false)// x += (1-gamma)*h*rold
        }
        else
        {
          // Solve W(x-xfree) = hr
          scal(h, *fonlds->r(), *x); // x = h*r
          //      scal(h,*fonlds->gAlpha(),*x); // x = h*gApha
        }

        W->PLUForwardBackwardInPlace(*x); // x =h* W^{-1} *r

        *x += *(fonlds->workspace(DynamicalSystem::free)); //*workX[ds]; // x+=xfree
      }
      else
      {
        *x = *(fonlds->workspace(DynamicalSystem::free)); //*workX[ds]; // x=xfree
      }

      if (baux)
      {
        ds->subWorkVector(x, DynamicalSystem::local_buffer);
        double aux = ((ds->workspace(DynamicalSystem::local_buffer))->norm2()) / (ds->normRef());
        if (aux > RelativeTol)
          simulationLink->setRelativeConvergenceCriterionHeld(false);
      }


      //  }else if (dsType == Type::FirstOrderLinearTIDS){
      //    SP::FirstOrderNonLinearDS fonlds = std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      //    SP::SiconosVector x = ds->x();
      //    // Solve W(x-xfree) = hr
      //    *x=*fonlds->r();
      //    W->PLUForwardBackwardInPlace(*x); // x = W^{-1} *r
      //    scal(h,*x,*x); // x = h*W^{-1}*r
      //    *x +=*(fonlds->xfree());//*workX[ds]; // x+=xfree
      //    //    std::cout<<"X alpha+1"<<endl;
      //    //    x->display();
    }
    else RuntimeException::selfThrow("EulerMoreauOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}


// bool EulerMoreauOSI::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
// {
//   DEBUG_PRINT("addInteractionInIndexSet(SP::Interaction inter, unsigned int i)\n");

//   assert(i == 1);
//   double h = simulationLink->timeStep();
//   double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
//   double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity

//   double gamma = 1.0 / 2.0;
//   if (_useGamma)
//   {
//     gamma = _gamma;
//   }
//   DEBUG_PRINTF("EulerMoreauOSI::addInteractionInIndexSet of level = %i yref=%e, yDot=%e, y_estimated=%e.\n", i,  y, yDot, y + gamma * h * yDot);
//   y += gamma * h * yDot;
//   assert(!isnan(y));
//   DEBUG_EXPR(
//     if (y <= 0)
//     DEBUG_PRINT("EulerMoreauOSI::addInteractionInIndexSet ACTIVATE.\n");
//   );
//   return (y <= 0.0);
// }


// bool EulerMoreauOSI::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
// {
//   assert(i == 1);
//   double h = simulationLink->timeStep();
//   double y = (inter->y(i - 1))->getValue(0); // for i=1 y(i-1) is the position
//   double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity
//   double gamma = 1.0 / 2.0;
//   if (_useGamma)
//   {
//     gamma = _gamma;
//   }
//   DEBUG_PRINTF("EulerMoreauOSI::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y + gamma * h * yDot);
//   y += gamma * h * yDot;
//   assert(!isnan(y));

//   DEBUG_EXPR(
//     if (y > 0)
//     DEBUG_PRINT("EulerMoreauOSI::removeInteractionInIndexSet DEACTIVATE.\n");
//   );
//   return (y > 0.0);
// }


void EulerMoreauOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== EulerMoreauOSI OSI display ======" <<std::endl;
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    std::cout << "--------------------------------" <<std::endl;
    std::cout << "--> W of dynamical system number " << (*it)->number() << ": " <<std::endl;
    if (WMap[(*it)->number()]) WMap[(*it)->number()]->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "--> and corresponding theta is: " << _theta <<std::endl;
  }
  std::cout << "================================" <<std::endl;
}

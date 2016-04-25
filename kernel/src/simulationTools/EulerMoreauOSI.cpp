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
#include "OneStepNSProblem.hpp"

//#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>

using namespace RELATION;

// --- constructor with theta parameter value  ---
EulerMoreauOSI::EulerMoreauOSI(double theta):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _gamma(1.0), _useGamma(false), _useGammaForRelation(false)
{
  _theta = theta;
}

// --- constructor from a set of data ---
EulerMoreauOSI::EulerMoreauOSI(double theta, double gamma):
  OneStepIntegrator(OSI::EULERMOREAUOSI), _useGammaForRelation(false)
{
  _theta = theta;
  _gamma = gamma;
  _useGamma = true;
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
  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - ds does not belong to the OSI.");


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


void EulerMoreauOSI::initialize(Model& m)
{
  OneStepIntegrator::initialize(m);
  // Get initial time
  double t0 = _simulation->startingTime();

  std::cout << std::endl;
  // Compute W(t0) for all ds
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    // W initialization
    initW(t0, ds);
    ds->allocateWorkVector(DynamicalSystem::local_buffer, WMap[ds->number()]->size(0));
  }
}
void EulerMoreauOSI::initW(double t, SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix W
  // - insert this matrix into WMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - ds == NULL");

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - ds does not belong to the OSI.");

  unsigned int dsN = ds->number();
  if (WMap.find(dsN) != WMap.end())
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - W(ds) is already in the map and has been initialized.");


  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  // Memory allocation for W
  //  WMap[ds].reset(new SimpleMatrix(sizeW,sizeW));
  //   SP::SiconosMatrix W = WMap[ds];

  double h = _simulation->timeStep();
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

  if (!(checkOSI(_dynamicalSystemsGraph->descriptor(ds))))
    RuntimeException::selfThrow("EulerMoreauOSI::initW(t,ds) - ds does not belong to the OSI.");

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


void EulerMoreauOSI::computeW(double t, DynamicalSystem& ds, DynamicalSystemsGraph::VDescriptor& dsgVD)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, WMap[ds] is supposed to exist and not to be null
  // Memory allocation has been done during initW.

  unsigned int dsN = ds.number();
  assert((WMap.find(dsN) != WMap.end()) &&
         "EulerMoreauOSI::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

  double h = _simulation->timeStep();
  Type::Siconos dsType = Type::value(ds);

  SiconosMatrix& W = *WMap[dsN];

  // 1 - First order non linear systems
  if (dsType == Type::FirstOrderNonLinearDS)
  {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    FirstOrderNonLinearDS& d = static_cast<FirstOrderNonLinearDS&> (ds);

    // Copy M or I if M is Null into W
    if (d.M())
      W = *d.M();
    else
      W.eye();

    d.computeJacobianfx(t);
    // Add -h*_theta*jacobianfx to W
    scal(-h * _theta, *d.jacobianfx(), W, false);
  }
  // 2 - First order linear systems
  else if (dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
  {
    FirstOrderLinearDS& d = static_cast<FirstOrderLinearDS&> (ds);
    if (dsType == Type::FirstOrderLinearDS)
      d.computeA(t);

    if (d.M())
      W = *d.M();
    else
      W.eye();

    scal(-h * _theta, *d.A(), W, false);
  }
  else RuntimeException::selfThrow("EulerMoreauOSI::computeW - not yet implemented for Dynamical system type :" + dsType);

//  if (_useGamma)
  {
    Topology& topo = *_simulation->nonSmoothDynamicalSystem()->topology();
    DynamicalSystemsGraph& DSG0 = *topo.dSG(0);
    InteractionsGraph& indexSet = *topo.indexSet(0);
    DynamicalSystemsGraph::OEIterator oei, oeiend;
    InteractionsGraph::VDescriptor ivd;
    SP::SiconosMatrix K;
    SP::Interaction inter;
    for (std11::tie(oei, oeiend) = DSG0.out_edges(dsgVD); oei != oeiend; ++oei)
    {
      inter = DSG0.bundle(*oei);
      ivd = indexSet.descriptor(inter);
      FirstOrderR& rel = static_cast<FirstOrderR&>(*inter->relation());
      K = rel.K();
      if (!K) K = (*indexSet.properties(ivd).workMatrices)[FirstOrderR::mat_K];
      if (K)
      {
        scal(-h * _gamma, *K, W, false);
      }
    }
  }
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

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  DEBUG_PRINTF("nextTime %f\n", t);
  DEBUG_PRINTF("startingTime %f\n", told);
  DEBUG_PRINTF("time step size %f\n", h);


  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  // XXX TMP hack -- xhub
  Topology& topo = *_simulation->nonSmoothDynamicalSystem()->topology();
  DynamicalSystemsGraph& DSG0 = *topo.dSG(0);


  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds); // Its type
    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
    // Maurice will do that with subgraph :)
    DynamicalSystemsGraph::VDescriptor dsgVD = topo.getDSG0Descriptor(ds);
    VectorOfVectors& workVectors = *DSG0.properties(dsgVD).workVectors;
    SiconosVector& residuFree = *workVectors[FirstOrderDS::residuFree];
    SiconosVector& residu = *workVectors[FirstOrderDS::residu];
    // 1 - First Order Non Linear Systems
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS)
    {
      // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
      //  $\mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
      //  $\mathcal R_{free}(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

      // Note: indices k/k+1 corresponds to value at the beginning/end of the time step.
      // Newton iterate are x and r

      FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // 1 - Compute the free residu (purely on the "smooth" dynamics)

      residuFree = *d.x(); // last saved value for x: could be x_k or x_{k+1}^alpha
      SiconosVector& xold = *d.xMemory()->getSiconosVector(0);
      residuFree -= xold; // state x_k (at previous time step)

      SP::SiconosMatrix M = d.M();
      if (M)
        prod(*M, residuFree, residuFree, true);
      // at this step, we have residuFree = M(x - x_k)
      DEBUG_PRINT("EulerMoreauOSI::computeResidu residuFree = M(x - x_k)\n");
      DEBUG_EXPR(residuFree.display());

      if (d.f())
      {
        double coef = -h * (1 - _theta);
        if (dsType == Type::FirstOrderLinearDS)
        {
          // computes f(t_k,x_k)
          //This computation is done since fold not  is up to date.
          d.computef(told, xold);
          // residuFree += -h * (1 - _theta) * f(t_k,x_k)
          scal(coef, *d.f(), residuFree, false);
        }
        else // FirstOrderNonLinearDS
        {
          // residuFree += -h * (1 - _theta) * f(t_k,x_k)
          scal(coef, *d.fold(), residuFree, false);
        }

         // computes f(t_{x+1}, x_{k+1}^alpha)
        d.computef(t);
        coef = -h * _theta;
        // residuFree += -h * _theta * f(t_{x+1}, x_{k+1}^alpha)
        scal(coef, *d.f(), residuFree, false);
      }

      // now compute the residu = residuFree - h*gamma*r - h*(1-gamma)r_k
      residu = residuFree;

      if (!_useGamma) // no gamma
      {
        scal(-h, *d.r(), residu, false); // residu = residu - h*r
      }
      else
      {
        scal(-h*_gamma, *d.r(), residu, false);
        scal(-h*(1-_gamma), *d.rMemory()->getSiconosVector(0), residu, false);
      }

      normResidu = residu.norm2();
      DEBUG_EXPR(residu.display());
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderLinearTIDS& d = *std11::static_pointer_cast<FirstOrderLinearTIDS>(ds);
      //Don't use W because it is LU factorized
      //Residu : R_{free} = M(x^{\alpha}_{k+1} - x_{k}) -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b_{k+1})
      SP::SiconosVector b = d.b();
      if (b)
        residuFree = *b;
      else
        residuFree.zero();

      SiconosVector& xBuffer = *workVectors[FirstOrderDS::xBuffer];
      SP::SiconosMatrix A = d.A();

      if (A) // residuFree += -h( A (\theta x_{k+1}^{\alpha} + (1-\theta) x_k)
      {

        prod(*A, *d.xMemory()->getSiconosVector(0), xBuffer, true);
        double coef = -h * (1 - _theta);
        scal(coef, xBuffer, residuFree, false);

        prod(*A, *d.x(), xBuffer, true);
        coef = -h * _theta;
        scal(coef, xBuffer, residuFree, false);
      }

      // final touch, residuFree += M(x_{k+1}^{\alpha} - x_k)
      xBuffer = *d.x() - *d.xMemory()->getSiconosVector(0);
      SP::SiconosMatrix M = d.M();
      if (M)
      {
        prod(*M, xBuffer, residuFree, false);
      }
      else
      {
        residuFree += xBuffer;
      }


    }
    else
      RuntimeException::selfThrow("EulerMoreauOSI::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    DEBUG_PRINT("EulerMoreauOSI::computeResidu final residuFree\n");
    DEBUG_EXPR(residuFree.display());


    if (normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void EulerMoreauOSI::computeFreeState()
{
  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.
  DEBUG_PRINT("EulerMoreauOSI::computeFreeState() starts\n");

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d.rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //

  SP::DynamicalSystem ds; // Current Dynamical System.
  Type::Siconos dsType ; // Type of the current DS.


  
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    ds = _dynamicalSystemsGraph->bundle(*dsi);

    // XXX TMP hack -- xhub
    // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
    // Maurice will do that with subgraph :)
    
    DynamicalSystemsGraph::VDescriptor dsgVD = _dynamicalSystemsGraph->descriptor(ds);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(dsgVD).workVectors;


    dsType = Type::value(*ds); // Its type
    SiconosMatrix& W = *WMap[ds->number()]; // Its W EulerMoreauOSI matrix of iteration.

    // 1 - First Order Non Linear Systems
    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      // xfree =  x - W^{-1} (ResiduFree - h(1-gamma)*rold)
      // with ResiduFree = = M(x - x_k) - h*theta*f(t_{k+1}, x) - h*(1-theta)*f(t_k, x_k)

      // to be updated at current time: W, f
      // fold is f at t_k
      // not time dependant: M
      FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      //    SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

      SiconosVector& x = *d.x(); // x = x_k or x = x_{k+1}^{\alpha}
      // xfree gets ResiduFree at first
      SiconosVector& xfree = *workVectors[FirstOrderDS::xfree];
      xfree = *workVectors[FirstOrderDS::residuFree];

      DEBUG_PRINT("EulerMoreauOSI::computeFreeState xfree <- residuFree\n");
      DEBUG_EXPR(xfree.display());

      if (_useGamma)
      {
        SiconosVector& rold = *d.rMemory()->getSiconosVector(0);
        double coeff = -h * (1 - _gamma);
        scal(coeff, rold, xfree, false); //  xfree += -h(1-gamma)*rold
      }


      // At this point xfree = (ResiduFree - h(1-gamma)*rold)
      // -> Solve WX = xfree and set xfree = X
      W.PLUForwardBackwardInPlace(xfree);

      // at this point, xfree = W^{-1} (ResiduFree - h(1-gamma)*rold)
      // -> compute real xfree = x - W^{-1} (ResiduFree - h(1-gamma)*rold)
      xfree *= -1.0;
      xfree += x;

      DEBUG_EXPR(xfree.display());

      // now the crazy intermediate variables
      // xPartialNS was updated before this fonction call
      // It constains either 0 (first Newton iterate)
      // or g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
      SiconosVector& xPartialNS = *workVectors[FirstOrderDS::xPartialNS];
      DEBUG_PRINT("EulerMoreauOSI::computeFreeState xPartialNS from Interaction\n");
      DEBUG_EXPR(xPartialNS.display());

      // -> Solve WX = g(x, \lambda, t_{k+1}) - B_{k+1}^{\alpha} \lambda - K_{k+1}^{\alpha} x
      // and set xPartialNS = X
      W.PLUForwardBackwardInPlace(xPartialNS);
      scal(h, xPartialNS, xPartialNS);

      // compute real xPartialNS = xfree + ...
      xPartialNS += xfree;
      DEBUG_PRINT("EulerMoreauOSI::computeFreeState xPartialNS real value\n");
      DEBUG_EXPR(xPartialNS.display());

      // deltaxForRelation = (\widetilde{K}_{k+1}^{\alpha})^{-1} xPartialNS - x
      SiconosVector& deltaxForRelation = *workVectors[FirstOrderDS::deltaxForRelation];
      deltaxForRelation = xPartialNS;

      deltaxForRelation -= x;

      DEBUG_EXPR(deltaxForRelation.display());

      // have a look at the end of the DevNotes for this part
      if (_useGammaForRelation)
      {
        if (!(dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS))
          RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - _useGammaForRelation == true is only implemented for FirstOrderLinearDS or FirstOrderLinearTIDS");

        deltaxForRelation = xfree;

        scal(_gamma, deltaxForRelation, deltaxForRelation);
        SiconosVector& xold = *d.xMemory()->getSiconosVector(0);

        scal(1.0 - _gamma, xold, deltaxForRelation, false);
      }

      // some output
      DEBUG_EXPR(xfree.display(););
      DEBUG_EXPR(xPartialNS.display(););
      DEBUG_EXPR(deltaxForRelation.display(););

    }
    else
      RuntimeException::selfThrow("EulerMoreauOSI::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }
}

void EulerMoreauOSI::prepareNewtonIteration(double time)
{
  // XXX TMP hack -- xhub
  // we have to iterate over the edges of the DSG0 -> the following won't be necessary anymore
  // Maurice will do that with subgraph :)
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    DynamicalSystemsGraph::VDescriptor dsgVD = _dynamicalSystemsGraph->descriptor(ds);
    computeW(time, *ds, dsgVD);
  }
}

/// @cond

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

/// @endcond

void EulerMoreauOSI::computeFreeOutput(InteractionsGraph::VDescriptor& vertex_inter, OneStepNSProblem* osnsp)
{
  /** \warning: ensures that it can also work with two different osi for two different ds ?
   */

  SP::OneStepNSProblems  allOSNS  = _simulation->oneStepNSProblems();
  SP::InteractionsGraph indexSet = osnsp->simulation()->indexSet(osnsp->indexSetLevel());
  SP::Interaction inter = indexSet->bundle(vertex_inter);

  VectorOfBlockVectors& DSlink = *indexSet->properties(vertex_inter).DSlink;
  VectorOfVectors& workV = *indexSet->properties(vertex_inter).workVectors;
  VectorOfSMatrices& workM = *indexSet->properties(vertex_inter).workMatrices;
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
  SP::BlockVector Xfree;

  SP::SiconosVector H_alpha;

  deltax = DSlink[FirstOrderR::deltax];
  SiconosVector& yForNSsolver = *inter->yForNSsolver();

  Xfree = DSlink[FirstOrderR::xfree];

  assert(Xfree);


  SP::Interaction mainInteraction = inter;
  assert(mainInteraction);
  assert(mainInteraction->relation());

  if (relationType == FirstOrder && (relationSubType == Type2R || relationSubType == NonLinearR))
  {
    SiconosVector& lambda = *inter->lambda(0);
    FirstOrderR& rel = *std11::static_pointer_cast<FirstOrderR>(mainInteraction->relation());
    C = rel.C();
    if (!C) C = workM[FirstOrderR::mat_C];
    D = rel.D();
    if (!D) D = workM[FirstOrderR::mat_D];

    if (D)
    {
      coord[3] = D->size(1);
      coord[5] = D->size(1);
      subprod(*D, lambda, yForNSsolver, coord, true);

      yForNSsolver *= -1.0;
    }
    else
    {
      subscal(0, yForNSsolver, yForNSsolver, coord, true);
    }

    if (C)
    {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *deltax, yForNSsolver, coord, false);

    }

    if (_useGammaForRelation)
    {
      RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Type2R and H_alpha->getValue() should return the mid-point value");
    }

    H_alpha = inter->Halpha();
    assert(H_alpha);
    yForNSsolver += *H_alpha;
  }
  else if (relationType == FirstOrder && relationSubType == Type1R)
  {
    FirstOrderType1R& rel = *std11::static_pointer_cast<FirstOrderType1R>(mainInteraction->relation());
    C = rel.C();
    if (!C) C = workM[FirstOrderR::mat_C];
    F = rel.F();
    if (!F) F = workM[FirstOrderR::mat_F];
    assert(Xfree);
    assert(deltax);

    if (F)
    {
      coord[3] = F->size(1);
      coord[5] = F->size(1);
      subprod(*F, *DSlink[FirstOrderR::z], yForNSsolver, coord, true);

    }
    if (C)
    {
      coord[3] = C->size(1);
      coord[5] = C->size(1);
      subprod(*C, *Xfree, yForNSsolver, coord, false);

    }

    if (_useGammaForRelation)
    {
      RuntimeException::selfThrow("EulerMoreauOSI::ComputeFreeOutput not yet implemented with useGammaForRelation() for FirstorderR and Typ2R and H_alpha->getValue() should return the mid-point value");
    }
    H_alpha = inter->Halpha();
    if(H_alpha)
    {
      yForNSsolver += *H_alpha;
    }
  }
  else // First Order Linear Relation
  {
    C = mainInteraction->relation()->C();
    if (!C) C = workM[FirstOrderR::mat_C];

    if (C)
    {

      assert(Xfree);
      assert(deltax);

      coord[3] = C->size(1);
      coord[5] = C->size(1);

      if (_useGammaForRelation)
      {
        subprod(*C, *deltax, yForNSsolver, coord, true);
      }
      else
      {
        subprod(*C, *Xfree, yForNSsolver, coord, true);
      }
    }

    if (relationType == FirstOrder && (relationSubType == LinearTIR || relationSubType == LinearR))
    {
      // In the first order linear case it may be required to add e + FZ to y.
      // y = CXfree + e + FZ
      SP::SiconosVector e;
      if (relationSubType == LinearTIR)
      {
        e = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->e();
        F = std11::static_pointer_cast<FirstOrderLinearTIR>(mainInteraction->relation())->F();
      }
      else
      {
        e = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->e();
        if (!e) e = workV[FirstOrderR::e];
        F = std11::static_pointer_cast<FirstOrderLinearR>(mainInteraction->relation())->F();
        if (!F) F = workM[FirstOrderR::mat_F];
      }

      if (e)
        yForNSsolver += *e;

      if (F)
      {
        coord[3] = F->size(1);
        coord[5] = F->size(1);
        subprod(*F, *DSlink[FirstOrderR::z], yForNSsolver, coord, false);
      }
    }

  }
}

void EulerMoreauOSI::integrate(double& tinit, double& tend, double& tout, int&)
{
  // Last parameter is not used (required for LsodarOSI but not for EulerMoreauOSI).

  //double h = tend - tinit;
  tout = tend;


  SP::SiconosMatrix W;
  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    W = WMap[ds->number()];
    Type::Siconos dsType = Type::value(*ds);
    RuntimeException::selfThrow("EulerMoreauOSI::integrate - not yet implemented for Dynamical system type :" + dsType);
  }
}

void EulerMoreauOSI::updateState(const unsigned int level)
{

  DEBUG_PRINT("EulerMoreauOSI::updateState\n");

  double h = _simulation->timeStep();

  double RelativeTol = _simulation->relativeConvergenceTol();
  bool useRCC = _simulation->useRelativeConvergenceCriteron();
  if (useRCC)
    _simulation->setRelativeConvergenceCriterionHeld(true);

  DynamicalSystemsGraph::VIterator dsi, dsend;

  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);

    SiconosMatrix& W = *WMap[ds->number()];

    // Get the DS type
    Type::Siconos dsType = Type::value(*ds);
    DynamicalSystemsGraph::VDescriptor dsgVD = _dynamicalSystemsGraph->descriptor(ds);
    VectorOfVectors& workVectors = *_dynamicalSystemsGraph->properties(dsgVD).workVectors;

    if (dsType == Type::FirstOrderNonLinearDS || dsType == Type::FirstOrderLinearDS || dsType == Type::FirstOrderLinearTIDS)
    {
      FirstOrderNonLinearDS& d = *std11::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      SiconosVector& x = *ds->x();
      DEBUG_PRINT("EulerMoreauOSI::updateState Old value of x\n");
      DEBUG_EXPR(x.display());
      DEBUG_PRINT("EulerMoreauOSI::updateState residu value\n");
      DEBUG_EXPR(d.r()->display());

      // TODO ???
      bool baux = (useRCC && dsType == Type::FirstOrderNonLinearDS && _simulation->relativeConvergenceCriterionHeld());
      if (level != LEVELMAX)
      {

        //    SP::SiconosVector xFree = d->xFree();

        // Save value of q in local_buffer for relative convergence computation
        if (baux)
          *workVectors[FirstOrderDS::xBuffer] = x;

        //        std::cout <<boolalpha << _useGamma << std::endl;
        //        std::cout <<_gamma << std::endl;
        if (_useGamma)
        {
          // XXX UseGamma broken ? -- xhub
          scal(_gamma * h, *d.r(), x); // x = gamma*h*r
        }
        else
        {
          scal(h, *d.r(), x); // x = h*r
        }

        W.PLUForwardBackwardInPlace(x); // x = h* W^{-1} *r

        x += *workVectors[FirstOrderDS::xfree]; // x+=xfree
      }
      else
      {
        RuntimeException::selfThrow("EulerMoreauOSI::updateState - level != LEVELMAX is not supposed to happen !");
        x = *workVectors[FirstOrderDS::xfree]; // x = xfree
      }

      if (baux)
      {
        *workVectors[FirstOrderDS::xBuffer] -= x;
        double aux = (workVectors[FirstOrderDS::xBuffer]->norm2()) / (ds->normRef());
        if (aux > RelativeTol)
          _simulation->setRelativeConvergenceCriterionHeld(false);
      }
      DEBUG_PRINT("EulerMoreauOSI::updateState New value of x\n");
      DEBUG_EXPR(x.display());
    }
    else RuntimeException::selfThrow("EulerMoreauOSI::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}

void EulerMoreauOSI::display()
{
  OneStepIntegrator::display();

  std::cout << "====== EulerMoreauOSI OSI display ======" <<std::endl;

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;
    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    std::cout << "--------------------------------" <<std::endl;
    std::cout << "--> W of dynamical system number " << ds->number() << ": " <<std::endl;
    if (WMap[ds->number()]) WMap[ds->number()]->display();
    else std::cout << "-> NULL" <<std::endl;
    std::cout << "--> and corresponding theta is: " << _theta <<std::endl;
  }
  std::cout << "================================" <<std::endl;
}

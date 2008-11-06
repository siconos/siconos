/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */
#include "Moreau.h"
#include "MoreauXML.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "LagrangianLinearTIDS.h"
#include "FirstOrderLinearTIDS.h"

using namespace std;
using namespace DS;
using namespace RELATION;

// --- xml constructor ---
Moreau::Moreau(SP::OneStepIntegratorXML osiXML, SP::DynamicalSystemsSet dsList):
  OneStepIntegrator("Moreau")
{
  // Note: we do not call xml constructor of OSI, but default one, since we need to download theta and DS at the same time.

  if (!osiXML)
    RuntimeException::selfThrow("Moreau::xml constructor - OneStepIntegratorXML object == NULL.");

  integratorXml = osiXML;
  SP::MoreauXML moreauXml = boost::static_pointer_cast<MoreauXML>(osiXML);

  // Required inputs: a list of DS and one theta per DS.
  // No xml entries at the time for sizeMem and W.

  if (!osiXML->hasDSList())
    RuntimeException::selfThrow("Moreau::xml constructor - DS list is missing in xml input file.");
  if (!moreauXml->hasThetaList())
    RuntimeException::selfThrow("Moreau::xml constructor - theta list is missing in xml input file.");

  vector<double> thetaXml;     // list of theta values
  // thetaXml[i] will correspond to the ds number i in the xml list. If "all" attribute is true in ds,
  // then theta values are sorted so to correspond to growing ds numbers order.

  if (moreauXml->hasAllTheta()) // if one single value for all theta
    thetaXml.push_back(moreauXml->getSingleTheta());
  else
    moreauXml->getTheta(thetaXml);

  if (osiXML->hasAllDS()) // if flag all=true is present -> get all ds from the nsds
  {
    unsigned int i = 0;
    for (DSIterator it = dsList->begin(); it != dsList->end(); ++it)
    {
      OSIDynamicalSystems->insert(*it);
      // get corresponding theta. In xml they must be sorted in an order that corresponds to growing DS-numbers order.
      if (moreauXml->hasAllTheta()) // if one single value for all theta
        thetaMap[*it] = thetaXml[0];
      else
        thetaMap[*it] = thetaXml[i++];
    }
  }
  else
  {
    // get list of ds numbers implicated in the OSI
    vector<int> dsNumbers;
    osiXML->getDSNumbers(dsNumbers);
    unsigned int i = 0;
    // get corresponding ds and insert them into the set.
    vector<int>::iterator it;
    SP::DynamicalSystem ds;
    for_each(dsList->begin(), dsList->end(), boost::bind(&DynamicalSystem::getNumber, _1));
    for (it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
    {
      ds = dsList->getPtr(*it);
      OSIDynamicalSystems->insert(ds);
      if (moreauXml->hasAllTheta()) // if one single value for all theta
        thetaMap[ds] = thetaXml[0];
      else
        thetaMap[ds] = thetaXml[i++];
    }
  }
  // W loading: not yet implemented
  if (moreauXml->hasWList())
    RuntimeException::selfThrow("Moreau::xml constructor - W matrix loading not yet implemented.");
}

// --- constructor from a minimum set of data ---
Moreau::Moreau(SP::DynamicalSystem newDS, double newTheta): OneStepIntegrator("Moreau")
{
  OSIDynamicalSystems->insert(newDS);
  thetaMap[newDS] = newTheta;
}

// --- constructor from a set of data ---
Moreau::Moreau(DynamicalSystemsSet& newDS, double newTheta): OneStepIntegrator("Moreau", newDS)
{
  for (DSIterator itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
    thetaMap[*itDS] = newTheta;
}

Moreau::Moreau(DynamicalSystemsSet& newDS, const MapOfDouble& newTheta): OneStepIntegrator("Moreau", newDS)
{
  thetaMap = newTheta;
}

void Moreau::setWMap(const MapOfDSMatrices& newMap)
{
  // check sizes.
  if (newMap.size() != OSIDynamicalSystems->size())
    RuntimeException::selfThrow("Moreau::setWMap(newMap): number of W matrices is different from number of DS.");

  // pointer links! No reallocation
  ConstMatIterator it;

  WMap.clear();

  for (it = newMap.begin(); it != newMap.end(); ++it)
  {
    WMap[(*it).first] = (*it).second;
  }
}

const SimpleMatrix Moreau::getW(SP::DynamicalSystem ds)
{
  assert(ds &&
         "Moreau::getW(ds): ds == NULL.");
  //    return *(WMap[0]);
  assert(WMap[ds] &&
         "Moreau::getW(ds): W[ds] == NULL.");
  return *(WMap[ds]); // Copy !!
}

SP::SiconosMatrix Moreau::getWPtr(SP::DynamicalSystem ds)
{
  assert(ds && "Moreau::getWPtr(ds): ds == NULL.");
  //  return WMap[0];
  //  if(WMap[ds]==NULL)
  //    RuntimeException::selfThrow("Moreau::getWPtr(ds): W[ds] == NULL.");
  return WMap[ds];
}

void Moreau::setW(const SiconosMatrix& newValue, SP::DynamicalSystem ds)
{
  // Check if ds is in the OSI
  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("Moreau::setW(newVal,ds) - ds does not belong to this Integrator ...");

  // Check dimensions consistency
  unsigned int line = newValue.size(0);
  unsigned int col  = newValue.size(1);

  if (line != col) // Check that newValue is square
    RuntimeException::selfThrow("Moreau::setW(newVal,ds) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("Moreau::setW(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("Moreau::setW(newVal,ds) - unconsistent dimension between newVal and dynamical system to be integrated ");

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
      RuntimeException::selfThrow("Moreau - setW: inconsistent dimensions with problem size for given input matrix W");
  }
}

void Moreau::setWPtr(SP::SiconosMatrix newPtr, SP::DynamicalSystem ds)
{
  unsigned int line = newPtr->size(0);
  unsigned int col  = newPtr->size(1);
  if (line != col) // Check that newPtr is square
    RuntimeException::selfThrow("Moreau::setWPtr(newVal) - newVal is not square! ");

  if (!ds)
    RuntimeException::selfThrow("Moreau::setWPtr(newVal,ds) - ds == NULL.");

  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  if (line != sizeW) // check consistency between newValue and dynamical system size
    RuntimeException::selfThrow("Moreau::setW(newVal) - unconsistent dimension between newVal and dynamical system to be integrated ");

  WMap[ds] = newPtr;                  // link with new pointer
}

void Moreau::setThetaMap(const MapOfDouble& newMap) // useless function ?
{
  thetaMap = newMap;
}

const double Moreau::getTheta(SP::DynamicalSystem ds)
{
  if (!ds)
    RuntimeException::selfThrow("Moreau::getTheta(ds) - ds == NULL");
  return thetaMap[ds];
}

void Moreau::setTheta(double newTheta, SP::DynamicalSystem ds)
{
  if (!ds)
    RuntimeException::selfThrow("Moreau::setTheta(val,ds) - ds == NULL");
  thetaMap[ds] = newTheta;
}

void Moreau::initialize(SP::Simulation sim)
{
  OneStepIntegrator::initialize(sim);
  // Get initial time
  double t0 = simulationLink->getModelPtr()->getT0();
  // Compute W(t0) for all ds
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    // Memory allocation for workX. workX[ds*] corresponds to xfree (or vfree in lagrangian case).
    workX[*itDS].reset(new SimpleVector((*itDS)->getDim()));

    // W initialization
    initW(t0, *itDS);

    // Allocate memory for a work vector used in ds and in Newton process convergence evaluation.
    if ((*itDS)->getType() == LNLDS || (*itDS)->getType() == FONLDS)
      (*itDS)->allocateWorkVector(DynamicalSystem::NewtonSave, WMap[*itDS]->size(0));

    if ((*itDS)->getType() == FOLDS)
    {
      SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS>(*itDS);
      if (d->getBPtr())
        d->computeB(t0); // bi
    }
  }
}

void Moreau::initW(double t, SP::DynamicalSystem ds)
{
  // This function:
  // - allocate memory for a matrix W
  // - insert this matrix into WMap with ds as a key

  if (!ds)
    RuntimeException::selfThrow("Moreau::initW(t,ds) - ds == NULL");

  if (!OSIDynamicalSystems->isIn(ds))
    RuntimeException::selfThrow("Moreau::initW(t,ds) - ds does not belong to the OSI.");

  if (WMap.find(ds) != WMap.end())
    RuntimeException::selfThrow("Moreau::initW(t,ds) - W(ds) is already in the map and has been initialized.");

  // Memory allocation for W
  unsigned int sizeW = ds->getDim(); // n for first order systems, ndof for lagrangian.
  WMap[ds].reset(new SimpleMatrix(sizeW, sizeW));

  SP::SiconosMatrix W = WMap[ds];
  double h = simulationLink->getTimeStep();
  double theta = thetaMap[ds];
  DS::TYPES dsType = ds->getType();

  // 1 - First order non linear systems
  if (dsType == FONLDS)
  {
    // W =  M - h*theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Copy M or I if M is Null into W
    if (d->getMPtr())
      *W = *d->getMPtr();
    else
      W->eye();

    // d->computeJacobianXF(t); // Computation of JacXF is not required here
    // since it must have been done in OSI->initialize, before a call to this function.

    // Add -h*theta*jacobian_XF to W
    scal(-h * theta, *d->getJacobianXFPtr(), *W, false);
  }
  // 2 - First order linear systems
  else if (dsType == FOLDS || dsType == FOLTIDS)
  {
    SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS> (ds);
    if (d->getMPtr())
      *W = *d->getMPtr();
    else
      W->eye();

    scal(-h * theta, *d->getAPtr(), *W, false);
  }
  // 3 - Lagrangian non linear systems
  else if (dsType == LNLDS)
  {
    SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    SP::SiconosMatrix K = d->getJacobianFLPtr(0); // jacobian according to q
    SP::SiconosMatrix C = d->getJacobianFLPtr(1); // jacobian according to velocity

    *W = *d->getMassPtr();

    if (C)
      scal(-h * theta, *C, *W, false); // W -= h*theta*C

    if (K)
      scal(-h * h * theta * theta, *K, *W, false); //*W -= h*h*theta*theta**K;
  }
  // 4 - Lagrangian linear systems
  else if (dsType == LLTIDS)
  {
    SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);
    SP::SiconosMatrix K = d->getKPtr();
    SP::SiconosMatrix C = d->getCPtr();

    *W = *d->getMassPtr();

    if (C)
      scal(h * theta, *C, *W, false); // W += h*theta *C

    if (K)
      scal(h * h * theta * theta, *K, *W, false); // W = h*h*theta*theta*K
  }

  // === ===
  else RuntimeException::selfThrow("Moreau::computeW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized nor inversed here.
  // Function PLUForwardBackward will do that if required.
}

void Moreau::computeW(double t, SP::DynamicalSystem ds)
{
  // Compute W matrix of the Dynamical System ds, at time t and for the current ds state.

  // When this function is called, WMap[ds] is supposed to exist and not to be null
  // Memory allocation has been done during initW.

  assert(ds &&
         "Moreau::computeW(t,ds) - ds == NULL");

  assert((WMap.find(ds) != WMap.end()) &&
         "Moreau::computeW(t,ds) - W(ds) does not exists. Maybe you forget to initialize the osi?");

  double h = simulationLink->getTimeStep();
  double theta = thetaMap[ds];
  DS::TYPES dsType = ds->getType();

  SP::SiconosMatrix W = WMap[ds];

  // 1 - First order non linear systems
  if (dsType == FONLDS)
  {
    // W =  M - h*theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Copy M or I if M is Null into W
    if (d->getMPtr())
      *W = *d->getMPtr();
    else
      W->eye();

    d->computeJacobianXF(t);
    // Add -h*theta*jacobian_XF to W
    scal(-h * theta, *d->getJacobianXFPtr(), *W, false);
  }
  // 2 - First order linear systems
  else if (dsType == FOLDS || dsType == FOLTIDS)
  {
    SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS> (ds);
    if (dsType == FOLDS)
      d->computeA(t);

    if (d->getMPtr())
      *W = *d->getMPtr();
    else
      W->eye();
    scal(-h * theta, *d->getAPtr(), *W, false);
  }
  // 3 - Lagrangian non linear systems
  else if (dsType == LNLDS)
  {
    SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    SP::SiconosMatrix K = d->getJacobianFLPtr(0); // jacobian according to q
    SP::SiconosMatrix C = d->getJacobianFLPtr(1); // jacobian according to velocity

    d->computeMass();
    *W = *d->getMassPtr();

    if (C)
    {
      d->computeJacobianFL(1, t);
      scal(-h * theta, *C, *W, false); // W -= h*theta*C
    }

    if (K)
    {
      d->computeJacobianFL(0, t);
      scal(-h * h * theta * theta, *K, *W, false); //*W -= h*h*theta*theta**K;
    }
  }
  // 4 - Lagrangian linear systems
  else if (dsType == LLTIDS)
  {
    // Nothing: W does not depend on time.
  }

  // === ===
  else RuntimeException::selfThrow("Moreau::computeW - not yet implemented for Dynamical system type :" + dsType);

  // Remark: W is not LU-factorized here.
  // Function PLUForwardBackward will do that if required.
}

double Moreau::computeResidu()
{

  // This function is used to compute the residu for each "Moreau-discretized" dynamical system.
  // It then computes the norm of each of them and finally return the maximum
  // value for those norms.
  //
  // The state values used are those saved in the DS, ie the last computed ones.

  double t = simulationLink->getNextTime(); // End of the time step
  double told = simulationLink->getStartingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it;
  SP::DynamicalSystem ds; // Current Dynamical System.
  DS::TYPES dsType ; // Type of the current DS.
  double theta; // Theta parameter of the current ds.

  double maxResidu = 0;
  double normResidu = maxResidu;

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    theta = thetaMap[ds]; // Its theta parameter
    dsType = ds->getType(); // Its type
    // 1 - First Order Non Linear Systems
    if (dsType == FONLDS)
    {
      // Residu = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi) - h*r^k_i+1

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0); // xi

      // xFree pointer is used to compute and save Residu.
      SP::SiconosVector xfree = workX[d];

      SP::SiconosVector x = d->getXPtr(); // last saved value for x
      SP::SiconosMatrix M = d->getMPtr();
      if (M)
      {
        prod(*M, *xold, *xfree, true);
        *xfree *= -1.0;
        prod(*M, *x, *xfree, false); // xfree = M(x - xold)
      }
      else // if M == Identity
      {
        *xfree = *x;
        *xfree -= *xold;
      }

      if (d->getFPtr())  // if f exists
      {
        // computes f(ti,xi)
        d->computeF(told, xold);
        double coef = -h * (1 - theta);
        // xfree += coef * f_i
        scal(coef, *d->getFPtr(), *xfree, false);

        // computes f(ti+1, x_k,i+1) = f(t,x)
        d->computeF(t);
        coef = -h * theta;
        // xfree += coef * fL_k,i+1
        scal(coef, *d->getFPtr(), *xfree, false);
      }

      scal(-h, *d->getRPtr(), *xfree, false); // xfree = xfree - h*r
      normResidu = xfree->norm2();
    }

    // 2 - First Order Linear Systems
    else if (dsType == FOLDS)
    {
      // Residu = M(x_i+1 - x_i) - h*theta*(A_i+1*x_i+1 + b_i+1)  - h*(1-theta)*(Ai*xi + bi) - h*r_i+1

      SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS>(ds);

      SP::SiconosVector xfree = workX[d];

      SP::SiconosVector x = d->getXPtr(); // last saved value for x
      // x value at told
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0);

      SP::SiconosMatrix  M = d->getMPtr();
      // If M not equal to identity matrix
      if (M)
      {
        prod(*M, *xold, *xfree, true);
        *xfree *= -1.0;
        prod(*M, *x, *xfree, false); // xfree = M(x - xold)
      }
      else // if M == Identity
      {
        *xfree = *x;
        *xfree -= *xold;
      }

      SP::SiconosMatrix A = d->getAPtr();
      if (A)
      {
        // xFree += -h(1-theta)A_i*x_i
        d->computeA(told);
        double coeff = -h * (1 - theta);
        prod(coeff, *A, *xold, *xfree, false);
        // xFree += -h*theta*A_i+1*x_i+1
        d->computeA(t);
        coeff = -h * theta;
        prod(coeff, *A, *x, *xfree, false);

      }
      SP::SiconosVector b = d->getBPtr();
      if (b)
      {
        // xFree -= h(1-theta)*bi + h*theta*bi+1
        //        d->computeB(told); // bi
        double coeff = -h * (1 - theta);
        scal(coeff, *d->getBPtr(), *xfree, false);
        coeff = -h * theta;
        d->computeB(t); // bi+1
        scal(coeff, *d->getBPtr(), *xfree, false);
      }

      scal(-h, *d->getRPtr(), *xfree, false); // xfree = xfree - h*r
      normResidu = xfree->norm2();
    }
    // 2bis - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == FOLTIDS)
    {
      // Residu = W*(x_i+1 - xi) -hb -hA*xi- hr_i+1
      SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearTIDS>(ds);

      SP::SiconosVector xfree = workX[d];
      SP::SiconosVector x = d->getXPtr(); // last saved value for x
      // x value at told
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0);

      *xfree = *x;
      *xfree -= *xold;
      computeW(t, d); // W is recomputed since it has been LU-factorized during solving.
      SP::SiconosMatrix W = WMap[ds]; // Its W Moreau matrix of iteration.
      prod(*W, *xfree, *xfree, true);

      SP::SiconosMatrix A = d->getAPtr();
      if (A)
        prod(-h, *A, *xold, *xfree, false); // xfree -= h*A*xi

      SP::SiconosVector b = d->getBPtr();
      if (b)
        scal(-h, *b, *xfree, false); // xfree -= hb

      scal(-h, *d->getRPtr(), *xfree, false); // xfree = xfree - h*r
      normResidu = xfree->norm2();
    }
    // 3 - Lagrangian Non Linear Systems
    else if (dsType == LNLDS)
    {
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0);
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0);

      // vFree pointer is used to compute and save the residu.
      SP::SiconosVector vfree = workX[d];

      d->computeMass();
      SP::SiconosMatrix M = d->getMassPtr();
      SP::SiconosVector v = d->getVelocityPtr(); // v = v_k,i+1
      prod(*M, (*v - *vold), *vfree); // vfree = M(v - vold)

      if (d->getFLPtr())  // if fL exists
      {
        // computes fL(ti,vi,qi)
        d->computeFL(told, qold, vold);
        double coef = -h * (1 - theta);
        // vfree += coef * fL_i
        scal(coef, *d->getFLPtr(), *vfree, false);

        // computes fL(ti+1, v_k,i+1, q_k,i+1) = fL(t,v,q)
        d->computeFL(t);
        coef = -h * theta;
        // vfree += coef * fL_k,i+1
        scal(coef, *d->getFLPtr(), *vfree, false);
      }

      *vfree -= *d->getPPtr(2);
      normResidu = vfree->norm2();
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == LLTIDS)
    {
      // Residu = M(vi+1 - v_i) - h*theta*( Fext_i+1 - Kqi+1 - Cvi+1) - h*(1-theta)*(Fext_i -Cvi -kqi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---

      // vFree pointer is used to compute and save residu.
      SP::SiconosVector vfree = workX[d];

      SP::SiconosMatrix M = d->getMassPtr();
      SP::SiconosVector v = d->getVelocityPtr(); // v = v_k,i+1
      prod(*M, (*v - *vold), *vfree); // vfree = M(v - vold)

      double coeff;

      SP::SiconosMatrix C = d->getCPtr();
      if (C)
      {
        coeff = h * (1 - theta);
        prod(coeff, *C, *vold, *vfree, false); // vfree += h(1-theta)*C*vi
        coeff = h * theta;
        prod(coeff, *C, *v, *vfree, false); // vfree += h*theta*C*vi+1
      }
      SP::SiconosMatrix K = d->getKPtr();
      if (K)
      {
        SP::SiconosVector q = d->getQPtr();
        coeff = h * (1 - theta);
        prod(coeff, *K, *qold, *vfree, false); // vfree += h(1-theta)*K*qi
        coeff = h * theta;
        prod(coeff, *K, *q, *vfree, false); // vfree += h*theta*K*qi+1
      }

      SP::SiconosVector Fext = d->getFExtPtr();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = -h * (1 - theta);
        scal(coeff, *Fext, *vfree, false); // vfree -= h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = -h * theta;
        scal(coeff, *Fext, *vfree, false); // vfree -= h*theta * fext(ti+1)
      }

      *vfree -= *d->getPPtr(2);
      normResidu = vfree->norm2();
    }
    else
      RuntimeException::selfThrow("Moreau::computeResidu - not yet implemented for Dynamical system type: " + dsType);

    if (normResidu > maxResidu) maxResidu = normResidu;

  }
  return maxResidu;
}

void Moreau::computeFreeState()
{

  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  double t = simulationLink->getNextTime(); // End of the time step
  double told = simulationLink->getStartingTime(); // Beginning of the time step
  double h = t - told; // time step length
  //h=0.0100000;
  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SimpleVector *rold = static_cast<SimpleVector*>(d->getRMemoryPtr()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it; // Iterator through the set of DS.

  SP::DynamicalSystem ds; // Current Dynamical System.
  SP::SiconosMatrix W; // W Moreau matrix of the current DS.
  DS::TYPES dsType ; // Type of the current DS.
  double theta; // Theta parameter of the current ds.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    theta = thetaMap[ds]; // Its theta parameter
    dsType = ds->getType(); // Its type
    W = WMap[ds]; // Its W Moreau matrix of iteration.

    // 1 - First Order Non Linear Systems
    if (dsType == FONLDS)
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
      SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0); // xi

      // --- ResiduFree computation ---
      // ResiduFree = M(x-xold) - h*[theta*f(t) + (1-theta)*f(told)]
      //
      // xFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector xfree = workX[d];

      SP::SiconosVector x = d->getXPtr(); // last saved value for x

      SP::SiconosMatrix M = d->getMPtr();
      if (M)
        prod(*M, (*x - *xold), *xfree, true); // vfree = M(v - vold)
      else // if M == Identity
      {
        *xfree = *x;
        *xfree -= *xold;
      }

      if (d->getFPtr())  // if f exists
      {
        // computes f(ti,xi)
        d->computeF(told, xold);
        double coef = -h * (1 - theta);
        // xfree += coef * f_i
        scal(coef, *d->getFPtr(), *xfree, false);

        // computes f(ti+1, x_k,i+1) = f(t,x)
        d->computeF(t);
        coef = -h * theta;
        // vfree += coef * fL_k,i+1
        scal(coef, *d->getFPtr(), *xfree, false);
      }

      // -- xfree =  x - W^{-1} ResiduFree --
      // At this point xfree = residuFree
      // -> Solve WX = xfree and set xfree = X
      // -- Update W --
      computeW(t, d);

      W->PLUForwardBackwardInPlace(*xfree);
      // -> compute real xfree
      *xfree *= -1.0;
      *xfree += *x;
    }

    // 2 - First Order Linear Systems
    else if (dsType == FOLDS)
    {
      // xFree = [W]^{-1} * [  (M + h(1-theta)A_i )x_i + h*theta*b_i+1 + h*(1-theta)*b_i ]

      // IN to be updated at current time: W, b
      // IN at told: A, b, x
      // IN, not time dependant: M

      SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS>(ds);

      SP::SiconosVector xfree = workX[d];

      // x value at told
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0);

      // If M not equal to identity matrix
      SP::SiconosMatrix M = d->getMPtr();
      if (M)
        prod(*M, *xold, *xfree); // xFree = M*xi
      else
        *xfree = *xold;

      SP::SiconosMatrix A = d->getAPtr();
      if (A)
      {
        d->computeA(told);
        double coeff = h * (1 - theta);
        prod(coeff, *A, *xold, *xfree, false);
        // xFree += h(1-theta)A_i*x_i
      }
      SP::SiconosVector b = d->getBPtr();
      if (b)
      {
        // xFree += h(1-theta)*bi + h*theta*bi+1
        //        d->computeB(told); // bi
        scal(h * (1.0 - theta), *d->getBPtr(), *xfree, false);
        d->computeB(t); // bi+1
        scal(h * theta, *d->getBPtr(), *xfree, false);
      }

      // -- Update W --
      computeW(t, d);

      W->PLUForwardBackwardInPlace(*xfree); // Solve WX = xfree and set xfree = X.
    }
    // 2bis - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == FOLTIDS)
    {
      // xFree = [W]^{-1} * [  (W + hA )x_i + h*b ] = xi + hW^-1(Axi + b)

      // IN to be updated at current time: none
      // IN at told: x
      // IN, not time dependant: A,b,W

      SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearTIDS>(ds);

      SP::SiconosVector xfree = workX[d];
      // x value at told
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0);

      SP::SiconosMatrix A = d->getAPtr();
      if (A)
        prod(h, *A, *xold, *xfree, true); // xfree = h*A*xi
      else
        xfree->zero();

      SP::SiconosVector b = d->getBPtr();
      if (b)
        scal(h, *b, *xfree, false); // xfree += hb

      W->PLUForwardBackwardInPlace(*xfree); // Solve WX = xfree and set xfree = X.

      *xfree += *xold; // xfree = xi + ...
    }
    // 3 - Lagrangian Non Linear Systems
    else if (dsType == LNLDS)
    {
      // IN to be updated at current time: W, M, q, v, fL
      // IN at told: qi,vi, fLi

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      // vFree = v_k,i+1 - W^{-1} ResiduFree
      // with
      // ResiduFree = M(q_k,i+1)(v_k,i+1 - v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi)

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0);
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0);

      // --- ResiduFree computation ---
      // ResFree = M(v-vold) - h*[theta*fL(t) + (1-theta)*fL(told)]
      //
      // vFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector vfree = workX[d];

      // -- Update W --
      // Note: during computeW, mass and jacobians of fL will be computed/
      computeW(t, d);

      SP::SiconosMatrix M = d->getMassPtr();
      SP::SiconosVector v = d->getVelocityPtr(); // v = v_k,i+1
      prod(*M, (*v - *vold), *vfree); // vfree = M(v - vold)

      if (d->getFLPtr())  // if fL exists
      {
        // computes fL(ti,vi,qi)
        d->computeFL(told, qold, vold);
        double coef = -h * (1 - theta);
        // vfree += coef * fL_i
        scal(coef, *d->getFLPtr(), *vfree, false);

        // computes fL(ti+1, v_k,i+1, q_k,i+1) = fL(t,v,q)
        d->computeFL(t);
        coef = -h * theta;
        // vfree += coef * fL_k,i+1
        scal(coef, *d->getFLPtr(), *vfree, false);
      }

      // -- vfree =  v - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      W->PLUForwardBackwardInPlace(*vfree);
      // -> compute real vfree
      *vfree *= -1.0;
      *vfree += *v;

      // -- Computes *qfree = (*qold) + h * (theta * (*vfree) + (1.0 - theta) * (*vold))  --
      //    // May be useless?
      //    SiconosVector *qfree = d->getQFreePtr();
      //    double coeff = h*(1-theta);
      //    scal(coeff,*vold,*qfree); // qfree = coeff*vold
      //    coeff = h*theta;
      //    scal(coeff,*vfree,*qfree,false); // qfree += coeff*vfree
      //    *qfree  += *qold;
      //  cout << "vfree :" << endl;
      //    vfree->display();
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == LLTIDS)
    {
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // vFree = v_i + W^{-1} ResiduFree
      // with
      // ResiduFree = (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 + h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---

      // vFree pointer is used to compute and save ResiduFree in this first step.

      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SP::SiconosVector vfree = workX[d];
      vfree->zero();
      double coeff;
      // -- No need to update W --
      SP::SiconosMatrix C = d->getCPtr();
      if (C)
        prod(-h, *C, *vold, *vfree, false); // vfree += -h*C*vi

      SP::SiconosMatrix K = d->getKPtr();
      if (K)
      {
        coeff = -h * h * theta;
        prod(coeff, *K, *vold, *vfree, false); // vfree += -h^2*theta*K*vi
        prod(-h, *K, *qold, *vfree, false); // vfree += -h*K*qi
      }

      SP::SiconosVector Fext = d->getFExtPtr();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - theta);
        scal(coeff, *Fext, *vfree, false); // vfree += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * theta;
        scal(coeff, *Fext, *vfree, false); // vfree += h*theta * fext(ti+1)
      }

      // -- vfree =  vi + W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      W->PLUForwardBackwardInPlace(*vfree);
      *vfree += *vold;
      // -- Computes *qfree = (*qold) + h * (theta * (*vfree) + (1.0 - theta) * (*vold))  --
      // May be useless?
      //    SiconosVector *qfree = d->getQFreePtr();
      //    coeff = h*(1-theta);
      //    scal(coeff, *vold,*qfree); // qfree = coeff*vold
      //    coeff = h*theta;
      //    scal(coeff,*vfree,*qfree,false); // qfree += coeff*vfree
      //    *qfree  += *qold;
    }
    else
      RuntimeException::selfThrow("Moreau::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
  }
}


void Moreau::integrate(double& tinit, double& tend, double& tout, int&)
{
  // Last parameter is not used (required for Lsodar but not for Moreau).

  double h = tend - tinit;
  tout = tend;

  DSIterator it;
  SP::SiconosMatrix W;
  double theta;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds];
    theta = thetaMap[ds];
    DS::TYPES dsType = ds->getType();

    if (dsType == LLTIDS)
    {
      // get the ds
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);
      // get velocity pointers for current time step
      SP::SiconosVector v = d->getVelocityPtr();
      // get q and velocity pointers for previous time step
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0);
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0);
      // get p pointer
      SP::SiconosVector p = d->getPPtr(2);

      // velocity computation :
      //
      // v = vi + W^{-1}[ -h*C*vi - h*h*theta*K*vi - h*K*qi + h*theta*Fext(t) + h*(1-theta) * Fext(ti) ] + W^{-1}*pi+1
      //

      *v = *p; // v = p

      double coeff;
      // -- No need to update W --
      SP::SiconosMatrix C = d->getCPtr();
      if (C)
        prod(-h, *C, *vold, *v, false); // v += -h*C*vi

      SP::SiconosMatrix K = d->getKPtr();
      if (K)
      {
        coeff = -h * h * theta;
        prod(coeff, *K, *vold, *v, false); // v += -h^2*theta*K*vi
        prod(-h, *K, *qold, *v, false); // v += -h*K*qi
      }

      SP::SiconosVector Fext = d->getFExtPtr();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(tinit);
        coeff = h * (1 - theta);
        scal(coeff, *Fext, *v, false); // v += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(tout);
        coeff = h * theta;
        scal(coeff, *Fext, *v, false); // v += h*theta * fext(ti+1)
      }
      // -> Solve WX = v and set v = X
      W->PLUForwardBackwardInPlace(*v);
      *v += *vold;
    }
    else RuntimeException::selfThrow("Moreau::integrate - not yet implemented for Dynamical system type :" + dsType);
  }
}

void Moreau::updateState(unsigned int level)
{
  double h = simulationLink->getTimeStep();

  DSIterator it;
  SP::SiconosMatrix W;
  double theta;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds];
    theta = thetaMap[ds];
    // Get the DS type

    DS::TYPES dsType = ds->getType();

    // 1 - First Order Systems
    if (dsType == FONLDS || dsType == FOLDS || dsType == FOLTIDS)
    {
      SP::FirstOrderNonLinearDS fonlds = boost::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      SP::SiconosVector x = ds->getXPtr();
      //    SP::SiconosVector xFree = fonlds->getXFreePtr();

      //          integration of r with theta method removed
      //      *x = h * theta * *fonlds->getRPtr();
      // Save value of q in stateTmp for future convergence computation
      if (dsType == FONLDS)
        ds->addWorkVector(x, DynamicalSystem::NewtonSave);

      // Solve W(x-xfree) = hr
      scal(h, *fonlds->getRPtr(), *x); // x = h*r
      W->PLUForwardBackwardInPlace(*x); // x =h* W^{-1} *r
      *x += *workX[ds]; // x+=xfree
    }
    // 3 - Lagrangian Systems
    else if (dsType == LNLDS || dsType == LLTIDS)
    {
      // get dynamical system
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      //    SiconosVector *vfree = d->getVelocityFreePtr();
      SP::SiconosVector v = d->getVelocityPtr();

      // To compute v, we solve W(v - vfree) = p
      *v = *d->getPPtr(level); // v = p
      W->PLUForwardBackwardInPlace(*v);
      *v += *workX[ds];

      // Compute q
      SP::SiconosVector q = d->getQPtr();
      // Save value of q in stateTmp for future convergence computation
      if (dsType == LNLDS)
        ds->addWorkVector(q, DynamicalSystem::NewtonSave);

      //  -> get previous time step state
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0);
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0);
      // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
      double coeff = h * theta;
      scal(coeff, *v, *q) ; // q = h*theta*v
      coeff = h * (1 - theta);
      scal(coeff, *vold, *q, false); // q += h(1-theta)*vold
      *q += *qold;
    }
    else RuntimeException::selfThrow("Moreau::updateState - not yet implemented for Dynamical system type: " + dsType);
  }
}


void Moreau::display()
{
  OneStepIntegrator::display();

  cout << "====== Moreau OSI display ======" << endl;
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    cout << "--------------------------------" << endl;
    cout << "--> W of dynamical system number " << (*it)->getNumber() << ": " << endl;
    if (WMap[*it]) WMap[*it]->display();
    else cout << "-> NULL" << endl;
    cout << "--> and corresponding theta is: " << thetaMap[*it] << endl;
  }
  cout << "================================" << endl;
}

void Moreau::saveWToXML()
{
  //   if(integratorXml != NULL)
  //     {
  //       (static_cast<MoreauXML*>(integratorXml))->setW(W);
  //     }
  //   else RuntimeException::selfThrow("Moreau::saveIntegratorToXML - IntegratorXML object not exists");
  RuntimeException::selfThrow("Moreau::saveWToXML -  not yet implemented.");
}

Moreau* Moreau::convert(OneStepIntegrator* osi)
{
  Moreau* moreau = dynamic_cast<Moreau*>(osi);
  return moreau;
}

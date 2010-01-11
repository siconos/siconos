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
#include "Moreau.hpp"
#include "MoreauXML.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "NewtonEulerDS.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "FirstOrderLinearTIDS.hpp"

using namespace std;
using namespace DS;
using namespace RELATION;

// --- xml constructor ---
Moreau::Moreau(SP::OneStepIntegratorXML osiXML, SP::DynamicalSystemsSet dsList):
  OneStepIntegrator(OSI::MOREAU)
{
  // Note: we do not call xml constructor of OSI, but default one, since we need to download _theta and DS at the same time.

  if (!osiXML)
    RuntimeException::selfThrow("Moreau::xml constructor - OneStepIntegratorXML object == NULL.");

  integratorXml = osiXML;
  SP::MoreauXML moreauXml = boost::static_pointer_cast<MoreauXML>(osiXML);

  // Required inputs: a list of DS and one _theta per DS.
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
        _theta = thetaXml[0];
      else
      {
        // should not happen
        RuntimeException::selfThrow("Multiples theta values for Moreau integrator are not valid anymore : use several Moreau instantiation instead.");
      }
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
    for_each(dsList->begin(), dsList->end(), boost::bind(&DynamicalSystem::number, _1));
    for (it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
    {
      ds = dsList->getPtr(*it);
      OSIDynamicalSystems->insert(ds);
      if (moreauXml->hasAllTheta()) // if one single value for all theta
        _theta = thetaXml[0];
      else
      {
        RuntimeException::selfThrow("Multiples theta values for Moreau integrator are not valid anymore : use several Moreau instantiation instead.");
      }
    }
  }
  // W loading: not yet implemented
  if (moreauXml->hasWList())
    RuntimeException::selfThrow("Moreau::xml constructor - W matrix loading not yet implemented.");
}

// --- constructor from a minimum set of data ---
Moreau::Moreau(SP::DynamicalSystem newDS, double newTheta) :
  OneStepIntegrator(OSI::MOREAU)
{
  OSIDynamicalSystems->insert(newDS);
  _theta = newTheta;
}

// --- constructor from a set of data ---
Moreau::Moreau(DynamicalSystemsSet& newDS, double newTheta):
  OneStepIntegrator(OSI::MOREAU, newDS)
{
  _theta = newTheta;
}

// Note: OSIDynamicalSystems and thetaMap must disappear
void Moreau::insertDynamicalSystem(SP::DynamicalSystem ds)
{
  OSIDynamicalSystems->insert(ds);
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

SP::SiconosMatrix Moreau::W(SP::DynamicalSystem ds)
{
  assert(ds && "Moreau::W(ds): ds == NULL.");
  //  return WMap[0];
  //  if(WMap[ds]==NULL)
  //    RuntimeException::selfThrow("Moreau::W(ds): W[ds] == NULL.");
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

inline const double Moreau::theta()
{
  return _theta;
}

void Moreau::setTheta(double newTheta)
{
  _theta = newTheta;
}

void Moreau::initialize(SP::Simulation sim)
{
  OneStepIntegrator::initialize(sim);
  // Get initial time
  double t0 = simulationLink->model()->t0();
  // Compute W(t0) for all ds
  ConstDSIterator itDS;
  for (itDS = OSIDynamicalSystems->begin(); itDS != OSIDynamicalSystems->end(); ++itDS)
  {
    // Memory allocation for workX. workX[ds*] corresponds to xfree (or vfree in lagrangian case).
    // workX[*itDS].reset(new SimpleVector((*itDS)->getDim()));

    // W initialization
    initW(t0, *itDS);

    //      if ((*itDS)->getType() == LNLDS || (*itDS)->getType() == FONLDS)
    (*itDS)->allocateWorkVector(DynamicalSystem::local_buffer, WMap[*itDS]->size(0));

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
  double h = simulationLink->timeStep();
  DS::TYPES dsType = ds->getType();

  // 1 - First order non linear systems
  if (dsType == FONLDS || dsType == FOLDS || dsType == FOLTIDS)
  {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Copy M or I if M is Null into W
    if (d->M())
      *W = *d->M();
    else
      W->eye();

    // d->computeJacobianXF(t); // Computation of JacxF is not required here
    // since it must have been done in OSI->initialize, before a call to this function.

    // Add -h*_theta*jacobian_XF to W
    scal(-h * _theta, *d->jacobianXF(), *W, false);
  }
  // 2 - First order linear systems
  //   else if (dsType == FOLDS || dsType == FOLTIDS)
  //     {
  //       SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS> (ds);
  //       if( d->M() )
  //  *W = *d->M();
  //       else
  //  W->eye();

  //       scal(-h*_theta, *d->A(),*W,false);
  //     }
  // 3 - Lagrangian non linear systems
  else if (dsType == LNLDS)
  {
    SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    SP::SiconosMatrix K = d->jacobianqFL(); // jacobian according to q
    SP::SiconosMatrix C = d->jacobianqDotFL(); // jacobian according to velocity

    *W = *d->mass();

    if (C)
      scal(-h * _theta, *C, *W, false); // W -= h*_theta*C

    if (K)
      scal(-h * h * _theta * _theta, *K, *W, false); //*W -= h*h*_theta*_theta**K;
  }
  // 4 - Lagrangian linear systems
  else if (dsType == LLTIDS)
  {
    SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);
    SP::SiconosMatrix K = d->K();
    SP::SiconosMatrix C = d->C();

    *W = *d->mass();

    if (C)
      scal(h * _theta, *C, *W, false); // W += h*_theta *C

    if (K)
      scal(h * h * _theta * _theta, *K, *W, false); // W = h*h*_theta*_theta*K
  }

  // === ===
  else if (dsType == NENLDS)
  {
    SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

    //    scal(1.0/h,*(d->luW()),*(d->luW()));
    cout << "Moreau::initW luW before LUFact\n";
    d->luW()->display();
    d->luW()->PLUFactorizationInPlace();

  }
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

  double h = simulationLink->timeStep();
  DS::TYPES dsType = ds->getType();

  SP::SiconosMatrix W = WMap[ds];

  // 1 - First order non linear systems
  if (dsType == FONLDS)
  {
    // W =  M - h*_theta* [jacobian_x f(t,x,z)]
    SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS> (ds);

    // Copy M or I if M is Null into W
    if (d->M())
      *W = *d->M();
    else
      W->eye();

    d->computeJacobianXF(t);
    // Add -h*_theta*jacobian_XF to W
    scal(-h * _theta, *d->jacobianXF(), *W, false);
  }
  // 2 - First order linear systems
  else if (dsType == FOLDS || dsType == FOLTIDS)
  {
    SP::FirstOrderLinearDS d = boost::static_pointer_cast<FirstOrderLinearDS> (ds);
    if (dsType == FOLDS)
      d->computeA(t);

    if (d->M())
      *W = *d->M();
    else
      W->eye();
    scal(-h * _theta, *d->A(), *W, false);
  }
  // 3 - Lagrangian non linear systems
  else if (dsType == LNLDS)
  {
    SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);
    SP::SiconosMatrix K = d->jacobianqFL(); // jacobian according to q
    SP::SiconosMatrix C = d->jacobianqDotFL(); // jacobian according to velocity

    d->computeMass();
    *W = *d->mass();

    if (C)
    {
      d->computeJacobianqDotFL(t);
      scal(-h * _theta, *C, *W, false); // W -= h*_theta*C
    }

    if (K)
    {
      d->computeJacobianqFL(t);
      scal(-h * h * _theta * _theta, *K, *W, false); //*W -= h*h*_theta*_theta**K;
    }
  }
  // 4 - Lagrangian linear systems
  else if (dsType == LLTIDS)
  {
    // Nothing: W does not depend on time.
  }

  // === ===
  else if (dsType == NENLDS)
    ;
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
  DS::TYPES dsType ; // Type of the current DS.

  double maxResidu = 0;
  double normResidu = maxResidu;

  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    dsType = ds->getType(); // Its type
    SP::SiconosVector residuFree = ds->residuFree();
    // 1 - First Order Non Linear Systems
    if (dsType == FONLDS || dsType == FOLDS)
    {
      // ResiduFree = M(x_k,i+1 - x_i) - h*theta*f(t,x_k,i+1) - h*(1-theta)*f(ti,xi)
      // Residu = Residu - h*r^k_i+1
      //  $\mathcal R(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) - h r$
      //  $\mathcal R_{free}(x,r) = M(x - x_{k}) -h\theta f( x , t_{k+1}) - h(1-\theta)f(x_k,t_k) $

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.

      SP::FirstOrderNonLinearDS d = boost::static_pointer_cast<FirstOrderNonLinearDS>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

      SP::SiconosVector x = d->x(); // last saved value for x
      SP::SiconosMatrix M = d->M();

      *residuFree = *x;
      *residuFree -= *xold;
      if (M)
        prod(*M, *residuFree, *residuFree, true);

      if (d->f())
      {
        // computes f(ti,xi)
        //        d->computeF(told,xold);
        double coef = -h * (1 - _theta);
        // residuFree += coef * f_i
        scal(coef, *d->fold(), *residuFree, false);

        // computes f(ti+1, x_k,i+1) = f(t,x)
        d->computeF(t);
        coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        scal(coef, *d->f(), *residuFree, false);
      }
      //      cout<<"Moreau: residu free"<<endl;
      //      (*residuFree).display();
      (*(d->workFree())) = *residuFree;
      scal(-h, *d->r(), (*d->workFree()), false); // residu = residu - h*r
      normResidu = d->workFree()->norm2();
      //    cout<<"Moreau: residu "<<endl;
      //    (workX[d])->display();
      //    cout<<"Moreau: norm residu :"<<normResidu<<endl;


      (*d->residur()) = (*d->r()) - (*d->gAlpha());

      //      cout<<"Moreau FONLDS: residu r"<<endl;
      //      (*d->residur()).display();
    }
    // 2 - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == FOLTIDS)
    {
      SP::FirstOrderLinearTIDS d = boost::static_pointer_cast<FirstOrderLinearTIDS>(ds);
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
      SP::SiconosVector xBuffer = d->getWorkVector(DynamicalSystem::local_buffer);
      *xBuffer = *(d->xMemory()->getSiconosVector(0));
      //    cout<<"Moreau TIDS::computeResidu: x_k"<<endl;
      //    xBuffer->display();

      SP::SiconosMatrix A = d->A();
      if (A)
        prod(*A, *xBuffer, *residuFree, false); // residuFree -= -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b_{k+1}


      *residuFree *= -h;


      //    cout<<"Moreau TIDS::computeResidu: residu free"<<endl;
      //    (*residuFree).display();


      /*SP::SiconosVector xBuffer = d->getWorkVector(DynamicalSystem::local_buffer);
      *xBuffer = theta * (*x_alpha);
      *xBuffer += (1-theta)*(*x_k);
      *xBuffer *=-h;
      SP::SiconosMatrix A = d->A();
      if( A )
        prod(*A,*xBuffer,*residuFree, false);  // residuFree -= -h( A (\theta x^{\alpha}_{k+1} + (1-\theta)  x_k) +b_{k+1}

      *xBuffer =  (* x_alpha);
      *xBuffer -= (*x_k);

       SP::SiconosMatrix M = d->M();
       if(M)
         prod(*M,*xBuffer,*residuFree, false);// residuFree += M(x^{\alpha}_{k+1} - x_{k})
       else
         *residuFree+=*xBuffer;

      (*workX[d])=*residuFree;
      scal(-h, *d->r(), (*workX[d]), false); // residu = residu - h*r
      normResidu = (workX[d])->norm2();*/
    }
    // 3 - Lagrangian Non Linear Systems
    else if (dsType == LNLDS)
    {
      // residu = M(q*)(v_k,i+1 - v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      SP::SiconosVector q = d->q();


      d->computeMass();
      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      prod(*M, (*v - *vold), *residuFree); // residuFree = M(v - vold)
      if (d->fL())  // if fL exists
      {
        // computes fL(ti,vi,qi)
        d->computeFL(told, qold, vold);
        double coef = -h * (1 - _theta);
        // residuFree += coef * fL_i
        scal(coef, *d->fL(), *residuFree, false);

        // computes fL(ti+1, v_k,i+1, q_k,i+1) = fL(t,v,q)

        // d->computeFL(t,q,v) ?
        d->computeFL(t);
        coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        scal(coef, *d->fL(), *residuFree, false);
      }
      *(d->workFree()) = *residuFree;
      *(d->workFree()) -= *d->p(2);
      normResidu = d->workFree()->norm2();
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == LLTIDS)
    {
      // ResiduFree = M(vi+1 - v_i) - h*theta*( Fext_i+1 - Kqi+1 - Cvi+1) - h*(1-theta)*(Fext_i -Cvi -kqi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---
      residuFree->zero();
      double coeff;
      // -- No need to update W --
      SP::SiconosMatrix C = d->C();
      if (C)
        prod(-h, *C, *vold, *residuFree, false); // vfree += -h*C*vi

      SP::SiconosMatrix K = d->K();
      if (K)
      {
        coeff = -h * h * _theta;
        prod(coeff, *K, *vold, *residuFree, false); // vfree += -h^2*_theta*K*vi
        prod(-h, *K, *qold, *residuFree, false); // vfree += -h*K*qi
      }

      SP::SiconosVector Fext = d->fExt();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - _theta);
        scal(coeff, *Fext, *residuFree, false); // vfree += h*(1-_theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * _theta;
        scal(coeff, *Fext, *residuFree, false); // vfree += h*_theta * fext(ti+1)
      }

      (* d->workFree()) = *residuFree;
      (* d->workFree()) -= *d->p(2);
      normResidu = d->workFree()->norm2();

    }
    else if (dsType == NENLDS)
    {
      // residu = M(q*)(v_k,i+1 - v_i) - h*_theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-_theta)*fL(ti,vi,qi) - pi+1

      // -- Convert the DS into a Lagrangian one.
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      SP::SiconosVector q = d->q();


      SP::SiconosMatrix W = d->W();
      SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      prod(*W, (*v - *vold), *residuFree); // residuFree = M(v - vold)
      if (d->fL())  // if fL exists
      {
        // computes fL(ti,vi,qi)
        SP::SiconosVector fLold = d->fLMemory()->getSiconosVector(0);
        double coef = -h * (1 - _theta);
        // residuFree += coef * fL_i
        scal(coef, *fLold, *residuFree, false);

        // computes fL(ti+1, v_k,i+1, q_k,i+1) = fL(t,v,q)

        // d->computeFL(t,q,v) ?
        d->computeFL(t);
        coef = -h * _theta;
        // residuFree += coef * fL_k,i+1
        scal(coef, *d->fL(), *residuFree, false);
      }
      *(d->workFree()) = *residuFree;
      //    cout<<"Moreau::computeResidu :\n";
      //    residuFree->display();
      *(d->workFree()) -= *d->p(2);
      normResidu = d->workFree()->norm2();
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

  double t = simulationLink->nextTime(); // End of the time step
  double told = simulationLink->startingTime(); // Beginning of the time step
  double h = t - told; // time step length
  //h=0.0100000;
  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SimpleVector *rold = static_cast<SimpleVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it; // Iterator through the set of DS.

  SP::DynamicalSystem ds; // Current Dynamical System.
  SP::SiconosMatrix W; // W Moreau matrix of the current DS.
  DS::TYPES dsType ; // Type of the current DS.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    dsType = ds->getType(); // Its type
    W = WMap[ds]; // Its W Moreau matrix of iteration.

    // 1 - First Order Non Linear Systems
    if (dsType == FONLDS || dsType == FOLDS || dsType == FOLTIDS)
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
      //    SP::SiconosVector xold = d->xMemory()->getSiconosVector(0); // xi

      // --- ResiduFree computation ---
      // ResiduFree = M(x-xold) - h*[theta*f(t) + (1-theta)*f(told)]
      //
      // xFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector xfree = d->workFree();//workX[d];
      *xfree = *(d->residuFree());


      SP::SiconosVector x = d->x(); // last saved value for x


      // -- xfree =  x - W^{-1} ResiduFree --
      // At this point xfree = residuFree
      // -> Solve WX = xfree and set xfree = X
      // -- Update W --
      if (dsType != FOLTIDS)
        computeW(t, d);

      W->PLUForwardBackwardInPlace(*xfree);


      // -> compute real xfree
      *xfree *= -1.0;
      *xfree += *x;
      //    cout<<" moreau::computefreestate xfree"<<endl;
      //    xfree->display();

      if (!simulationLink->model()->nonSmoothDynamicalSystem()->isLinear())
      {
        SP::SiconosVector xp = d->xp();
        //      cout<<"before moreau::computefreestate xp"<<endl;
        //      xp->display();
        W->PLUForwardBackwardInPlace(*xp);
        scal(h, *xp, *xp);
        *xp += *xfree;
        //      cout<<"after moreau::computefreestate xp"<<endl;
        //      xp->display();
        SP::SiconosVector xq = d->xq();
        *xq = *xp;
        *xq -= *x;
      }

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
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      // --- ResiduFree computation ---
      // ResFree = M(v-vold) - h*[theta*fL(t) + (1-theta)*fL(told)]
      //
      // vFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector vfree = d->workFree();//workX[d];
      (*vfree) = *(d->residuFree());

      // -- Update W --
      // Note: during computeW, mass and jacobians of fL will be computed/
      computeW(t, d);
      SP::SiconosVector v = d->velocity(); // v = v_k,i+1



      // -- vfree =  v - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      W->PLUForwardBackwardInPlace(*vfree);
      // -> compute real vfree
      *vfree *= -1.0;
      *vfree += *v;

    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == LLTIDS)
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
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---

      // vFree pointer is used to compute and save ResiduFree in this first step.

      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SP::SiconosVector vfree = d->workFree();//workX[d];
      (*vfree) = *(d->residuFree());


      W->PLUForwardBackwardInPlace(*vfree);

      *vfree += *vold;
    }
    else if (dsType == NENLDS)
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
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);

      // --- ResiduFree computation ---
      // ResFree = M(v-vold) - h*[theta*fL(t) + (1-theta)*fL(told)]
      //
      // vFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector vfree = d->workFree();//workX[d];
      (*vfree) = *(d->residuFree());

      // -- Update W --
      // Note: during computeW, mass and jacobians of fL will be computed/
      SP::SiconosVector v = d->velocity(); // v = v_k,i+1



      // -- vfree =  v - W^{-1} ResiduFree --
      // At this point vfree = residuFree
      // -> Solve WX = vfree and set vfree = X
      //    cout<<"Moreau::computeFreeState residu free"<<endl;
      //    vfree->display();
      d->luW()->PLUForwardBackwardInPlace(*vfree);
      //    cout<<"Moreau::computeFreeState -WRfree"<<endl;
      //    vfree->display();
      //    scal(h,*vfree,*vfree);
      // -> compute real vfree
      *vfree *= -1.0;
      *vfree += *v;

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
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds];
    DS::TYPES dsType = ds->getType();

    if (dsType == LLTIDS)
    {
      // get the ds
      SP::LagrangianLinearTIDS d = boost::static_pointer_cast<LagrangianLinearTIDS> (ds);
      // get velocity pointers for current time step
      SP::SiconosVector v = d->velocity();
      // get q and velocity pointers for previous time step
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      // get p pointer
      SP::SiconosVector p = d->p(2);

      // velocity computation :
      //
      // v = vi + W^{-1}[ -h*C*vi - h*h*theta*K*vi - h*K*qi + h*theta*Fext(t) + h*(1-theta) * Fext(ti) ] + W^{-1}*pi+1
      //

      *v = *p; // v = p

      double coeff;
      // -- No need to update W --
      SP::SiconosMatrix C = d->C();
      if (C)
        prod(-h, *C, *vold, *v, false); // v += -h*C*vi

      SP::SiconosMatrix K = d->K();
      if (K)
      {
        coeff = -h * h * _theta;
        prod(coeff, *K, *vold, *v, false); // v += -h^2*theta*K*vi
        prod(-h, *K, *qold, *v, false); // v += -h*K*qi
      }

      SP::SiconosVector Fext = d->fExt();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(tinit);
        coeff = h * (1 - _theta);
        scal(coeff, *Fext, *v, false); // v += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(tout);
        coeff = h * _theta;
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
  double h = simulationLink->timeStep();

  const double& RelativeTol = simulationLink->getRelativeConvergenceTol();
  bool useRCC = simulationLink->getUseRelativeConvergenceCriteron();
  if (useRCC)
    simulationLink->setRelativeConvergenceCriterionHeld(true);

  DSIterator it;
  SP::SiconosMatrix W;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    SP::DynamicalSystem ds = *it;
    W = WMap[ds];
    // Get the DS type

    DS::TYPES dsType = ds->getType();

    // 1 - First Order Systems
    if (dsType == FONLDS || dsType == FOLDS || dsType == FOLTIDS)
    {
      SP::FirstOrderNonLinearDS fonlds = boost::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      SP::SiconosVector x = ds->x();
      bool baux = (useRCC && dsType == FONLDS && simulationLink->getRelativeConvergenceCriterionHeld());
      //    SP::SiconosVector xFree = fonlds->xFree();

      // Save value of q in local_buffer for relative convergence computation
      if (baux)
        ds->addWorkVector(x, DynamicalSystem::local_buffer);


      // Solve W(x-xfree) = hr
      scal(h, *fonlds->r(), *x); // x = h*r
      W->PLUForwardBackwardInPlace(*x); // x =h* W^{-1} *r
      *x += *(fonlds->workFree()); //*workX[ds]; // x+=xfree
      if (baux)
      {
        ds->subWorkVector(x, DynamicalSystem::local_buffer);
        double aux = ((ds->getWorkVector(DynamicalSystem::local_buffer))->norm2()) / (ds->normRef());
        if (aux > RelativeTol)
          simulationLink->setRelativeConvergenceCriterionHeld(false);
      }

      //  }else if (dsType == FOLTIDS){
      //    SP::FirstOrderNonLinearDS fonlds = boost::static_pointer_cast<FirstOrderNonLinearDS>(ds);
      //    SP::SiconosVector x = ds->x();


      //    // Solve W(x-xfree) = hr
      //    *x=*fonlds->r();
      //    W->PLUForwardBackwardInPlace(*x); // x = W^{-1} *r
      //    scal(h,*x,*x); // x = h*W^{-1}*r
      //    *x +=*(fonlds->xfree());//*workX[ds]; // x+=xfree
      //    //    cout<<"X alpha+1"<<endl;
      //    //    x->display();
    }
    // 3 - Lagrangian Systems
    else if (dsType == LNLDS || dsType == LLTIDS)
    {
      // get dynamical system
      SP::LagrangianDS d = boost::static_pointer_cast<LagrangianDS> (ds);

      //    SiconosVector *vfree = d->velocityFree();
      SP::SiconosVector v = d->velocity();
      bool baux = dsType == LNLDS && useRCC && simulationLink->getRelativeConvergenceCriterionHeld();

      // To compute v, we solve W(v - vfree) = p
      *v = *d->p(level); // v = p
      W->PLUForwardBackwardInPlace(*v);
      *v +=  * ds->workFree();

      // Compute q
      SP::SiconosVector q = d->q();
      // Save value of q in stateTmp for future convergence computation
      if (baux)
        ds->addWorkVector(q, DynamicalSystem::local_buffer);

      //  -> get previous time step state
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0);
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
      double coeff = h * _theta;
      scal(coeff, *v, *q) ; // q = h*theta*v
      coeff = h * (1 - _theta);
      scal(coeff, *vold, *q, false); // q += h(1-theta)*vold
      *q += *qold;

      if (baux)
      {
        ds->subWorkVector(q, DynamicalSystem::local_buffer);
        double aux = ((ds->getWorkVector(DynamicalSystem::local_buffer))->norm2()) / (ds->normRef());
        if (aux > RelativeTol)
          simulationLink->setRelativeConvergenceCriterionHeld(false);
      }

    }
    else if (dsType == NENLDS)
    {
      // get dynamical system
      SP::NewtonEulerDS d = boost::static_pointer_cast<NewtonEulerDS> (ds);
      SP::SiconosVector v = d->velocity();

      //  cout<<"Moreau::updatestate prev v"<<endl;
      //  v->display();

      /*d->p has been fill by the Relation->computeInput, it constains B \lambda _{k+1}*/
      *v = *d->p(level); // v = p
      d->luW()->PLUForwardBackwardInPlace(*v);

      //  cout<<"Moreau::updatestate hWB lambda"<<endl;
      //  v->display();

      *v +=  * ds->workFree();

      //  cout<<"Moreau::updatestate work free"<<endl;
      //  ds->workFree()->display();
      //  cout<<"Moreau::updatestate new v"<<endl;
      //  v->display();

      //compute q
      //first step consists in computing  \dot q.
      //second step consists in updating q.
      //
      SP::SiconosMatrix T = d->T();
      SP::SiconosVector dotq = d->dotq();
      prod(*T, *v, *dotq, true);
      //  cout<<"Moreau::updateState v"<<endl;
      //  v->display();
      //  cout<<"Moreau::updateState dotq"<<endl;
      //  dotq->display();




      SP::SiconosVector q = d->q();

      //  -> get previous time step state
      SP::SiconosVector dotqold = d->dotqMemory()->getSiconosVector(0);
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      // *q = *qold + h*(theta * *v +(1.0 - theta)* *vold)
      double coeff = h * _theta;
      scal(coeff, *dotq, *q) ; // q = h*theta*v
      coeff = h * (1 - _theta);
      scal(coeff, *dotqold, *q, false); // q += h(1-theta)*vold
      *q += *qold;
      //  cout<<"new q before normalizing"<<endl;
      //  q->display();

      //q[3:6] must be normalized
      double normq = sqrt(q->getValue(3) * q->getValue(3) + q->getValue(4) * q->getValue(4) + q->getValue(5) * q->getValue(5) + q->getValue(6) * q->getValue(6));
      //  printf("-> normq : %f\n",normq);
      normq = 1 / normq;
      q->setValue(3, q->getValue(3) * normq);
      q->setValue(4, q->getValue(4) * normq);
      q->setValue(5, q->getValue(5) * normq);
      q->setValue(6, q->getValue(6) * normq);
      d->updateT();
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
    cout << "--> W of dynamical system number " << (*it)->number() << ": " << endl;
    if (WMap[*it]) WMap[*it]->display();
    else cout << "-> NULL" << endl;
    cout << "--> and corresponding theta is: " << _theta << endl;
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

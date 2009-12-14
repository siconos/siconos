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

#include "Lsodar.hpp"
#include "EventDriven.hpp"
#include "LagrangianLinearTIDS.hpp"
#include "BlockVector.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "Model.hpp"
#include "Topology.hpp"

using namespace std;

// ===== Out of class objects and functions =====

// global object and wrapping functions -> required for function plug-in and call in fortran routine.
SP::Lsodar global_object;

// This first function must have the same signature as argument F (arg 1) in DLSODAR (see opkdmain.f in Numerics)
extern "C" void Lsodar_f_wrapper(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot)
{
  return global_object->f(sizeOfX, time, x, xdot);
}

// Function to wrap g: same signature as argument G (arg 18) in DLSODAR (see opkdmain.f in Numerics)
extern "C" void Lsodar_g_wrapper(integer* nEq, doublereal* time, doublereal* x, integer* ng, doublereal* gOut)
{
  return global_object->g(nEq, time, x, ng, gOut);
}

// Function to wrap jacobianF: same signature as argument JAC (arg 16) in DLSODAR (see opkdmain.f in Numerics)
extern "C" void Lsodar_jacobianF_wrapper(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml, integer* mu,  doublereal* jacob, integer* nrowpd)
{
  return global_object->jacobianF(sizeOfX, time, x, ml, mu, jacob, nrowpd);
}

Lsodar::Lsodar(SP::OneStepIntegratorXML osiXML, SP::DynamicalSystemsSet dsList, SP::InteractionsSet interactionsList):
  OneStepIntegrator(OSI::LSODAR, osiXML, dsList, interactionsList)
{
  // local time discretisation is set by default to those of the simulation.
  intData.resize(9);
  sizeMem = 2;
}

Lsodar::Lsodar(SP::DynamicalSystem ds):
  OneStepIntegrator(OSI::LSODAR)
{
  // add ds in the set
  OSIDynamicalSystems->insert(ds);

  intData.resize(9);
  sizeMem = 2;
}

Lsodar::Lsodar(DynamicalSystemsSet& newDS):
  OneStepIntegrator(OSI::LSODAR, newDS)
{
  intData.resize(9);
  sizeMem = 2;
}

void Lsodar::setTol(integer newItol, SA::doublereal newRtol, SA::doublereal newAtol)
{
  //            The input parameters ITOL, RTOL, and ATOL determine
  //         the error control performed by the solver.  The solver will
  //         control the vector E = (E(i)) of estimated local errors
  //         in y, according to an inequality of the form
  //                     max-norm of ( E(i)/EWT(i) )   .le.   1,
  //         where EWT = (EWT(i)) is a vector of positive error weights.
  //         The values of RTOL and ATOL should all be non-negative.
  //         The following table gives the types (scalar/array) of
  //         RTOL and ATOL, and the corresponding form of EWT(i).
  //
  //            ITOL    RTOL       ATOL          EWT(i)
  //             1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  //             2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  //             3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  //             4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)

  intData[2] = newItol; // itol

  rtol = newRtol;
  atol = newAtol;
}

void Lsodar::updateData()
{
  // Used to update some data (iwork ...) when intData is modified.
  // Warning: it only checks sizes and possibly reallocate memory, but no values are set.

  unsigned int sizeTol = intData[0]; // size of rtol, atol ... If itol (intData[0]) = 1 => scalar else, vector of size neq (intData[0]).
  //  if(intData[0]==1) sizeTol = 1;
  //  else sizeTol = intData[0];

  rtol.reset(new doublereal[sizeTol]) ;    // rtol, relative tolerance

  atol.reset(new doublereal[sizeTol]) ;  // atol, absolute tolerance

  iwork.reset(new integer[intData[7]]);

  rwork.reset(new doublereal[intData[6]]);

  jroot.reset(new integer[intData[1]]);
}

void Lsodar::fillXWork(integer* sizeOfX, doublereal* x)
{
  unsigned int sizeX = (unsigned int)(*sizeOfX);
  for (unsigned int i = 0; i < sizeX ; ++i)
    (*xWork)(i) = x[i];
}

void Lsodar::computeRhs(double t)
{
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    (*it)->computeRhs(t);
}

void Lsodar::computeJacobianRhs(double t)
{
  DSIterator it;
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    (*it)->computeJacobianXRhs(t);
}

void Lsodar::f(integer* sizeOfX, doublereal* time, doublereal* x, doublereal* xdot)
{
  boost::static_pointer_cast<EventDriven>(simulationLink)->computeF(shared_from_this(), sizeOfX, time, x, xdot);
}

void Lsodar::g(integer* nEq, doublereal*  time, doublereal* x, integer* ng, doublereal* gOut)
{
  boost::static_pointer_cast<EventDriven>(simulationLink)->computeG(shared_from_this(), nEq, time, x, ng, gOut);
}

void Lsodar::jacobianF(integer* sizeOfX, doublereal* time, doublereal* x, integer* ml, integer* mu,  doublereal* jacob, integer* nrowpd)
{
  boost::static_pointer_cast<EventDriven>(simulationLink)->computeJacobianF(shared_from_this(), sizeOfX, time, x, jacob);
}

void Lsodar::initialize(SP::Simulation sim)
{
  OneStepIntegrator::initialize(sim);
  xWork.reset(new BlockVector());
  DSIterator it;
  string type;
  // initialize xWork with x values of the dynamical systems present in the set.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    xWork->insertPtr((*it)->x());

  //   Integer parameters for LSODAR are saved in vector intParam.
  //   The link with variable names in opkdmain.f is indicated in comments

  // 1 - Neq; x vector size.
  intData[0] = xWork->size();

  // 2 - Ng, number of constraints:
  intData[1] =  simulationLink->model()->nonSmoothDynamicalSystem()->topology()->numberOfConstraints();

  // 3 - Itol, itask, iopt
  intData[2] = 1; // itol, 1 if ATOL is a scalar, else 2 (ATOL array)
  intData[3] = 1; // itask, an index specifying the task to be performed. 1: normal computation.
  intData[5] = 0; // iopt: 0 if no optional input else 1.

  // 4 - Istate
  intData[4] = 1; // istate, an index used for input and output to specify the state of the calculation.
  // On input:
  //                 1: first call for the problem (initializations will be done).
  //                 2: means this is not the first call, and the calculation is to continue normally, with no change in any input
  //                    parameters except possibly TOUT and ITASK.
  //                 3:  means this is not the first call, and the calculation is to continue normally, but with
  //                     a change in input parameters other than TOUT and ITASK.
  // On output:
  //                 1: means nothing was done; TOUT = t and ISTATE = 1 on input.
  //                 2: means the integration was performed successfully, and no roots were found.
  //                 3: means the integration was successful, and one or more roots were found before satisfying the stop condition specified by ITASK. See JROOT.
  //                 <0: error. See table below, in integrate function output message.


  // 5 - lrw, size of rwork
  intData[6] = 22 + intData[0] * max(16, (int)intData[0] + 9) + 3 * intData[1];

  // 6 - liw, size of iwork
  intData[7] = 20 + intData[0];

  // 7 - JT, Jacobian type indicator
  intData[8] = 2;   // jt, Jacobian type indicator.
  //           1 means a user-supplied full (NEQ by NEQ) Jacobian.
  //           2 means an internally generated (difference quotient) full Jacobian (using NEQ extra calls to f per df/dx value).
  //           4 means a user-supplied banded Jacobian.
  //           5 means an internally generated banded Jacobian (using ML+MU+1 extra calls to f per df/dx evaluation).

  // memory allocation for doublereal*, according to intData values ...
  updateData();

  // Set atol and rtol values ...
  rtol[0] = MACHINE_PREC; // rtol
  atol[0] = MACHINE_PREC;  // atol

  // === Error handling in LSODAR===

  //   parameters: itol, rtol, atol.
  //   Control vector E = (E(i)) of estimated local errors in y:
  //   max-norm of ( E(i)/EWT(i) )< 1
  //   EWT = (EWT(i)) vector of positive error weights.
  //   The values of RTOL and ATOL should all be non-negative.
  //
  //  ITOL    RTOL       ATOL          EWT(i)
  //   1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL
  //   2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i)
  //   3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL
  //   4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i)

}

void Lsodar::integrate(double& tinit, double& tend, double& tout, int& istate)
{
  // For details on DLSODAR parameters, see opkdmain.f in Numerics/src/odepack

  doublereal tend_DR = tend  ;       // next point where output is desired (different from t!)
  doublereal tinit_DR = tinit;       // current (starting) time

  // === Pointers to function ===
  //  --> definition and initialisation thanks to wrapper:
  global_object = boost::static_pointer_cast<Lsodar>(shared_from_this()); // Warning: global object must be initialized to current one before pointers to function initialisation.

  // function to compute the righ-hand side of xdot = f(x,t) + Tu
  fpointer pointerToF = Lsodar_f_wrapper;

  // function to compute the Jacobian/x of the rhs.
  jacopointer pointerToJacobianF = Lsodar_jacobianF_wrapper; // function to compute the Jacobian/x of the rhs.

  // function to compute the constraints
  gpointer pointerToG;
  pointerToG = Lsodar_g_wrapper; // function to compute the constraints

  // === LSODAR CALL ===

  SP::SimpleVector xtmp(new SimpleVector(*xWork)); // A copy of xWork is required since at the time, there are no contiguous values in memory for BlockVectors.
  if (istate == 3)
  {
    istate = 1; // restart TEMPORARY
  }

  intData[4] = istate;

  F77NAME(dlsodar)(pointerToF, &(intData[0]), &(*xtmp)(0), &tinit_DR, &tend_DR, &(intData[2]), rtol.get(), atol.get(), &(intData[3]), &(intData[4]), &(intData[5]), rwork.get(), &(intData[6]), iwork.get(), &(intData[7]), pointerToJacobianF, &(intData[8]), pointerToG, &(intData[1]), jroot.get());

  // jroot: jroot[i] = 0 if g(i) has a root at t, else jroot[i] = 0.

  // === Post ===
  if (intData[4] < 0) // if istate < 0 => LSODAR failed
  {
    cout << "LSodar::integrate(...) failed - Istate = " << intData[4] << endl;
    cout << " -1 means excess work done on this call (perhaps wrong JT)." << endl;
    cout << " -2 means excess accuracy requested (tolerances too small)." << endl;
    cout << " -3 means illegal input detected (see printed message)." << endl;
    cout << " -4 means repeated error test failures (check all inputs)." << endl;
    cout << " -5 means repeated convergence failures (perhaps bad Jacobian supplied or wrong choice of JT or tolerances)." << endl;
    cout << " -6 means error weight became zero during problem. (Solution component i vanished, and ATOL or ATOL(i) = 0.)" << endl;
    cout << " -7 means work space insufficient to finish (see messages)." << endl;
    RuntimeException::selfThrow("Lsodar, integration failed");
  }

  *xWork = *xtmp;
  istate = intData[4];
  tout  = tinit_DR; // real ouput time
  tend  = tend_DR; // necessary for next start of DLSODAR


  if (istate == 3)
  {
    //      std:: cout << "ok\n";
    assert(true);
  }

  //  tinit = tinit_DR;
}


void Lsodar::updateState(unsigned int level)
{
  // Compute all required (ie time-dependent) data for the DS of the OSI.
  DSIterator it;

  if (level == 1) // ie impact case: compute velocity
  {
    for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
    {
      SP::LagrangianDS lds = boost::static_pointer_cast<LagrangianDS>(*it);
      lds->computePostImpactVelocity();
    }
  }
  else if (level == 2)// compute acceleration ie RHS and its jacobian.
  {
    double time = simulationLink->model()->currentTime();
    for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
      (*it)->update(time);
  }
  else RuntimeException::selfThrow("Lsodar::updateState(index), index is out of range. Index = " + level);
}

void Lsodar::display()
{
  OneStepIntegrator::display();
  cout << " --- > Lsodar specific values: " << endl;
  cout << "Number of equations: " << intData[0] << endl;
  cout << "Number of constraints: " << intData[1] << endl;
  cout << "itol, itask, istate, iopt, lrw, liw, jt: (for details on what are these variables see opkdmain.f)" << endl;
  cout << intData[2] << ", " << intData[3] << ", " << intData[4] << ", " << intData[5] << ", " << intData[6]  << ", " << intData[7]  << ", " << intData[8] << endl;
  cout << "====================================" << endl;

}

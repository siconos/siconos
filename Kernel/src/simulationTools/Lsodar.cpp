/* Siconos version 1.0, Copyright INRIA 2005.
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

#include "Lsodar.h"
using namespace std;

// ===== Out of class objects and functions =====

// global object and wrapping functions -> required for function plug-in and call in fortran routine.

Lsodar* global_object;

// This first function must have the same signature as argument F (arg 1) in DSLODAR (see opkdmain.f in Numerics)
extern "C" void Lsodar_f_wrapper(integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{
  return global_object->f(sizeOfX, time, x, xdot);
}

// Function to wrap g: same signature as argument G (arg 18) in DSLODAR (see opkdmain.f in Numerics)
extern "C" void Lsodar_g_wrapper(integer * nEq, doublereal * time, doublereal* x, integer* ng, doublereal * gOut)
{
  return global_object->g(nEq, time, x, ng, gOut);
}

// Function to wrap jacobianF: same signature as argument JAC (arg 16) in DSLODAR (see opkdmain.f in Numerics)
extern "C" void Lsodar_jacobianF_wrapper(integer * sizeOfX, doublereal * time, doublereal * x, integer* ml, integer * mu,  doublereal * jacob, integer * nrowpd)
{
  return global_object->jacobianF(sizeOfX, time, x, ml, mu, jacob, nrowpd);
}

// ===== Lsodar methods =====

Lsodar::Lsodar(): OneStepIntegrator(), localTimeDiscretisation(NULL), isLocalTimeDiscretisationAllocatedIn(false), iwork(NULL)
{
  integratorType = LSODAR_INTEGRATOR;
  intData.resize(9);
  doubleData.resize(4);
}

Lsodar::Lsodar(OneStepIntegratorXML* osiXML):
  OneStepIntegrator(osiXML), localTimeDiscretisation(NULL), isLocalTimeDiscretisationAllocatedIn(false), iwork(NULL)
{
  integratorType = LSODAR_INTEGRATOR;
  // local time discretisasation is set by default to those of the strategy.
  localTimeDiscretisation = timeDiscretisation;
  intData.resize(9);
  doubleData.resize(4);
}

Lsodar::Lsodar(TimeDiscretisation* td, DynamicalSystem* ds):
  OneStepIntegrator(td, ds), localTimeDiscretisation(NULL), isLocalTimeDiscretisationAllocatedIn(false), iwork(NULL)
{
  integratorType = LSODAR_INTEGRATOR;
  // local time discretisasation is set by default to those of the strategy.
  localTimeDiscretisation = timeDiscretisation;
  intData.resize(9);
  doubleData.resize(4);
}

Lsodar::~Lsodar()
{
  global_object = NULL;
  if (isLocalTimeDiscretisationAllocatedIn) delete localTimeDiscretisation;
  localTimeDiscretisation = NULL;
  vector<doublereal*>::iterator it;
  for (it = doubleData.begin(); it != doubleData.end(); ++it)
  {
    if (*it != NULL) delete *it;
    *it = NULL;
  }
  if (iwork != NULL) delete iwork;
  iwork = NULL;
}

void Lsodar::setTimeDiscretisationPtr(TimeDiscretisation* td)
{
  if (isLocalTimeDiscretisationAllocatedIn) delete localTimeDiscretisation;
  localTimeDiscretisation = td;
  isLocalTimeDiscretisationAllocatedIn = false;
}

void Lsodar::setIntData(const std::vector<integer>& newVector)
{
  if (newVector.size() != intData.size())
    RuntimeException::selfThrow("Lsodar::setIntData(vector), inconsistent sizes between input vector and class object member intData.");
  intData = newVector;
  // update doubleData and iwork according to new value of intData.
  updateData();
}

void Lsodar::setDoubleData(const std::vector<doublereal*>& newVector)
{
  // Check size coherence ...
  if (newVector.size() != doubleData.size())
    RuntimeException::selfThrow("Lsodar::setDoubleData(vector), inconsistent sizes between input vector and class object member doubleData.");
  long int sizeTol; // size of rtol, atol ... If itol (intData[2]) = 1 => scalar else vector of size neq (intData[0]).
  if (intData[2] == 1) sizeTol = 1;
  else sizeTol = intData[0];
  long int size1 = sizeof(newVector[0]) / sizeof(doublereal);
  long int size2 = sizeof(newVector[1]) / sizeof(doublereal);
  long int size3 = sizeof(newVector[2]) / sizeof(doublereal);
  long int size4 = sizeof(newVector[3]) / sizeof(doublereal);
  if (size1 != size2 || size1 != sizeTol || size3 != intData[6] || size4 != intData[1])
    RuntimeException::selfThrow("Lsodar::setDoubleData(vector), inconsistent sizes between an element of input vector and one of the class object member doubleData.");
  // set values by copy ...
  for (unsigned int i = 0; i < doubleData.size(); ++i)
    *(doubleData[i]) = *(newVector[i]);
}

void Lsodar::setIwork(integer* newValue)
{
  long int size = sizeof(newValue) / sizeof(integer);
  if (size != intData[7])
    RuntimeException::selfThrow("Lsodar::setIwork(integer*), inconsistent size between input and iwork.");
  // set by copy:
  for (long int i = 0; i < size; ++i)
    iwork[i] = newValue[i];
}

void Lsodar::updateData()
{
  // Used to update some data (ie doubleData, iwork) when intData is modified.
  // Warning: it only checks sizes and possibly reallocate memory, but no values are set for doubleData or iwork.

  unsigned int sizeTol; // size of rtol, atol ... If itol (intData[0]) = 1 => scalar else, vector of size neq (intData[0]).
  if (intData[0] == 1) sizeTol = 1;
  else sizeTol = intData[0];
  // Allocate memory for iwork
  if (iwork != NULL) delete iwork;
  iwork = new integer[intData[7]];
  // Allocate memory for doubleData ...
  if (doubleData[0] != NULL) delete doubleData[0];
  doubleData[0] = new doublereal[sizeTol] ;    // rtol, relative tolerance
  if (doubleData[1] != NULL) delete doubleData[1];
  doubleData[1] = new doublereal[sizeTol];  // atol, absolute tolerance
  if (doubleData[2] != NULL) delete doubleData[2];
  doubleData[2] = new doublereal[intData[6]]; // rwork
  if (doubleData[3] != NULL) delete doubleData[3];
  doubleData[3] = new doublereal[intData[1]]; // jroot
}

void Lsodar::f(integer * sizeOfX, doublereal * time, doublereal * x, doublereal * xdot)
{
  unsigned int size = *sizeOfX;
  SimpleVector *xtmp = new SimpleVector(size) ;

  // copy x in a temporary SimpleVector, to set x in Dynamical system.
  for (unsigned int i = 0; i < size; i++) /// Warning: copy !!
    (*xtmp)(i) = x[i];
  ds->setX(*xtmp);

  // Compute the vector field (=f) for the current ds
  double t = *time;
  ds->computeVectorField(t);

  // Save xdot values from dynamical system into current xdot (in-out parameter)
  SiconosVector * xtmp2 = ds->getXDotPtr();
  for (unsigned int i = 0; i < size; i++) /// Warning: copy !!
    xdot[i] = (*xtmp2)(i);

  delete xtmp;
  xtmp2 = NULL;
}

void Lsodar::g(integer * nEq, doublereal * time, doublereal* x, integer * ng, doublereal * gOut)
{}

void Lsodar::jacobianF(integer *sizeOfX, doublereal *time, doublereal *x, integer* ml, integer *mu,  doublereal *jacob, integer *nrowpd)
{
  unsigned int size = *sizeOfX;
  SimpleVector *xtmp = new SimpleVector(size) ;;

  // copy x in a temporary SimpleVector, to set x in Dynamical system.
  for (unsigned int i = 0; i < size; i++) /// Warning: copy !!
    (*xtmp)(i) = x[i];
  ds->setX(*xtmp);

  // Compute the vector field (=f) for the current ds
  double t = *time;
  ds->computeJacobianX(t);

  // Save jacobianX values from dynamical system into current jacob (in-out parameter)
  SiconosMatrix * jacotmp = ds->getJacobianXPtr();

  unsigned int k = 0;
  for (unsigned int j = 0; j < size; j++) /// Warning: copy !!
  {
    for (unsigned i = 0 ; i < size ; i++)
    {
      jacob[k] = (*jacotmp)(i, j);
      k++;
    }
  }
  delete xtmp;
}

void Lsodar::initialize()
{
  ds->initMemory(sizeMem);
  // check that all data (int, double and iwork) have been filled in.  ??

}

void Lsodar::computeFreeState() // useless??
{
  IN("Lsodar::computeFreeState\n");
  integrate();
  OUT("Lsodar::computeFreeState\n");
}

void Lsodar::integrate()
{
  IN("Lsodar::integrate\n");
  SiconosVector * y = ds->getXPtr();

  // get current LOCAL time discretisation vector;
  SimpleVector * tk = localTimeDiscretisation->getTkPtr();
  // get current step number
  unsigned int k = localTimeDiscretisation->getK();

  doublereal tout = (*tk)(k);             // next point where output is desired (different from t!)
  doublereal t = (*tk)(k - 1);            // current time

  //   Integer parameters for LSODAR are saved in vector intParam.
  //   The link with variable names in opkdmain.f is indicated in comments
  intData[0] =  y->size();  // neq, number of equations, ie dim of y
  intData[1] = 0 ;  // ng, number of constraints
  intData[2] = 1; // itol, 1 or 2 according as ATOL (below) is a scalar or an array.
  intData[3] = 1; // itask
  intData[4] = 1; // istate
  intData[5] = 0; // iopt
  intData[6] = 22 + intData[0] * max(16, (int)intData[0] + 9) + 3 * intData[1]; // lrw
  intData[7] = 20 + intData[0];  // liw
  intData[8] = 2;   // jt
  // update memory size for doubleData and iwork according to intData values ...
  updateData();

  //   Doublereal parameters for LSODAR are saved in vector doubleData.
  //   The link with variable names in opkdmain.f is indicated in comments
  *(doubleData[0]) = 0.0;
  *(doubleData[1]) = 1.0e-6;

  // Pointers to function definition and initialisation thanks to wrapper:
  global_object = this; // Warning: global object must be initialized to current one before pointers to function initialisation.
  fpointer pointerToF = Lsodar_f_wrapper;
  jacopointer pointerToJacobianF = Lsodar_jacobianF_wrapper;
  gpointer pointerToG;
  pointerToG = Lsodar_g_wrapper;

  F77NAME(dlsodar)(pointerToF, &(intData[0]), &(*y)(0), &t, &tout, &(intData[2]), doubleData[0], doubleData[1], &(intData[3]), &(intData[4]), &(intData[5]), doubleData[2],
                   &(intData[6]), iwork, &(intData[7]), pointerToJacobianF, &(intData[8]), pointerToG, &(intData[1]), doubleData[3]);
  //   integer nqu;
  //   doublereal hu;
  //   hu  = doubleData[2][10];
  //   nqu = iwork[13];
  //   cout << t << "     " << (*y)(0)  << "     " << (*y)(1) << "     " << nqu << "     " << hu << "     " << endl ;
  //   if (istate<0)
  //     break;
  //   iopar = iout%2;
  //   if (iopar!=0)
  //     tout = tout + dt;
  //   else
  //     {
  //       er = abs((*y)(0))/atol;
  //       ero = max(ero,er);
  //       if (er>1000)
  //  {
  //    cout <<" Warning: error exceeds 1000 * tolerance" << endl;
  //    nerr = nerr + 1;
  //  }
  //  tout = tout + (*tk)(k+2)-(*tk)(k+1);

  // update local time discretisation
  //      localTimeDiscretisation ->setK(k+1);

  if (intData[2] < 0) RuntimeException::selfThrow("Lsodar, integration failed (see opkdmain.f for details about istate value), istate = " + intData[2]);

  OUT("Lsodar::integrate\n");
}

void Lsodar::updateState()
{}


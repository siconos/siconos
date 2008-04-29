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
#include "Moreau2.h"
#include "Simulation.h"
#include "Model.h"
#include "NonSmoothDynamicalSystem.h"
#include "TimeDiscretisation.h"

#include "LagrangianLinearTIDS.h"
#include "FirstOrderLinearTIDS.h"

using namespace std;

// --- constructor from a minimum set of data ---
Moreau2::Moreau2(DynamicalSystem* newDS, double newTheta, Simulation* newS): Moreau(newDS, newTheta, newS)
{
}


Moreau2::~Moreau2()
{
}

// SiconosVector * Moreau2::getFfree(FirstOrderLinearDS *d){
//   return workX[d];
// }

void Moreau2::computeFreeState()
{

  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  double t = simulationLink->getNextTime(); // End of the time step
  double told = simulationLink->getStartingTime(); // Beginning of the time step
  double h = t - told; // time step length

  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SimpleVector *rold = static_cast<SimpleVector*>(d->getRMemoryPtr()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it; // Iterator through the set of DS.

  DynamicalSystem* ds; // Current Dynamical System.
  SiconosMatrix * W; // W Moreau matrix of the current DS.
  SiconosMatrix * M; // W Moreau matrix of the current DS.
  string dsType ; // Type of the current DS.
  double theta; // Theta parameter of the current ds.
  for (it = OSIDynamicalSystems->begin(); it != OSIDynamicalSystems->end(); ++it)
  {
    ds = *it; // the considered dynamical system
    theta = thetaMap[ds]; // Its theta parameter
    dsType = ds->getType(); // Its type
    W = WMap[ds]; // Its W Moreau matrix of iteration.

    if (dsType == FOLDS)
    {
      // fFree = h(1-theta)Ax_i +Mx_i +h*theta*b_i+1 + h(1-theta)*b_i

      // IN to be updated at current time: W, b
      // IN at told: A, b, x
      // IN, not time dependant: M

      FirstOrderLinearDS *d = static_cast<FirstOrderLinearDS*>(ds);

      SiconosVector *xfree = workX[d];

      // x value at told
      //SiconosVector *xold = d->getXMemoryPtr()->getSiconosVector(0);
      SiconosVector *xold = d->getXPtr();

      // If M not equal to identity matrix
      SiconosMatrix * M = d->getMPtr();
      if (M != NULL)
        prod(*M, *xold, *xfree); // fFree = M*xi
      else
        *xfree = *xold;

      SiconosMatrix *A = d->getAPtr();
      if (A != NULL)
      {
        d->computeA(told);
        double coeff = h * (1 - theta);
        prod(coeff, *A, *xold, *xfree, false);
        // fFree += h(1-theta)A_i*x_i
      }
      SiconosVector *b = d->getBPtr();
      if (b != NULL)
      {
        // fFree += h(1-theta)*bi + h*theta*bi+1
        //        d->computeB(told); // bi
        scal(h * (1.0 - theta), *d->getBPtr(), *xfree, false);
        d->computeB(t); // bi+1
        scal(h * theta, *d->getBPtr(), *xfree, false);
      }

      // -- Update W --
      computeW(t, d);
    }
    // 2bis - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == FOLTIDS)
    {
      // fFree = h(1-theta)Ax_i + Mx_i +hbt(i+1)

      // IN to be updated at current time: none
      // IN at told: x
      // IN, not time dependant: A,b,W

      FirstOrderLinearDS *d = static_cast<FirstOrderLinearTIDS*>(ds);
      M = d->getMPtr();
      SiconosVector *xfree = workX[d];
      // x value at told
      SiconosVector *xold = d->getXMemoryPtr()->getSiconosVector(0);

      SiconosMatrix *A = d->getAPtr();
      if (A != NULL)
        prod(h * (1 - theta), *A, *xold, *xfree, true); // ffree = h*(1-theta)*A*xi
      else
        xfree->zero();

      SiconosVector *b = d->getBPtr();
      if (b != NULL)
        scal(h, *b, *xfree, false); // ffree += hb
      if (M != NULL)
        prod(*M, *xold, *xfree, false); // ffree += M*xi
    }
    else
    {
      RuntimeException::selfThrow("Moreau2::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
    }
  }
}


void Moreau2::updateState(unsigned int level)
{

}
Moreau2* Moreau2::convert(OneStepIntegrator* osi)
{
  Moreau2* moreau = dynamic_cast<Moreau2*>(osi);
  return moreau;
}

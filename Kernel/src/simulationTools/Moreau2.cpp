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
Moreau2::Moreau2(SP::DynamicalSystem newDS, double newTheta, SP::Simulation newS): Moreau(newDS, newTheta, newS)
{
  integratorType = "Moreau2";
}
Moreau2::Moreau2(DynamicalSystemsSet& newDS, double newTheta, SP::Simulation newS): Moreau(newDS, newTheta, newS)
{
  integratorType = "Moreau2";
}




Moreau2::~Moreau2()
{
}

SP::SiconosVector  Moreau2::getWorkX(DynamicalSystem *d)
{
  return workX[d];
}

void Moreau2::computeFreeState()
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
  SP::SiconosMatrix  W; // W Moreau matrix of the current DS.
  SP::SiconosMatrix  M; // W Moreau matrix of the current DS.
  DSTYPES dsType ; // Type of the current DS.
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

      SP::SiconosVector ffree = workX[d];

      // x value at told
      //SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0);
      SP::SiconosVector xold = d->getXPtr();

      // If M not equal to identity matrix
      SP::SiconosMatrix  M = d->getMPtr();
      if (M != NULL)
        prod(*M, *xold, *ffree); // fFree = M*xi
      else
        *ffree = *xold;

      SP::SiconosMatrix A = d->getAPtr();
      if (A != NULL)
      {
        d->computeA(told);
        double coeff = h * (1 - theta);
        prod(coeff, *A, *xold, *ffree, false);
        // fFree += h(1-theta)A_i*x_i
      }
      SP::SiconosVector b = d->getBPtr();
      if (b != NULL)
      {
        // fFree += h(1-theta)*bi + h*theta*bi+1
        //        d->computeB(told); // bi
        scal(h * (1.0 - theta), *d->getBPtr(), *ffree, false);
        d->computeB(t); // bi+1
        scal(h * theta, *d->getBPtr(), *ffree, false);
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
      SP::SiconosVector ffree = workX[d];
      // x value at told
      SP::SiconosVector xold = d->getXMemoryPtr()->getSiconosVector(0);

      SP::SiconosMatrix A = d->getAPtr();
      if (A != NULL)
        prod(h * (1 - theta), *A, *xold, *ffree, true); // ffree = h*(1-theta)*A*xi
      else
        ffree->zero();
      SP::SiconosVector b = d->getBPtr();
      if (b != NULL)
        scal(h, *b, *ffree, false); // ffree += hb
      if (M != NULL)
        prod(*M, *xold, *ffree, false); // ffree += M*xi
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


      // fFree =  W v_k,i+1 -M(q_k,i+1)(v_k,i+1- v_i) - h*theta*fL(t,v_k,i+1, q_k,i+1) - h*(1-theta)*fL(ti,vi,qi)

      // -- Convert the DS into a Lagrangian one.
      LagrangianDS* d = static_cast<LagrangianDS*>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0);
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0); // vol =v_i

      // --- ResiduFree computation ---
      // vFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector ffree = workX[d];

      // -- Update W --
      // Note: during computeW, mass and jacobians of fL will be computed/
      computeW(t, d);

      SP::SiconosMatrix M = d->getMassPtr();
      SP::SiconosVector v = d->getVelocityPtr(); // v = v_k,i+1
      prod(*M, (*v - *vold), *ffree); // ffree = M(v - vold)

      *ffree *= -1.0;
      if (d->getFLPtr() != NULL) // if fL exists
      {
        // computes fL(ti,vi,qi)
        d->computeFL(told, qold, vold);
        double coef = h * (1 - theta);
        // ffree += coef * fL_i
        scal(coef, *d->getFLPtr(), *ffree, false);

        // computes fL(ti+1, v_k,i+1, q_k,i+1) = fL(t,v,q)
        d->computeFL(t);
        coef = h * theta;
        // ffree += coef * fL_k,i+1
        scal(coef, *d->getFLPtr(), *ffree, false);
      }

      SP::SiconosVector  ftmp = new SimpleVector(*ffree);
      prod(*W, (*v), *ftmp);
      *ffree += *ftmp;
      delete(ftmp);
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == LLTIDS)
    {
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // fFree = +W v_i + (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 + h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      LagrangianLinearTIDS* d = static_cast<LagrangianLinearTIDS*>(ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->getQMemoryPtr()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->getVelocityMemoryPtr()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---
      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SP::SiconosVector ffree = workX[d];

      SP::SiconosMatrix M = d->getMassPtr();
      SP::SiconosVector v = d->getVelocityPtr(); // v = v_k,i+1
      prod(*W, (*v - *vold), *ffree); // ffree = W (vold)

      // vFree pointer is used to compute and save ResiduFree in this first step.
      double coeff;

      SP::SiconosMatrix  C = d->getCPtr();
      if (C != NULL)
        prod(-h, *C, *vold, *ffree, false); // ffree += -h*C*vi

      SP::SiconosMatrix  K = d->getKPtr();
      if (K != NULL)
      {
        coeff = -h * h * theta;
        prod(coeff, *K, *vold, *ffree, false); // ffree += -h^2*theta*K*vi
        prod(-h, *K, *qold, *ffree, false); // ffree += -h*K*qi
      }

      SP::SiconosVector  Fext = d->getFExtPtr();
      if (Fext != NULL)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - theta);
        scal(coeff, *Fext, *ffree, false); // ffree += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * theta;
        scal(coeff, *Fext, *ffree, false); // ffree += h*theta * fext(ti+1)
      }

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

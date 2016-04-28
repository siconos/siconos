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
#include "MoreauJeanOSI2.hpp"
#include "Simulation.hpp"
#include "Model.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "TimeDiscretisation.hpp"

#include "LagrangianLinearTIDS.hpp"
#include "FirstOrderLinearTIDS.hpp"


using namespace RELATION;

// --- constructor from a minimum set of data ---
MoreauJeanOSI2::MoreauJeanOSI2(double theta): MoreauJeanOSI(theta)
{
  _integratorType = OSI::MOREAUJEANOSI2;
}

MoreauJeanOSI2::~MoreauJeanOSI2()
{
}

void MoreauJeanOSI2::computeFreeState()
{

  // This function computes "free" states of the DS belonging to this Integrator.
  // "Free" means without taking non-smooth effects into account.

  double t = _simulation->nextTime(); // End of the time step
  double told = _simulation->startingTime(); // Beginning of the time step
  double h = t - told; // time step length
  //h=0.0100000;
  // Operators computed at told have index i, and (i+1) at t.

  //  Note: integration of r with a theta method has been removed
  //  SiconosVector *rold = static_cast<SiconosVector*>(d->rMemory()->getSiconosVector(0));

  // Iteration through the set of Dynamical Systems.
  //
  DSIterator it; // Iterator through the set of DS.

  SP::DynamicalSystem ds; // Current Dynamical System.
  SP::SiconosMatrix  W; // W MoreauJeanOSI matrix of the current DS.
  SP::SiconosMatrix  M; // W MoreauJeanOSI matrix of the current DS.
  Type::Siconos dsType ; // Type of the current DS.

  DynamicalSystemsGraph::VIterator dsi, dsend;
  for (std11::tie(dsi, dsend) = _dynamicalSystemsGraph->vertices(); dsi != dsend; ++dsi)
  {
    if (!checkOSI(dsi)) continue;

    SP::DynamicalSystem ds = _dynamicalSystemsGraph->bundle(*dsi);
    dsType = Type::value(*ds); // Its type
    W = _dynamicalSystemsGraph->properties(*dsi).W; // Its W MoreauJeanOSI matrix of iteration.

    if (dsType == Type::FirstOrderLinearDS)
    {
      // fFree = h(1-theta)Ax_i +Mx_i +h*theta*b_i+1 + h(1-theta)*b_i

      // IN to be updated at current time: W, b
      // IN at told: A, b, x
      // IN, not time dependant: M

      SP::FirstOrderLinearDS d = std11::static_pointer_cast<FirstOrderLinearDS>(ds);

      SP::SiconosVector ffree = d->workspace(DynamicalSystem::free);

      // x value at told
      //SP::SiconosVector xold = d->xMemory()->getSiconosVector(0);
      SP::SiconosVector xold = d->x();

      // If M not equal to identity matrix
      SP::SiconosMatrix  M = d->M();
      if (M)
        prod(*M, *xold, *ffree); // fFree = M*xi
      else
        *ffree = *xold;

      SP::SiconosMatrix A = d->A();
      if (A)
      {
        d->computeA(told);
        double coeff = h * (1 - _theta);
        prod(coeff, *A, *xold, *ffree, false);
        // fFree += h(1-theta)A_i*x_i
      }
      SP::SiconosVector b = d->b();
      if (b)
      {
        // fFree += h(1-theta)*bi + h*theta*bi+1
        //        d->computeb(told); // bi
        scal(h * (1.0 - _theta), *d->b(), *ffree, false);
        d->computeb(t); // bi+1
        scal(h * _theta, *d->b(), *ffree, false);
      }

      // -- Update W --
      computeW(t, d, *W);
    }
    // 2bis - First Order Linear Systems with Time Invariant coefficients
    else if (dsType == Type::FirstOrderLinearTIDS)
    {
      // fFree = h(1-theta)Ax_i + Mx_i +hbt(i+1)

      // IN to be updated at current time: none
      // IN at told: x
      // IN, not time dependant: A,b,W

      SP::FirstOrderLinearDS d = std11::static_pointer_cast<FirstOrderLinearTIDS>(ds);
      M = d->M();
      SP::SiconosVector ffree = d->workspace(DynamicalSystem::free);
      // x value at told
      SP::SiconosVector xold = d->xMemory()->getSiconosVector(0);

      SP::SiconosMatrix A = d->A();
      if (A)
        prod(h * (1 - _theta), *A, *xold, *ffree, true); // ffree = h*(1-theta)*A*xi
      else
        ffree->zero();
      SP::SiconosVector b = d->b();
      if (b)
        scal(h, *b, *ffree, false); // ffree += hb
      if (M)
        prod(*M, *xold, *ffree, false); // ffree += M*xi
    }
    // 3 - Lagrangian Non Linear Systems
    else if (dsType == Type::LagrangianDS)
    {
      // IN to be updated at current time: W, M, q, v, fL
      // IN at told: qi,vi, fLi

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // Index k stands for Newton iteration and thus corresponds to the last computed
      // value, ie the one saved in the DynamicalSystem.
      // "i" values are saved in memory vectors.


      // fFree =  W v_k,i+1 -M(q_k,i+1)(v_k,i+1- v_i) - h*theta*forces(t,v_k,i+1, q_k,i+1) - h*(1-theta)*forces(ti,vi,qi)

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianDS d = std11::static_pointer_cast<LagrangianDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0);
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); // vol =v_i

      // --- ResiduFree computation ---
      // vFree pointer is used to compute and save ResiduFree in this first step.
      SP::SiconosVector ffree = d->workspace(DynamicalSystem::free);

      // -- Update W --
      // Note: during computeW, mass and jacobians of fL will be computed/
      computeW(t, d, *W);

      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      prod(*M, (*v - *vold), *ffree); // ffree = M(v - vold)

      *ffree *= -1.0;
      if (d->forces()) // if fL exists
      {
        // computes forces(ti,vi,qi)
        d->computeForces(told, qold, vold);
        double coef = h * (1 - _theta);
        // ffree += coef * fL_i
        scal(coef, *d->forces(), *ffree, false);

        // computes forces(ti+1, v_k,i+1, q_k,i+1) = forces(t,v,q)
        d->computeForces(t);
        coef = h * _theta;
        // ffree += coef * fL_k,i+1
        scal(coef, *d->forces(), *ffree, false);
      }

      SP::SiconosVector  ftmp(new SiconosVector(*ffree));
      prod(*W, (*v), *ftmp);
      *ffree += *ftmp;
    }
    // 4 - Lagrangian Linear Systems
    else if (dsType == Type::LagrangianLinearTIDS)
    {
      // IN to be updated at current time: Fext
      // IN at told: qi,vi, fext
      // IN constants: K,C

      // Note: indices i/i+1 corresponds to value at the beginning/end of the time step.
      // "i" values are saved in memory vectors.

      // fFree = +W v_i + (-h*C -h^2*theta*K)*vi - h*K*qi + h*theta * Fext_i+1 + h*(1-theta)*Fext_i

      // -- Convert the DS into a Lagrangian one.
      SP::LagrangianLinearTIDS d = std11::static_pointer_cast<LagrangianLinearTIDS> (ds);

      // Get state i (previous time step) from Memories -> var. indexed with "Old"
      SP::SiconosVector qold = d->qMemory()->getSiconosVector(0); // qi
      SP::SiconosVector vold = d->velocityMemory()->getSiconosVector(0); //vi

      // --- ResiduFree computation ---
      // Velocity free and residu. vFree = RESfree (pointer equality !!).
      SP::SiconosVector ffree = d->workspace(DynamicalSystem::free);

      SP::SiconosMatrix M = d->mass();
      SP::SiconosVector v = d->velocity(); // v = v_k,i+1
      prod(*W, (*v - *vold), *ffree); // ffree = W (vold)

      // vFree pointer is used to compute and save ResiduFree in this first step.
      double coeff;

      SP::SiconosMatrix  C = d->C();
      if (C)
        prod(-h, *C, *vold, *ffree, false); // ffree += -h*C*vi

      SP::SiconosMatrix  K = d->K();
      if (K)
      {
        coeff = -h * h * _theta;
        prod(coeff, *K, *vold, *ffree, false); // ffree += -h^2*theta*K*vi
        prod(-h, *K, *qold, *ffree, false); // ffree += -h*K*qi
      }

      SP::SiconosVector  Fext = d->fExt();
      if (Fext)
      {
        // computes Fext(ti)
        d->computeFExt(told);
        coeff = h * (1 - _theta);
        scal(coeff, *Fext, *ffree, false); // ffree += h*(1-theta) * fext(ti)
        // computes Fext(ti+1)
        d->computeFExt(t);
        coeff = h * _theta;
        scal(coeff, *Fext, *ffree, false); // ffree += h*theta * fext(ti+1)
      }

    }



    else
    {
      RuntimeException::selfThrow("MoreauJeanOSI2::computeFreeState - not yet implemented for Dynamical system type: " + dsType);
    }
  }
}



void MoreauJeanOSI2::updateState(const unsigned int level)
{

}

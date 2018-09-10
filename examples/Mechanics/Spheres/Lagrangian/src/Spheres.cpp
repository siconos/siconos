/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file Disks.cpp

  Some Disks (2D), friction, and walls.
  Direct description of the model.  Simulation with
  a Time-Stepping scheme.
*/

// Siconos
#include <SiconosBodies.hpp>
#include <SiconosKernel.hpp>

#include "Spheres.hpp"

#include <SphereLDS.hpp>
#include <SpaceFilter.hpp>

#include <iostream>

using namespace std;

/* do nothing if solver does not converge */
void localCheckSolverOuput(int, Simulation*)
{};


// ================= Creation of the model =======================
void Spheres::init()
{
  SP::TimeDiscretisation timedisc_;
  SP::FrictionContact osnspb_;
  SP::NonSmoothDynamicalSystem nsds;

  // User-defined main parameters

  double t0 = 0;                   // initial computation time

  double T = std::numeric_limits<double>::infinity();

  double h = 0.01;                // time step
  double g = 9.81;

  double theta = 0.5;              // theta for MoreauJeanOSI integrator

  std::string solverName = "NSGS";

  // -----------------------------------------
  // --- Dynamical systems && interactions ---
  // -----------------------------------------


  double R;
  double m;

  try
  {

    // ------------
    // --- Init ---
    // ------------

    std::cout << "====> Model loading ..." << std::endl << std::endl;

    _plans.reset(new SimpleMatrix("plans.dat", true));

    SP::SiconosMatrix Spheres;
    Spheres.reset(new SimpleMatrix("spheres.dat", true));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new MoreauJeanOSI(theta));

    // -- Time discretisation --
    timedisc_.reset(new TimeDiscretisation(t0, h));

    // -- Model --
    nsds.reset(new NonSmoothDynamicalSystem(t0, T));

    // -- Simulation (for associate()) --
    _sim.reset(new TimeStepping(nsds, timedisc_));

    for (unsigned int i = 0; i < Spheres->size(0); i++)
    {
      R = Spheres->getValue(i, 3);
      m = Spheres->getValue(i, 4);

      SP::SiconosVector qTmp;
      SP::SiconosVector vTmp;

      qTmp.reset(new SiconosVector(NDOF));
      vTmp.reset(new SiconosVector(NDOF));
      vTmp->zero();
      (*qTmp)(0) = (*Spheres)(i, 0);
      (*qTmp)(1) = (*Spheres)(i, 1);
      (*qTmp)(2) = (*Spheres)(i, 2);

      (*qTmp)(3) = M_PI / 2;
      (*qTmp)(4) = M_PI / 4;
      (*qTmp)(5) = M_PI / 2;

      (*vTmp)(0) = 0;
      (*vTmp)(1) = 0;
      (*vTmp)(2) = 0;


      (*vTmp)(3) = 0;
      (*vTmp)(4) = 0;
      (*vTmp)(5) = 0;


      SP::LagrangianDS body;
      body.reset(new SphereLDS(R, m, std11::shared_ptr<SiconosVector>(qTmp), std11::shared_ptr<SiconosVector>(vTmp)));

      // -- Set external forces (weight) --
      SP::SiconosVector FExt;
      FExt.reset(new SiconosVector(NDOF));
      FExt->zero();
      FExt->setValue(2, -m * g);
      body->setFExtPtr(FExt);

      // associate the dynamical system with the one step integrator
      _sim->associate(osi, body);

      // add the dynamical system in the non smooth dynamical system
      nsds->insertDynamicalSystem(body);
    }

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- OneStepNsProblem --
    osnspb_.reset(new FrictionContact(3));

    osnspb_->numericsSolverOptions()->iparam[0] = 100; // Max number of
    // iterations
    osnspb_->numericsSolverOptions()->iparam[1] = 20; // compute error
    // iterations

    osnspb_->numericsSolverOptions()->iparam[4] = 2; // projection

    osnspb_->numericsSolverOptions()->dparam[0] = 1e-6; // Tolerance
    osnspb_->numericsSolverOptions()->dparam[2] = 1e-8; // Local tolerance


    osnspb_->setMaxSize(16384);       // max number of interactions
    osnspb_->setMStorageType(1);      // Sparse storage
    osnspb_->setNumericsVerboseMode(0); // 0 silent, 1 verbose
    osnspb_->setKeepLambdaAndYState(true); // inject previous solution

    _sim->insertNonSmoothProblem(osnspb_);
    //     simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, 0.8, 3));

    _playground.reset(new SpaceFilter(3, 6, _plans, _moving_plans));

    _playground->insertNonSmoothLaw(nslaw, 0, 0);

    _sim->insertInteractionManager(_playground);
  }

  catch (SiconosException e)
  {
    std::cout << e.report() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cout << "Exception caught in Spheres::init()" << std::endl;
    exit(1);
  }
}



// =========================== End of model definition ===========================

// ================================= Computation =================================


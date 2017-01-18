/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*!\file Spheres.cpp

  Some Spheres (3D), friction, and walls.
  Direct description of the model.  Simulation with
  a Time-Stepping scheme.
*/

// Siconos
#include <SiconosBodies.hpp>
#include <SiconosKernel.hpp>
#include <SphereNEDS.hpp>
#include <SphereNEDSPlanR.hpp>
#include <TimeStepping.hpp>
#include <FrictionContact.hpp>
#include <MoreauJeanOSI.hpp>
#include <NewtonImpactFrictionNSL.hpp>

#include <SpaceFilter.hpp>

//#include <Siconos/io/SiconosRestart.hpp>
//#include <Siconos/io/SiconosVTKOutput.hpp>
#include "Spheres.hpp"

using namespace std;

/* do nothing if solver does not converge */
void localCheckSolverOuput(int, Simulation*)
{};
//#define WITH_GENERIC_SOLVER

// ================= Creation of the model =======================
void Spheres::init()
{

  SP::TimeDiscretisation timedisc_;
  SP::Simulation simulation_;
#ifdef WITH_GENERIC_SOLVER
  SP::GenericMechanical osnspb_;
#else
  SP::FrictionContact osnspb_;
#endif
  // User-defined main parameters

  double t0 = 0;                   // initial computation time

  double T = std::numeric_limits<double>::infinity();

  double h = 0.005;                // time step
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

    // -- Model --
    _model.reset(new Model(t0, T));

    for (unsigned int i = 0; i < Spheres->size(0); i++)
    {
      R = Spheres->getValue(i, 3);
      m = Spheres->getValue(i, 4);

      SP::SiconosVector qTmp;
      SP::SiconosVector vTmp;

      qTmp.reset(new SiconosVector(7));
      vTmp.reset(new SiconosVector(6));
      qTmp->zero();
      vTmp->zero();

      (*qTmp)(0) = (*Spheres)(i, 0);
      (*qTmp)(1) = (*Spheres)(i, 1);
      (*qTmp)(2) = (*Spheres)(i, 2);

      (*qTmp)(3) = 1.;
      (*qTmp)(4) = 0.;
      (*qTmp)(5) = 0.;
      (*qTmp)(6) = 0.;

      SP::SimpleMatrix IMat(new SimpleMatrix(3, 3));
      (*IMat)(0, 0) = (*IMat)(1, 1) = (*IMat)(2, 2) = m * R * R * 2. / 5.;

      SP::NewtonEulerDS body;
      body.reset(new SphereNEDS(R, m, IMat, boost::shared_ptr<SiconosVector>(qTmp),
                                boost::shared_ptr<SiconosVector>(vTmp)));

      // -- Set external forces (weight) --
      SP::SiconosVector FExt;
      FExt.reset(new SiconosVector(6)); //
      FExt->zero();
      //if (i)
      FExt->setValue(2, -m * g);
      body->setFExtPtr(FExt);

      // add the dynamical system in the non smooth dynamical system
      _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);
    }

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    timedisc_.reset(new TimeDiscretisation(t0, h));

    // -- OneStepNsProblem --
#ifdef WITH_GENERIC_SOLVER
    osnspb_.reset(new GenericMechanical(SICONOS_FRICTION_3D_ONECONTACT_QUARTIC));
#else
    osnspb_.reset(new FrictionContact(3));
    osnspb_->numericsSolverOptions()->iparam[0] = 1000; // Max number of
    // iterations
    osnspb_->numericsSolverOptions()->iparam[1] = 20; // compute error
    // iterations

    osnspb_->numericsSolverOptions()->iparam[4] = 2; // projection

    osnspb_->numericsSolverOptions()->dparam[0] = 1e-7; // Tolerance
    osnspb_->numericsSolverOptions()->dparam[2] = 1e-7; // Local tolerance


    osnspb_->setMaxSize(16384);       // max number of interactions
    osnspb_->setMStorageType(1);      // Sparse storage
    osnspb_->setNumericsVerboseMode(0); // 0 silent, 1 verbose
    osnspb_->setKeepLambdaAndYState(true); // inject previous solution
#endif
    simulation_.reset(new TimeStepping(timedisc_, osi, osnspb_));
    //     simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0.0, 0.0, 0.6, 3));

    _playground.reset(new SpaceFilter(3, 6, _model, _plans, _moving_plans));

    _playground->insert(nslaw, 0, 0);

    _model->setSimulation(simulation_);
    _model->initialize();

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


void Spheres::compute()
{

  try
  {

    _playground->buildInteractions(_model->currentTime());

    _model->simulation()->advanceToEvent();

    _model->simulation()->processEvents();

    //    output.update();

    //    output.write();

  }

  catch (SiconosException e)
  {
    std::cout << e.report() << std::endl;
  }
  catch (...)
  {
    std::cout << "Exception caught in SiconosBodies::compute()" << std::endl;
  }
}


// =========================== End of model definition ===========================

// ================================= Computation =================================


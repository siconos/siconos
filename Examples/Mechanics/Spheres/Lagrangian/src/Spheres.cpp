/* Siconos-sample , Copyright INRIA 2005-2011.
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

/*!\file Disks.cpp

  Some Disks (2D), friction, and walls.
  Direct description of the model.  Simulation with
  a Time-Stepping scheme.
*/

// Siconos
#include <SiconosKernel.hpp>

#include "Spheres.hpp"

#include <SphereLDS.hpp>
#include <SpaceFilter.hpp>

using namespace std;

/* do nothing if solver does not converge */
void localCheckSolverOuput(int, Simulation*)
{};


// ================= Creation of the model =======================
void Spheres::init()
{

  SP::TimeDiscretisation timedisc_;
  SP::Simulation simulation_;
  SP::FrictionContact osnspb_;


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

    // -- Model --
    _model.reset(new Model(t0, T));

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

      // add the dynamical system to the one step integrator
      osi->insertDynamicalSystem(body);

      // add the dynamical system in the non smooth dynamical system
      _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);

    }

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    timedisc_.reset(new TimeDiscretisation(t0, h));

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

    simulation_.reset(new TimeStepping(timedisc_));
    simulation_->insertIntegrator(osi);
    simulation_->insertNonSmoothProblem(osnspb_);
    //     simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, 0.8, 3));

    _playground.reset(new SpaceFilter(3, 6, _model, nslaw, _plans, _moving_plans));

    _model->initialize(simulation_);

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


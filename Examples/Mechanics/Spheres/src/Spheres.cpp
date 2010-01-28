/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
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
  Direct description of the model without XML input.  Simulation with
  a Time-Stepping scheme.
*/

// Siconos
#include <SiconosKernel.hpp>

#include "Spheres.hpp"

using namespace std;

/* do nothing if solver does not converge */
void localCheckSolverOuput(int, Simulation*)
{};


// ================= Creation of the model =======================
void Spheres::init()
{

  SP::TimeDiscretisation timedisc_;
  SP::Simulation simulation_;
  SP::NonSmoothDynamicalSystem nsds_;
  SP::FrictionContact osnspb_;

  DynamicalSystemsSet allDS_;
  InteractionsSet allInteractions_;



  // User-defined main parameters

  double t0 = 0;                   // initial computation time

  double T = 0.02;

  double h = 0.01;                // time step
  double g = 9.81;

  double theta = 0.5;              // theta for Moreau integrator

  std::string solverName = "NSGS";

  // -----------------------------------------
  // --- Dynamical systems && interactions ---
  // -----------------------------------------

  unsigned int j;
  int interCounter = 0;

  double R;
  double m;

  try
  {

    // ------------
    // --- Init ---
    // ------------

    std::cout << "====> Model loading ..." << std::endl << std::endl;

    plans_.reset(new SimpleMatrix("plans.dat", true));

    SP::SiconosMatrix Spheres;
    Spheres.reset(new SimpleMatrix("spheres.dat", true));

    for (unsigned int i = 0; i < Spheres->size(0); i++)
    {
      R = Spheres->getValue(i, 3);
      m = Spheres->getValue(i, 4);

      SP::SimpleVector qTmp;
      SP::SimpleVector vTmp;

      qTmp.reset(new SimpleVector(NDOF));
      vTmp.reset(new SimpleVector(NDOF));
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
      body.reset(new SphereLDS(R, m, boost::shared_ptr<SimpleVector>(qTmp), boost::shared_ptr<SimpleVector>(vTmp)));

      // -- Set external forces (weight) --
      SP::SimpleVector FExt;
      FExt.reset(new SimpleVector(NDOF));
      FExt->zero();
      FExt->setValue(2, -m * g);
      body->setFExtPtr(FExt);

      allDS_.insert(body);

    }



    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    nsds_.reset(new NonSmoothDynamicalSystem(allDS_, allInteractions_));

    // -------------
    // --- Model ---
    // -------------
    model_.reset(new Model(t0, T));
    model_->setNonSmoothDynamicalSystemPtr(nsds_); // set NonSmoothDynamicalSystem of this model

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    timedisc_.reset(new TimeDiscretisation(t0, h));

    simulation_.reset(new TimeStepping(timedisc_));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new Moreau(allDS_, theta));

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 100;// Max number of iterations
    iparam[1] = 20; // compute error iterations
    iparam[4] = 1; // compute error iterations

    DoubleParameters dparam(5);
    dparam[0] = 1e-6; // Tolerance
    dparam[2] = 1e-8; // Local Tolerance


    SP::NonSmoothSolver mySolver;
    mySolver.reset(new NonSmoothSolver(solverName, iparam, dparam));
    osnspb_.reset(new FrictionContact(3, mySolver));

    osnspb_->setMaxSize(16384);
    osnspb_->setMStorageType(1);
    osnspb_->setNumericsVerboseMode(0);
    osnspb_->setKeepLambdaAndYState(true);
    simulation_->insertIntegrator(osi);
    simulation_->insertNonSmoothProblem(osnspb_);
    //     simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, 0.8, 3));

    playground_.reset(new SpaceFilter(3, 6, nsds_, nslaw, plans_));
    playground_->buildInteractions();

    nsds_->topology()->initialize();

    model_->initialize(simulation_);

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


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
#include <SiconosPointers.hpp>


#include "Disks.hpp"

using namespace std;

/* do nothing if solver does not converge */
void localCheckSolverOuput(int, Simulation*)
{};


// ================= Creation of the model =======================
void Disks::init()
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
    if (plans_->size(0) == 0)
    {
      /* default plans */
      double A1 = P1A;
      double B1 = P1B;
      double C1 = P1C;

      double A2 = P2A;
      double B2 = P2B;
      double C2 = P2C;

      plans_.reset(new SimpleMatrix(6, 6));
      plans_->zero();
      (*plans_)(0, 0) = 0;
      (*plans_)(0, 1) = 1;
      (*plans_)(0, 2) = -GROUND;

      (*plans_)(1, 0) = 1;
      (*plans_)(1, 1) = 0;
      (*plans_)(1, 2) = WALL;

      (*plans_)(2, 0) = 1;
      (*plans_)(2, 1) = 0;
      (*plans_)(2, 2) = -WALL;

      (*plans_)(3, 0) = 0;
      (*plans_)(3, 1) = 1;
      (*plans_)(3, 2) = -TOP;

      (*plans_)(4, 0) = A1;
      (*plans_)(4, 1) = B1;
      (*plans_)(4, 2) = C1;

      (*plans_)(5, 0) = A2;
      (*plans_)(5, 1) = B2;
      (*plans_)(5, 2) = C2;

    }

    /* set center positions */
    for (unsigned int i = 0 ; i < plans_->size(0); ++i)
    {
      SP::DiskPlanR tmpr;
      tmpr.reset(new DiskPlanR(1, (*plans_)(i, 0), (*plans_)(i, 1), (*plans_)(i, 2),
                               (*plans_)(i, 3), (*plans_)(i, 4), (*plans_)(i, 5)));
      (*plans_)(i, 3) = tmpr->getXCenter();
      (*plans_)(i, 4) = tmpr->getYCenter();
    }

    SP::SiconosMatrix Disks;
    Disks.reset(new SimpleMatrix("disks.dat", true));

    for (unsigned int i = 0; i < Disks->size(0); i++)
    {
      R = Disks->getValue(i, 2);
      m = Disks->getValue(i, 3);

      SP::SiconosVector qTmp;
      SP::SiconosVector vTmp;

      qTmp.reset(new SimpleVector(NDOF));
      vTmp.reset(new SimpleVector(NDOF));
      vTmp->zero();
      (*qTmp)(0) = (*Disks)(i, 0);
      (*qTmp)(1) = (*Disks)(i, 1);

      SP::LagrangianDS body;
      if (R > 0)
        body.reset(new Disk(R, m, qTmp, vTmp));
      else
        body.reset(new Circle(-R, m, qTmp, vTmp));

      // -- Set external forces (weight) --
      SP::SimpleVector FExt;
      FExt.reset(new SimpleVector(NDOF));
      FExt->zero();
      FExt->setValue(1, -m * g);
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

    DoubleParameters dparam(5);
    dparam[0] = 1e-3; // Tolerance

    SP::NonSmoothSolver mySolver;
    mySolver.reset(new NonSmoothSolver(solverName, iparam, dparam));
    osnspb_.reset(new FrictionContact(2, mySolver));

    osnspb_->setMaxSize(6 * ((3 * Ll * Ll + 3 * Ll) / 2 - Ll));
    osnspb_->setMStorageType(1);
    osnspb_->setNumericsVerboseMode(0);
    osnspb_->setKeepLambdaAndYState(true);
    simulation_->recordIntegrator(osi);
    simulation_->recordNonSmoothProblem(osnspb_);
    //simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, 0.8, 2));

    playground_.reset(new SpaceFilter(3, 6, nsds_, nslaw, plans_));
    playground_->buildInteractions();

    nsds_->getTopologyPtr()->initialize();

    model_->initialize(simulation_);

  }

  catch (SiconosException e)
  {
    std::cout << e.report() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cout << "Exception caught in Disks::init()" << std::endl;
    exit(1);
  }
}



// =========================== End of model definition ===========================

// ================================= Computation =================================


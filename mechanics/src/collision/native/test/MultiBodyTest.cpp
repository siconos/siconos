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
#include "MultiBodyTest.hpp"

#pragma GCC diagnostic ignored "-Wmissing-declarations"

#ifndef Disks_h
#define Disks_h

// 2D
#define NDOF 3

// WALLS, TOP and GROUND
#define WALL 100
#define TOP 100
#define GROUND 0

// DEFAULT PLANS : a ground and two walls to support a crystal

// CRYSTAL SIZE
#ifndef Ll
#define Ll 7
#endif

#define Rr 1

#define COSPI6  0.866025403784439
#define SINPI6  0.5
#define TANPI6  0.577350269189626 // tan(pi/6)

#define SY 3.73205080756888  // ((cos(a)+1)/(cos(a)*sin(a)) - tan(a)) a=pi/6, R=1
#define SYL 1/TANPI6


// Plan1
#define P1A COSPI6
#define P1B -SINPI6
#define P1C (SY+(Ll-1)*SYL-Rr)*P1B

// Plan2
#define P2A COSPI6
#define P2B SINPI6
#define P2C (SY+(Ll-1)*SYL-Rr)*P2B


#define GROUND_ID -1
#define MAX_RADIUS INFINITY

#include "SiconosBodies.hpp"
#include "Disk.hpp"
#include "Circle.hpp"
#include "DiskPlanR.hpp"
#include "SpaceFilter.hpp"

class Disks : public SiconosBodies, public std11::enable_shared_from_this<Disks>
{
public:
  void init()
  {
    assert(false);
  };
  void init(std::string);
};

TYPEDEF_SPTR(Disks)

#endif //Disks_h

// Siconos
#include <SiconosKernel.hpp>
#include <SiconosPointers.hpp>

using namespace std;

/* do nothing if solver does not converge */
void localCheckSolverOuput(int info, Simulation*)
{
  if (info) exit(1);
}

double A(double t)
{
  return 0. ;
}
double B(double t)
{
  return 1. ;
}
double C(double t)
{
  return 0.0;//1.1*cos(32.*M_PI*t) ;
}
double DA(double t)
{
  return 0. ;
}
double DB(double t)
{
  return 0. ;
}
double DC(double t)
{
  return 0.0;//-1.1*32.*M_PI*sin(32.*M_PI*t) ;
}


// ================= Creation of the model =======================
void Disks::init(std::string disks_input)
{

  SP::TimeDiscretisation timedisc_;
  SP::TimeStepping simulation_;
  SP::FrictionContact osnspb_;

  // User-defined main parameters

  double t0 = 0;                   // initial computation time

  double T = 10;

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
    if (_plans->size(0) == 0)
    {
      /* default plans */
      double A1 = P1A;
      double B1 = P1B;
      double C1 = P1C;

      double A2 = P2A;
      double B2 = P2B;
      double C2 = P2C;

      _plans.reset(new SimpleMatrix(6, 6));
      _plans->zero();
      (*_plans)(0, 0) = 0;
      (*_plans)(0, 1) = 1;
      (*_plans)(0, 2) = -GROUND;

      (*_plans)(1, 0) = 1;
      (*_plans)(1, 1) = 0;
      (*_plans)(1, 2) = WALL;

      (*_plans)(2, 0) = 1;
      (*_plans)(2, 1) = 0;
      (*_plans)(2, 2) = -WALL;

      (*_plans)(3, 0) = 0;
      (*_plans)(3, 1) = 1;
      (*_plans)(3, 2) = -TOP;

      (*_plans)(4, 0) = A1;
      (*_plans)(4, 1) = B1;
      (*_plans)(4, 2) = C1;

      (*_plans)(5, 0) = A2;
      (*_plans)(5, 1) = B2;
      (*_plans)(5, 2) = C2;

    }

    /* set center positions */
    for (unsigned int i = 0 ; i < _plans->size(0); ++i)
    {
      SP::DiskPlanR tmpr;
      tmpr.reset(new DiskPlanR(1, (*_plans)(i, 0), (*_plans)(i, 1), (*_plans)(i, 2),
                               (*_plans)(i, 3), (*_plans)(i, 4), (*_plans)(i, 5)));
      (*_plans)(i, 3) = tmpr->getXCenter();
      (*_plans)(i, 4) = tmpr->getYCenter();
    }

    /*    _moving_plans.reset(new FMatrix(1,6));
        (*_moving_plans)(0,0) = &A;
        (*_moving_plans)(0,1) = &B;
        (*_moving_plans)(0,2) = &C;
        (*_moving_plans)(0,3) = &DA;
        (*_moving_plans)(0,4) = &DB;
        (*_moving_plans)(0,5) = &DC;*/



    SP::SiconosMatrix Disks;
    Disks.reset(new SimpleMatrix(disks_input, true));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new MoreauJeanOSI(theta));

    // -- Model --
    _model.reset(new Model(t0, T));

    for (unsigned int i = 0; i < Disks->size(0); i++)
    {
      R = Disks->getValue(i, 2);
      m = Disks->getValue(i, 3);

      SP::SiconosVector qTmp;
      SP::SiconosVector vTmp;

      qTmp.reset(new SiconosVector(NDOF));
      vTmp.reset(new SiconosVector(NDOF));
      vTmp->zero();
      (*qTmp)(0) = (*Disks)(i, 0);
      (*qTmp)(1) = (*Disks)(i, 1);

      SP::LagrangianDS body;
      if (R > 0)
        body.reset(new Disk(R, m, qTmp, vTmp));
      else
        body.reset(new Circle(-R, m, qTmp, vTmp));

      // -- Set external forces (weight) --
      SP::SiconosVector FExt;
      FExt.reset(new SiconosVector(NDOF));
      FExt->zero();
      FExt->setValue(1, -m * g);
      body->setFExtPtr(FExt);

      // add the dynamical system in the non smooth dynamical system
      _model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);
    }


    _model->nonSmoothDynamicalSystem()->setSymmetric(true);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    timedisc_.reset(new TimeDiscretisation(t0, h));

    // -- OneStepNsProblem --
    osnspb_.reset(new FrictionContact(2));

    osnspb_->numericsSolverOptions()->iparam[0] = 100; // Max number of
    // iterations
    osnspb_->numericsSolverOptions()->iparam[1] = 20; // compute error
    // iterations
    osnspb_->numericsSolverOptions()->dparam[0] = 1e-3; // Tolerance


    osnspb_->setMaxSize(6 * ((3 * Ll * Ll + 3 * Ll) / 2 - Ll));
    osnspb_->setMStorageType(1);            // Sparse storage
    osnspb_->setNumericsVerboseMode(0);

    osnspb_->setKeepLambdaAndYState(true);  // inject previous solution

    // -- Simulation --
    simulation_.reset(new TimeStepping(timedisc_));

    std11::static_pointer_cast<TimeStepping>(simulation_)->setNewtonMaxIteration(3);

    simulation_->insertIntegrator(osi);
    simulation_->insertNonSmoothProblem(osnspb_);

    simulation_->setCheckSolverFunction(localCheckSolverOuput);

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0, 0, 0.3, 2));

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
    std::cout << "Exception caught in Disks::init()" << std::endl;
    exit(1);
  }
}


// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(MultiBodyTest);


void MultiBodyTest::setUp()
{
}

void MultiBodyTest::tearDown()
{
}


// multiples disks
void MultiBodyTest::t1()
{
  SP::Disks disks(new Disks());

  disks->init("disks.dat");


  // just try to run a simulation
  // if something is broken with SpaceFilter
  // an exception may occurs
  for (unsigned int i = 0; i < 20; ++i)
  {
    disks->compute();
  }

  CPPUNIT_ASSERT(1);

}

// one disk without interaction at the beginning
void MultiBodyTest::t2()
{
  SP::Disks disks(new Disks());

  disks->init("disks-nointer.dat");


  // just try to run a simulation
  // if something is broken with SpaceFilter
  // an exception may occurs
  // test fail with rev 3146
  for (unsigned int i = 0; i < 20; ++i)
  {
    disks->compute();
  }

  CPPUNIT_ASSERT(1);


}

void MultiBodyTest::t3()
{
}

void MultiBodyTest::t4()
{
}

void MultiBodyTest::t5()
{
}

void MultiBodyTest::t6()
{
}


void MultiBodyTest::t7()
{
}

void MultiBodyTest::t8()
{
}

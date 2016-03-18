
#include "testContact.hpp"

#include "Contactor.hpp"
#include "SiconosShape.hpp"
#include "BulletBroadphase.hpp"
#include "BodyDS.hpp"
#include "MechanicsTimeStepping.hpp"

#include "SiconosKernel.hpp"

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ContactTest);

void ContactTest::setUp() {}
void ContactTest::tearDown() {}

void ContactTest::t1()
{
  try
  {
    // Initial state
    SP::SiconosVector pos(new SiconosVector(7));
    SP::SiconosVector vel(new SiconosVector(6));
    pos->zero();
    vel->zero();

    // Set up a Siconos Mechanics environment:
    // A BodyDS with a contactor consisting of a single sphere.
    (*pos)(2) = 10.0;
    SP::BodyDS body(new BodyDS(pos, vel, 1.0));
    SP::Contactor contactor(new Contactor());
    SP::SiconosSphere sphere(new SiconosSphere(0,0,0,1.0));
    contactor->addShape(sphere);
    body->setContactor(contactor);

    // A BodyDS with a contactor consisting of a plane
    (*pos)(2) = 0.0;
    SP::BodyDS body2(new BodyDS(pos, vel, 1.0));
    SP::Contactor contactor2(new Contactor());
    SP::SiconosPlane plane(new SiconosPlane(0,0,1,0.0));
    contactor2->addShape(plane);
    body2->setContactor(contactor2);

    // Object to manage the Bullet implementation of broadphase
    SP::BulletBroadphase broadphase(new BulletBroadphase());

    // Build Bullet representational mirror of contactors
    broadphase->buildGraph(contactor);
    broadphase->buildGraph(contactor2);

    // Perform Bullet broadphase, generates IndexSet1
    broadphase->performBroadphase();

    // Update a property
    sphere->setRadius(0.5);

    // Check for dirty objects and update the graph
    broadphase->updateGraph();

    // Perform Bullet broadphase, generates IndexSet1
    broadphase->performBroadphase();
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

void ContactTest::t2()
{
  try
  {
    // User-defined main parameters
    double t0 = 0;                   // initial computation time
    double T = 20.0;                 // end of computation time
    double h = 0.005;                // time step
    double position_init = 10.0;     // initial position
    double velocity_init = 0.0;      // initial velocity

    double g = 9.81;
    double theta = 0.5;              // theta for MoreauJeanOSI integrator

    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new MoreauJeanOSI(theta));

    // -- Model --
    SP::Model model(new Model(t0, T));

    SP::SiconosVector q0(new SiconosVector(7));
    SP::SiconosVector v0(new SiconosVector(6));
    q0->zero();
    v0->zero();

    (*q0)(2) = position_init;
    (*q0)(3) = 1.0;
    (*v0)(2) = velocity_init;

    SP::SiconosVector q1(new SiconosVector(7));
    SP::SiconosVector v1(new SiconosVector(6));
    q1->zero();
    v1->zero();

    /////////

    // Bodies

    // Set up a Siconos Mechanics environment:
    // A BodyDS with a contactor consisting of a single sphere.
    SP::BodyDS body(new BodyDS(q0, v0, 1.0));
    SP::Contactor contactor(new Contactor());
    SP::SiconosSphere sphere(new SiconosSphere(0,0,10,1.0));
    contactor->addShape(sphere);
    body->setContactor(contactor);

    // A BodyDS with a contactor consisting of a plane
    SP::BodyDS body2(new BodyDS(q0, v0, 1.0));
    SP::Contactor contactor2(new Contactor());
    SP::SiconosPlane plane(new SiconosPlane(0,0,1,0.0));
    contactor2->addShape(plane);
    body2->setContactor(contactor2);

    /////////

    // -- Set external forces (weight) --
    float mass = 1.0;
    SP::SiconosVector FExt;
    FExt.reset(new SiconosVector(3)); //
    FExt->zero();
    FExt->setValue(2, - g * mass);
    body->setFExtPtr(FExt);

    // -- Add the dynamical system in the non smooth dynamical system
    osi->insertDynamicalSystem(body);
    model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);

    // -- Time discretisation --
    SP::TimeDiscretisation timedisc(new TimeDiscretisation(t0, h));

    // -- OneStepNsProblem --
    SP::FrictionContact osnspb(new FrictionContact(3));

    // -- Some configuration
    osnspb->numericsSolverOptions()->iparam[0] = 1000;
    osnspb->numericsSolverOptions()->dparam[0] = 1e-5;

    osnspb->setMaxSize(16384);
    osnspb->setMStorageType(1);
    osnspb->setNumericsVerboseMode(0);
    osnspb->setKeepLambdaAndYState(true);

    /////////

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    int N = ceil((T - t0) / h); // Number of time steps

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0.8, 0., 0.0, 3));

    // TODO pass nslaw to broadphase

    // -- MoreauJeanOSI Time Stepping with Bullet Dynamical Systems
    SP::MechanicsTimeStepping simulation(new MechanicsTimeStepping(timedisc));

    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

    model->initialize(simulation);

    // Object to manage the Bullet implementation of broadphase
    SP::BulletBroadphase broadphase(new BulletBroadphase());

    // Build Bullet representational mirror of contactors
    broadphase->buildGraph(model);

    std::cout << "====> End of initialisation ..." << std::endl << std::endl;

    ///////

    int k=0;
    while (simulation->hasNextEvent())
    {
      // Update a property at step 100
      if (k==100) {
        sphere->setRadius(0.5);
      }

      // Check for dirty objects and update the broadphase graph
      broadphase->updateGraph();

      // Perform Bullet broadphase, generates IndexSet1
      broadphase->performBroadphase();

      // Update integrator and solve constraints
      simulation->computeOneStep();

      // Advance simulation
      simulation->nextStep();
      k++;
    }

    std::cout << std::endl << "End of computation - Number of iterations done: " << k - 1 << std::endl;
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

// -*- compile-command: "make -C ~/projects/siconos/bld/mechanics && valgrind --leak-check=full --suppressions=$HOME/projects/siconos/cmake/valgrind.supp ~/projects/siconos/bld/mechanics/src/proposed/test/testContact ContactTest" -*-

#include "ContactTest.hpp"

#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"
#include "BulletBroadphase.hpp"
#include "BodyDS.hpp"
#include "BodyTimeStepping.hpp"

#include "SiconosKernel.hpp"

#include <string>

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ContactTest);

void ContactTest::setUp() {}
void ContactTest::tearDown() {}

void ContactTest::t1()
{
  try
  {
    printf("\n==== t1\n");

    // Initial state
    SP::SiconosVector pos(new SiconosVector(7));
    SP::SiconosVector vel(new SiconosVector(6));
    pos->zero();
    vel->zero();

    // Set up a Siconos Mechanics environment:
    // A BodyDS with a contactor consisting of a single sphere.
    (*pos)(2) = 10.0;
    SP::BodyDS body(new BodyDS(pos, vel, 1.0));
    SP::SiconosContactor contactor(new SiconosContactor());
    SP::SiconosSphere sphere(new SiconosSphere(0,0,0,1.0));
    contactor->addShape(sphere);
    body->setContactor(contactor);

    // A BodyDS with a contactor consisting of a plane
    (*pos)(2) = 0.0;
    SP::BodyDS body2(new BodyDS(pos, vel, 1.0));
    SP::SiconosContactor contactor2(new SiconosContactor());
    SP::SiconosPlane plane(new SiconosPlane(0,0,0));
    contactor2->addShape(plane);
    body2->setContactor(contactor2);

    // Object to manage the Bullet implementation of broadphase
    SP::SiconosBroadphase broadphase(new BulletBroadphase());

    // Build broadphase-specific mirror of contactor graph
    std::vector<SP::BodyDS> bodies;
    bodies.push_back(body);
    bodies.push_back(body2);
    broadphase->buildGraph(bodies);

    // Perform broadphase, generates IndexSet0
    broadphase->performBroadphase();

    // Update a property
    sphere->setRadius(0.5);

    // Check for dirty objects and update the graph
    broadphase->updateGraph();

    // Perform broadphase, generates IndexSet0
    broadphase->performBroadphase();
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

struct BounceParams
{
  bool trace;
  bool dynamic;
  double size;
  double mass;
  double margin;

  void dump() {
    printf("  trace:    %s\n", trace?"on":"off");
    printf("  dynamic:  %s\n", trace?"on":"off");
    printf("  size:     %.1g\n", size);
    printf("  mass:     %.1g\n", mass);
    printf("  margin:   %.1g\n", margin);
  }
};

struct BounceResult
{
  double final_position;
  double final_position_std;
};

BounceResult bounceTest(std::string moving,
                        std::string ground,
                        const BounceParams &params)
{
    // User-defined main parameters
    double t0 = 0;                   // initial computation time
    double T = 20.0;                 // end of computation time
    double h = 0.005;                // time step
    double position_init = 3.0;      // initial position
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
    (*q1)(3) = 1.0;
    v1->zero();

    /////////

    // Bodies

    printf("== Testing: %s falling on %s .. ",
           moving.c_str(), ground.c_str());
    if (params.trace) printf("==\n");
    fflush(stdout);

    // TODO: set the shape margin

    // Set up a Siconos Mechanics environment:
    // A BodyDS with a contactor consisting of a single sphere.
    SP::BodyDS body(new BodyDS(q0, v0, params.mass));
    SP::SiconosContactor contactor(new SiconosContactor());
    SP::SiconosSphere sphere;
    if (moving=="sphere")
    {
      sphere.reset(new SiconosSphere(0,0,0, params.size/2));
      contactor->addShape(sphere);
    }
    else if (moving=="box")
    {
      SP::SiconosBox box(
        new SiconosBox(0,0,0,params.size,params.size,params.size));
      contactor->addShape(box);
    }
    body->setContactor(contactor);

    // A contactor with no body (static contactor) consisting of a plane
    SP::SiconosContactor static_contactor(new SiconosContactor());
    if (ground=="plane")
    {
      SP::SiconosPlane plane(new SiconosPlane(0,0,0));
      static_contactor->addShape(plane);
    }
    else if (ground=="box")
    {
      SP::SiconosBox floorbox(new SiconosBox(0,0,-50,10,10,100));
      static_contactor->addShape(floorbox);
    }
    else if (ground=="sphere")
    {
      SP::SiconosSphere floorsphere(new SiconosSphere(0,0,-1.0,1.0));
      static_contactor->addShape(floorsphere);
    }

    /////////

    // -- Set external forces (weight) --
    SP::SiconosVector FExt;
    FExt.reset(new SiconosVector(3));
    FExt->zero();
    FExt->setValue(2, - g * params.mass);
    body->setFExtPtr(FExt);

    // -- Add the dynamical systems into the non smooth dynamical system
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

    int N = ceil((T - t0) / h); // Number of time steps

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0.8, 0., 0.0, 3));

    // TODO pass nslaw to broadphase

    // -- MoreauJeanOSI Time Stepping with Body-based Dynamical Systems
    SP::BodyTimeStepping simulation(new BodyTimeStepping(timedisc));

    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

    model->initialize(simulation);

    // Object to manage the Bullet implementation of broadphase
    SP::SiconosBroadphase broadphase(new BulletBroadphase());

    // Build broadphase-specific mirror of contactor graph
    broadphase->buildGraph(model);
    broadphase->buildGraph(static_contactor);

    ///////

    int k=0;
    double std=0, final_pos=0;
    while (simulation->hasNextEvent())
    {
      // Update a property at step 500
      if (params.dynamic && k==500 && moving=="sphere") {
        sphere->setRadius(0.3);
      }

      // Check for dirty objects and update the broadphase graph
      broadphase->updateGraph();

      // Perform broadphase, generates IndexSet0
      broadphase->performBroadphase();

      // Update integrator and solve constraints
      simulation->computeOneStep();

      if (params.trace && k<3999) {
        printf("pos, %f, (%f, %f, %f, %f)\n", (*body->q())(2),
               (*body->q())(3), (*body->q())(4),
               (*body->q())(5), (*body->q())(6));
      }

      double pos = (*body->q())(2);

      // Standard deviation (cheating by not calculating mean!)
      if (k==(3999-100))
        final_pos = pos;
      if (k>(3999-100)) {
        std += (pos-final_pos)*(pos-final_pos);
      }

      // Advance simulation
      simulation->nextStep();
      k++;
    }

    if (!params.trace)
      printf("(done, iterations=%d)\n", k - 1);

    BounceResult r;
    r.final_position =  (*body->q())(2);
    r.final_position_std = sqrt(std/100);
    return r;
}

void ContactTest::t2()
{
  try
  {
    printf("\n==== t2\n");

    BounceParams p;
    p.trace = true;
    p.dynamic = false;
    p.size = 0.1;
    p.mass = 1.0;
    p.margin = 0.1; // TODO

    BounceResult r = bounceTest("box", "plane", p);

    fprintf(stderr, "\nFinal position: %g  (std=%g)\n\n",
            r.final_position, r.final_position_std);
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

void ContactTest::t3()
{
  try
  {
    printf("\n==== t3\n");

    BounceParams params;
    params.trace = false;
    params.dynamic = false;
    params.size = 0.1;
    params.mass = 1.0;
    params.margin = 0.1; // TODO

    BounceResult results[2][3];
    results[0][0] = bounceTest("sphere", "plane",  params);
    results[1][0] = bounceTest("box",    "plane",  params);
    results[0][1] = bounceTest("sphere", "sphere", params);
    results[1][1] = bounceTest("box",    "sphere", params);
    results[0][2] = bounceTest("sphere", "box",    params);
    results[1][2] = bounceTest("box",    "box",    params);

    // Report
    printf("\nParams:\n\n");
    params.dump();
    printf("\nFinal resting positions:\n\n");
    printf("       | sphere | box    | plane\n");
    printf("-------+--------+--------+-------\n");
    printf("sphere | %#6.3f | %#6.3f | %#6.3f\n",
           results[0][0].final_position,
           results[0][1].final_position,
           results[0][2].final_position);
    printf("box    | %#6.3f | %#6.3f | %#6.3f\n",
           results[1][0].final_position,
           results[1][1].final_position,
           results[1][2].final_position);
    printf("\nFinal resting position standard deviation:\n\n");
    printf("       | sphere   | box      | plane\n");
    printf("-------+----------+----------+---------\n");
    printf("sphere | %#.2e | %#.2e | %#.2e\n",
           results[0][0].final_position_std,
           results[0][1].final_position_std,
           results[0][2].final_position_std);
    printf("box    | %#.2e | %#.2e | %#.2e\n",
           results[1][0].final_position_std,
           results[1][1].final_position_std,
           results[1][2].final_position_std);
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

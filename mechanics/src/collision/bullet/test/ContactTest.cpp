// -*- compile-command: "make -C ~/projects/siconos/bld/mechanics && valgrind --leak-check=full --suppressions=$HOME/projects/siconos/cmake/valgrind.supp ~/projects/siconos/bld/mechanics/src/collision/bullet/test/testContact ContactTest" -*-
// make --no-print-directory -C ~/projects/siconos/bld/mechanics && ~/projects/siconos/bld/mechanics/src/collision/bullet/test/testContact ContactTest::t2 | grep pos, | cut -d, -f2-8 | plot.py -s

#include "ContactTest.hpp"

#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"
#include "SiconosCollisionManager.hpp"
#include "SiconosBulletCollisionManager.hpp"
#include "BodyDS.hpp"

#include "SiconosKernel.hpp"

#include <string>
#include <sys/time.h>
#include <boost/make_shared.hpp>

// Experimental settings for SiconosBulletCollisionManager
extern double extra_margin;
extern double breaking_threshold;
extern double box_ch_added_dimension;
extern double box_convex_hull_margin;
extern double bullet_world_scaling;

// Experimental statistics from SiconosBulletCollisionManager
extern int new_interaction_counter;
extern int existing_interaction_counter;
extern int interaction_warning_counter;

// test suite registration
CPPUNIT_TEST_SUITE_REGISTRATION(ContactTest);

void ContactTest::setUp() {}
void ContactTest::tearDown() {}

struct BounceParams
{
  bool trace;
  bool dynamic;
  double size;
  double mass;
  double position;
  double timestep;
  double insideMargin;
  double outsideMargin;
  SiconosBulletOptions options;

  void dump() {
    printf("  trace:              %s\n", trace?"on":"off");
    printf("  dynamic:            %s\n", trace?"on":"off");
    printf("  size:               %.3g\n", size);
    printf("  mass:               %.3g\n", mass);
    printf("  position:           %.3g\n", position);
    printf("  insideMargin:       %.3g\n", insideMargin);
    printf("  outsideMargin:      %.3g\n", outsideMargin);
    printf("  breakingThreshold:  %.3g\n", options.breakingThreshold);
    printf("  worldScale:         %.3g\n", options.worldScale);
  }
};

struct BounceResult
{
  double bounce_error_sum;
  double bounce_error[6];
  int n_bounce_error;
  double final_position;
  double final_position_std;
  int num_interactions;
  int num_interaction_warnings;
  int max_simultaneous_contacts;
  double avg_simultaneous_contacts;
  double displacement_on_first_contact;
};

BounceResult bounceTest(std::string moving,
                        std::string ground,
                        const BounceParams &params)
{
    // User-defined main parameters
    double t0 = 0;                   // initial computation time
    double T = 20.0;                 // end of computation time
    double h = params.timestep;      // time step
    double position_init = params.position;  // initial position
    double velocity_init = 0.0;      // initial velocity

    double g = 9.81;
    double theta = 0.5;              // theta for MoreauJeanOSI integrator

    int steps = (T-t0)/h;

    // Statistics of this run
    int bounce_counter = 0;
    const int n_desired_bounces = 6;
    double desired_bounce_times[6]  = { 2.585, 4.645, 6.290, 7.610, 8.660, 9.505 };
    double desired_bounce_ratios[6] = { 0.6358, 0.4048, 0.2577, 0.1630, 0.1033, 0.0647 };
    double actual_bounce_times[6]   = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
    double actual_bounce_ratios[6]  = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

    int local_new_interaction_count=0;
    int max_simultaneous_contacts=0;
    double avg_simultaneous_contacts=0.0;

    bool first_contact = true;
    double displacement_on_first_contact=0.0;

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
    SP::SiconosContactorSet contactors(new SiconosContactorSet());
    SP::SiconosSphere sphere;
    if (moving=="sphere")
    {
      sphere.reset(new SiconosSphere(params.size/2));
      sphere->setInsideMargin(params.insideMargin);
      sphere->setOutsideMargin(params.outsideMargin);
      contactors->push_back(std11::make_shared<SiconosContactor>(sphere));
    }
    else if (moving=="box")
    {
      SP::SiconosBox box(new SiconosBox(params.size,params.size,params.size));
      box->setInsideMargin(params.insideMargin);
      box->setOutsideMargin(params.outsideMargin);
      contactors->push_back(std11::make_shared<SiconosContactor>(box));
    }
    else if (moving=="ch")
    {
      float siz = params.size;
      SP::SiconosMatrix pts(new SimpleMatrix(4,3));
      (*pts)(0,0) = 0.0; (*pts)(0,1) = 0.0; (*pts)(0,2) = 0.0;
      (*pts)(1,1) = siz; (*pts)(1,1) = 0.0; (*pts)(1,2) = 0.0;
      (*pts)(2,0) = 0.0; (*pts)(2,1) = siz; (*pts)(2,2) = 0.0;
      (*pts)(3,0) = 0.0; (*pts)(3,1) = 0.0; (*pts)(3,2) = siz;
      SP::SiconosConvexHull ch(new SiconosConvexHull(pts));
      ch->setInsideMargin(params.insideMargin);
      ch->setOutsideMargin(params.outsideMargin);
      contactors->push_back(std11::make_shared<SiconosContactor>(ch));
    }
    body->contactors() = contactors;

    // A contactor with no body (static contactor) consisting of a plane
    // Positioned such that bouncing and resting position = 0.0
    SP::SiconosContactorSet static_contactors;
    if (ground=="plane")
    {
      SP::SiconosPlane plane(new SiconosPlane());
      plane->setInsideMargin(params.insideMargin);
      plane->setOutsideMargin(params.outsideMargin);
      static_contactors->push_back(std11::make_shared<SiconosContactor>(plane));
    }
    else if (ground=="box")
    {
      SP::SiconosBox floorbox(new SiconosBox(100,100,100));
      floorbox->setInsideMargin(params.insideMargin);
      floorbox->setOutsideMargin(params.outsideMargin);
      SP::SiconosVector pos(new SiconosVector(7));
      (*pos)(2) = -50-params.size/2;
      (*pos)(3) = 1.0;
      static_contactors->push_back(std11::make_shared<SiconosContactor>(floorbox, pos));
    }
    else if (ground=="sphere")
    {
      SP::SiconosSphere floorsphere(new SiconosSphere(1.0));
      floorsphere->setInsideMargin(params.insideMargin);
      floorsphere->setOutsideMargin(params.outsideMargin);
      SP::SiconosVector pos(new SiconosVector(7));
      (*pos)(2) = -1-params.size/2;
      (*pos)(3) = 1.0;
      static_contactors->push_back(std11::make_shared<SiconosContactor>(floorsphere, pos));
    }

    /////////

    // -- Set external forces (weight) --
    SP::SiconosVector FExt;
    FExt.reset(new SiconosVector(3));
    FExt->zero();
    FExt->setValue(2, - g * params.mass);
    body->setFExtPtr(FExt);

    // -- Add the dynamical systems into the non smooth dynamical system
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

    // TODO pass nslaw to collisionMan

    // -- MoreauJeanOSI Time Stepping
    SP::TimeStepping simulation(new TimeStepping(timedisc));

    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

    model->setSimulation(simulation);
    model->initialize();

    // Object to manage the Bullet implementation of collisionMan
    SP::SiconosBulletCollisionManager collisionMan(
      new SiconosBulletCollisionManager(params.options));

    simulation->insertInteractionManager(collisionMan);

    // Add static shapes (centered at zero by default)
    collisionMan->insertStaticContactorSet(static_contactors);

    ///////

    int new_interaction_total = 0;
    collisionMan->resetStatistics();

    ///////

    int k=0;
    double std=0, final_pos=0;
    double last_pos=position_init, last_vel=0;
    while (simulation->hasNextEvent())
    {
      // Update a property at step 500
      if (params.dynamic && k==500 && moving=="sphere") {
        sphere->setRadius(0.3);
      }

      // Update integrator and solve constraints
      simulation->computeOneStep();

      double vel = (*body->velocity())(2);
      double pos = (*body->q())(2);

      if (params.trace && (k+1) < steps) {
        printf("pos, %f, %f, %f\n", pos, last_pos - pos, vel);
      }

      // Peaks (velocity crosses zero towards positive)
      if (vel <= 0 && last_vel > 0) {
        if (bounce_counter < n_desired_bounces)
        {
          actual_bounce_times[bounce_counter] = k*h;
          actual_bounce_ratios[bounce_counter] = pos / position_init;
          bounce_counter ++;
        }
      }

      // Interaction statistics
      if (collisionMan->statistics().new_interactions_created > 0 && first_contact) {
        first_contact = false;
        displacement_on_first_contact = last_pos - pos;
      }

      int interactions = collisionMan->statistics().new_interactions_created
        + collisionMan->statistics().existing_interactions_processed;

      local_new_interaction_count += collisionMan->statistics().new_interactions_created;

      if (interactions > max_simultaneous_contacts)
        max_simultaneous_contacts = interactions;

      avg_simultaneous_contacts += interactions;

      // Reset interaction counters for next iteration
      collisionMan->resetStatistics();

      // Standard deviation (cheating by not calculating mean!)
      if (k==(steps-100))
        final_pos = pos;
      if (k>(steps-100)) {
        std += (pos-final_pos)*(pos-final_pos);
      }

      // State
      last_pos = pos;
      last_vel = vel;

      // Advance simulation
      simulation->nextStep();
      k++;
    }

    if (!params.trace)
      printf("(done, iterations=%d)\n", k - 1);

    BounceResult r;
    r.bounce_error_sum = 0.0;
    r.n_bounce_error = 6;
    for (int i=0; i < r.n_bounce_error; i++) {
      double er = actual_bounce_ratios[i] - desired_bounce_ratios[i];
      r.bounce_error[i] = fabs(er);
      r.bounce_error_sum += er*er;
    }
    r.bounce_error_sum = sqrt(r.bounce_error_sum/r.n_bounce_error);
    r.final_position = (*body->q())(2);
    r.final_position_std = sqrt(std/100);

    r.num_interactions = local_new_interaction_count;
    r.num_interaction_warnings = collisionMan->statistics().interaction_warnings;
    r.max_simultaneous_contacts = max_simultaneous_contacts;
    r.avg_simultaneous_contacts = avg_simultaneous_contacts / (double)k;

    r.displacement_on_first_contact = displacement_on_first_contact;

    return r;
}

void ContactTest::t1()
{
  try
  {
    printf("\n==== t1\n");

    BounceParams params;
    params.trace = true;
    params.dynamic = false;
    params.size = 1.0;
    params.mass = 1.0;
    params.position = 3.0;
    params.timestep = 0.005;
    params.insideMargin = 0.3;
    params.outsideMargin = 0.3;

    BounceResult r = bounceTest("box", "box", params);

    fprintf(stderr, "\nSize: %g\n", params.size);
    fprintf(stderr, "Final position: %g  (std=%g)\n\n",
            r.final_position, r.final_position_std);
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
    printf("\n==== t2\n");

    BounceParams params;
    params.trace = false;
    params.dynamic = false;
    params.size = 0.01;
    params.mass = 1.0;
    params.position = 3.0;
    params.timestep = 0.005;
    params.insideMargin = 0.1;
    params.outsideMargin = 0.1;

    BounceResult results[2][3];
    results[0][0] = bounceTest("sphere", "sphere",  params);
    results[1][0] = bounceTest("box",    "sphere",  params);
    results[0][1] = bounceTest("sphere", "box", params);
    results[1][1] = bounceTest("box",    "box", params);
    results[0][2] = bounceTest("sphere", "plane",    params);
    results[1][2] = bounceTest("box",    "plane",    params);

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

void ContactTest::t3()
{
  try
  {
    printf("\n==== t3\n");

    BounceParams params[3];
    params[0].trace = true;
    params[0].dynamic = false;
    params[0].size = 0.02;
    params[0].mass = 0.02;
    params[0].position = 10.0;
    params[0].timestep = 0.005;
    params[0].insideMargin = 0.1;
    params[0].outsideMargin = 0.1;

    params[1] = params[0];
    params[1].size = 0.1;
    params[1].mass = 0.1;
    params[1].position = 10.0;

    params[2] = params[0];
    params[2].size = 1.0;
    params[2].mass = 1.0;
    params[2].position = 10.0;

    BounceResult results[3];
    int i;
    for (i=0; i<3; i++) {
      results[i] = bounceTest("box", "box", params[i]);
    }

    // Report
    fprintf(stderr, "\nFinal resting positions:\n\n");
    for (i=0; i<3; i++) {
      fprintf(stderr, "   pos %.2f: %f  (std=%f)\n",
              params[i].position,
              results[i].final_position,
              results[i].final_position_std);
    }
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

template<typename T>
class Var
{
public:
  T valmin, valmax;
  double prob;

  Var(T _valmin, T _valmax)
    : valmin(_valmin), valmax(_valmax) {}
  Var(double _prob) : prob(_prob) {}
  T sample() {
    return (T)(rand()*(double)(valmax-valmin)/(double)RAND_MAX) + valmin;
  }
};

template<>
bool Var<bool>::sample() { return (rand()/(double)RAND_MAX) > prob; }

void ContactTest::t4()
{
  printf("\n==== t4\n");

  FILE *fresults = fopen("results.json", "w");
  if (!fresults) {
    CPPUNIT_ASSERT(1);
  }

  fprintf(fresults, "[\n");

  int n_tests = 3; // Set this to a large number to gather proper statistics

  float t0 = clock();
  int test_count=0;

  Var<double> var_inside_margin(0, 0.5);
  Var<double> var_outside_margin(0, 0.5);
  Var<double> var_breaking_threshold(0, 1.0);
  Var<double> var_params_size(log(0.01), log(10.0));
  Var<double> var_params_mass(log(0.01), log(10.0));
  Var<double> var_params_position(3.0, 15.0);
  Var<int> var_world_scale(log(0.001), log(100.0));
  Var<int> var_params_timestep(0, 3);
  double possible_timesteps[4] = { 0.0005, 0.001, 0.005, 0.01 };

  while (test_count < n_tests)
  {
    long elapsed = (long)((clock() - t0) / (double)CLOCKS_PER_SEC);
    long remaining = (long)(elapsed / (test_count+1) / 60.0 * (n_tests - test_count));

    printf("Iteration %d (%0.2f%%, remaining %ld min): ",
           test_count,
           test_count*100/(double)n_tests,
           remaining);

    fprintf(fresults, "%s{\n", test_count==0 ? "" : ",\n");
    test_count ++;

    BounceParams params;
    params.trace = false;
    params.dynamic = false;
    params.size = exp(var_params_size.sample());
    params.mass = exp(var_params_mass.sample());
    params.position = var_params_position.sample();
    params.timestep = possible_timesteps[var_params_timestep.sample()];
    params.insideMargin = var_inside_margin.sample();
    params.outsideMargin = var_outside_margin.sample();

    // Experimental settings
    SiconosBulletOptions options;
    options.breakingThreshold = var_breaking_threshold.sample();
    options.worldScale = exp(var_world_scale.sample());
    params.options = options;

    bool success = false;

    fprintf(fresults, "  \"size\": %g,\n", params.size);
    fprintf(fresults, "  \"mass\": %g,\n", params.mass);
    fprintf(fresults, "  \"position\": %g,\n", params.position);
    fprintf(fresults, "  \"timestep\": %g,\n", params.timestep);
    fprintf(fresults, "  \"inside_margin\": %g,\n", params.insideMargin);
    fprintf(fresults, "  \"outside_margin\": %g,\n", params.outsideMargin);
    fprintf(fresults, "  \"breaking_threshold\": %g,\n", options.breakingThreshold);
    fprintf(fresults, "  \"world_scale\": %g,\n", options.worldScale);

    BounceResult results;
    try
    {
      results = bounceTest("box", "plane", params);
      success = true;
    }
    catch (SiconosException e)
    {
      std::cout << "SiconosException: " << e.report() << std::endl;
    }

    fprintf(fresults, "  \"success\": %s,\n", success ? "true" : "false");
    fprintf(fresults, "  \"bounce_error_sum\": %g,\n", results.bounce_error_sum);
    fprintf(fresults, "  \"bounce_error\": [");
    for (int i=0; i < results.n_bounce_error; i++)
      fprintf(fresults, "%s %g", i==0 ? "" : ",", results.bounce_error[i]);
    fprintf(fresults, " ],\n");
    fprintf(fresults, "  \"final_position\": %g,\n", results.final_position);
    fprintf(fresults, "  \"final_position_std\": %g,\n", results.final_position_std);

    fprintf(fresults, "  \"num_interactions\": %d,\n", results.num_interactions);
    fprintf(fresults, "  \"num_interaction_warnings\": %d,\n", results.num_interaction_warnings);
    fprintf(fresults, "  \"max_simultaneous_contacts\": %d,\n", results.max_simultaneous_contacts);
    fprintf(fresults, "  \"avg_simultaneous_contacts\": %g,\n", results.avg_simultaneous_contacts);
    fprintf(fresults, "  \"displacement_on_first_contact\": %g\n", results.displacement_on_first_contact);

    fprintf(fresults, "}");
    fflush(fresults);
  }

  fprintf(fresults, "\n]\n");
  fclose(fresults);

  CPPUNIT_ASSERT(1);
}

void ContactTest::t5()
{
  printf("\n==== t5\n");

  BounceParams params;
  params.trace = false;
  params.dynamic = false;
  params.size = 1.0;
  params.mass = 1.0;
  params.position = 3.0;
  params.timestep = 0.005;
  params.insideMargin = 0.1;
  params.outsideMargin = 0.1;

  SiconosBulletOptions options;
  options.breakingThreshold = 0.4;
  options.worldScale = 1.0;
  params.options = options;

  bool success = false;

  BounceResult results;
  try
  {
    results = bounceTest("ch", "plane", params);
    success = true;
  }
  catch (SiconosException e)
  {
    std::cout << "SiconosException: " << e.report() << std::endl;
    CPPUNIT_ASSERT(1);
  }
}

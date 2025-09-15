#include "Contact2dTest.hpp"

#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"
#include "SiconosCollisionManager.hpp"
#include "SiconosBulletCollisionManager.hpp"
#include "RigidBody2dDS.hpp"
#include "SolverOptions.h"
#include "SiconosKernel.hpp"

#include <string>
#include <sys/time.h>

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
CPPUNIT_TEST_SUITE_REGISTRATION(Contact2dTest);

void Contact2dTest::setUp() {}
void Contact2dTest::tearDown() {}

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

  void dump() const
  {
    printf("  trace:              %s\n", trace?"on":"off");
    printf("  dynamic:            %s\n", trace?"on":"off");
    printf("  size:               %.3g\n", size);
    printf("  mass:               %.3g\n", mass);
    printf("  position:           %.3g\n", position);
    printf("  insideMargin:       %.3g\n", insideMargin);
    printf("  outsideMargin:      %.3g\n", outsideMargin);
    printf("  breakingThreshold:  %.3g\n", options.contactBreakingThreshold);
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

static
BounceResult bounceTest(std::string moving,
                        std::string ground,
                        const BounceParams &params)
{

  params.dump();
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
  SP::NonSmoothDynamicalSystem nsds(new NonSmoothDynamicalSystem(t0, T));

  SP::SiconosVector q0(new SiconosVector(3));
  SP::SiconosVector v0(new SiconosVector(3));
  q0->zero();
  v0->zero();
  (*q0)(1) = position_init;
  (*v0)(1) = velocity_init;

  SP::SiconosVector q1(new SiconosVector(3));
  SP::SiconosVector v1(new SiconosVector(3));
  q1->zero();
  (*q1)(1) = 1.0;
  v1->zero();

  /////////

  // Bodies

  printf("== Testing: %s falling on %s .. ",
         moving.c_str(), ground.c_str());
  if(params.trace) printf("==\n");
  fflush(stdout);

  // Set up a Siconos Mechanics environment:
  // A RigidBody2dDS with a contactor consisting of a single disk.

  SP::SiconosMatrix mass(new SimpleMatrix(3,3));
  mass->setValue(0,0,params.mass);
  mass->setValue(1,1,params.mass);
  mass->setValue(2,2,params.mass);
  SP::RigidBody2dDS body(new RigidBody2dDS(q0, v0, mass));
  SP::SiconosContactorSet contactors(new SiconosContactorSet());
  SP::SiconosDisk disk;
  if(moving=="disk")
  {
    disk.reset(new SiconosDisk(params.size));
    disk->setInsideMargin(params.insideMargin);
    disk->setOutsideMargin(params.outsideMargin);
    contactors->push_back(std::make_shared<SiconosContactor>(disk));
  }
  else if(moving=="box")
  {
    SP::SiconosBox box(new SiconosBox(params.size,params.size,params.size));
    box->setInsideMargin(params.insideMargin);
    box->setOutsideMargin(params.outsideMargin);
    contactors->push_back(std::make_shared<SiconosContactor>(box));
  }
  else if(moving=="ch2d")
  {
    float siz = params.size;
    SP::SiconosMatrix pts(new SimpleMatrix(4,2));
    (*pts)(0,0) = 0.0;
    (*pts)(0,1) = 0.0;
    (*pts)(1,1) = siz;
    (*pts)(1,1) = 0.0;
    (*pts)(2,0) = 0.0;
    (*pts)(2,1) = siz;
    SP::SiconosConvexHull2d ch2d(new SiconosConvexHull2d(pts));
    ch2d->setInsideMargin(params.insideMargin);
    ch2d->setOutsideMargin(params.outsideMargin);
    contactors->push_back(std::make_shared<SiconosContactor>(ch2d));
  }
  body->setContactors(contactors);

  // A contactor with no body (static contactor) consisting of a plane
  // positioned such that bouncing and resting position = 0.0
  SP::SiconosContactorSet static_contactors(std::make_shared<SiconosContactorSet>());

  if(ground=="disk")
  {
    SP::SiconosDisk floordisk(new SiconosDisk(params.size));
    floordisk->setInsideMargin(params.insideMargin);
    floordisk->setOutsideMargin(params.outsideMargin);
    SP::SiconosVector pos(new SiconosVector(7));
    pos->zero();
    (*pos)(1) =  -2.0;//-params.size/2;

    (*pos)(3) = 1.0; //unit quaternion
    static_contactors->push_back(std::make_shared<SiconosContactor>(floordisk, pos));
  }
  else if(ground=="box")
  {
    SP::SiconosBox2d floorbox(new SiconosBox2d(100,1.));
    floorbox->setInsideMargin(params.insideMargin);
    floorbox->setOutsideMargin(params.outsideMargin);
    SP::SiconosVector pos(new SiconosVector(7));
    (*pos)(1) = -2.0;//-params.size/2;
    (*pos)(3) = 1.0;
    static_contactors->push_back(std::make_shared<SiconosContactor>(floorbox, pos));
  }
  else if(ground=="sphere")
  {
    SP::SiconosSphere floorsphere(new SiconosSphere(1.0));
    floorsphere->setInsideMargin(params.insideMargin);
    floorsphere->setOutsideMargin(params.outsideMargin);
    SP::SiconosVector pos(new SiconosVector(7));
    (*pos)(1) = -1.0 -params.size/2;
    (*pos)(3) = 1.0;
    static_contactors->push_back(std::make_shared<SiconosContactor>(floorsphere, pos));
  }

  /////////

  // -- Set external forces (weight) --
  SP::SiconosVector FExt;
  FExt.reset(new SiconosVector(3));
  FExt->zero();
  FExt->setValue(1, - g * params.mass);
  body->setFExtPtr(FExt);

  // -- Add the dynamical systems into the non smooth dynamical system
  nsds->insertDynamicalSystem(body);

  // -- Time discretisation --
  SP::TimeDiscretisation timedisc(new TimeDiscretisation(t0, h));

  // -- OneStepNsProblem --
  SP::FrictionContact osnspb(new FrictionContact(2));

  // -- Some configuration
  osnspb->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 1000;
  osnspb->numericsSolverOptions()->dparam[SICONOS_DPARAM_TOL] = 1e-5;

  osnspb->setMaxSize(16384);
  osnspb->setMStorageType(NM_SPARSE_BLOCK);
  osnspb->setNumericsVerboseMode(0);
  osnspb->setKeepLambdaAndYState(true);

  /////////

  // --- Simulation initialization ---

  // -- MoreauJeanOSI Time Stepping
  SP::TimeStepping simulation(new TimeStepping(nsds, timedisc));

  simulation->insertIntegrator(osi);
  simulation->insertNonSmoothProblem(osnspb);


  // Object to manage the Bullet implementation of collisionMan
  SP::SiconosBulletCollisionManager collisionMan(
    new SiconosBulletCollisionManager(params.options));

  simulation->insertInteractionManager(collisionMan);

  // Add static shapes (centered at zero by default)
  collisionMan->addStaticBody(static_contactors);

  // Add a non-smooth law
  SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0.8, 0., 0.0, 2));
  collisionMan->insertNonSmoothLaw(nslaw, 0, 0);

  ///////

  int k=0;
  double std=0, final_pos=0;
  double last_pos=position_init, last_vel=0;
  while(simulation->hasNextEvent() and k <= 1000)
  {
    // // Update a property at step 500
    // if (params.dynamic && k==500 && moving=="disk") {
    //   disk->setRadius(0.3);
    // }

    // Update integrator and solve constraints
    simulation->computeOneStep();
    // osnspb->setNumericsVerboseMode(true);
    // osnspb->setNumericsVerboseLevel(1);
    // if (nsds->getNumberOfInteractions()>0)
    //   osnspb->display();
    double vel = (*body->velocity())(1);
    double pos = (*body->q())(1);

    if(params.trace && (k+1) < steps)
    {
      printf("k, %i, pos, %f, %f, %f\n", k,  pos, last_pos - pos, vel);
    }

    // Peaks (velocity crosses zero towards positive)
    if(vel <= 0 && last_vel > 0)
    {
      if(bounce_counter < n_desired_bounces)
      {
        actual_bounce_times[bounce_counter] = k*h;
        actual_bounce_ratios[bounce_counter] = pos / position_init;
        bounce_counter ++;
      }
    }

    // Interaction statistics
    if(collisionMan->statistics().new_interactions_created > 0 && first_contact)
    {
      first_contact = false;
      displacement_on_first_contact = last_pos - pos;
    }

    int interactions = collisionMan->statistics().new_interactions_created
                       + collisionMan->statistics().existing_interactions_processed;

    local_new_interaction_count += collisionMan->statistics().new_interactions_created;

    if(interactions > max_simultaneous_contacts)
      max_simultaneous_contacts = interactions;

    avg_simultaneous_contacts += interactions;

    // Standard deviation (cheating by not calculating mean!)
    if(k==(steps-100))
      final_pos = pos;
    if(k>(steps-100))
    {
      std += (pos-final_pos)*(pos-final_pos);
    }

    // State
    last_pos = pos;
    last_vel = vel;

    // Advance simulation
    simulation->nextStep();
    k++;
  }

  if(!params.trace)
    printf("(done, iterations=%d)\n", k - 1);

  BounceResult r;
  r.bounce_error_sum = 0.0;
  r.n_bounce_error = 6;
  for(int i=0; i < r.n_bounce_error; i++)
  {
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


void Contact2dTest::t1()
{
  try
  {
    printf("\n==== t1\n");

    BounceParams params;
    params.trace = true;
    params.dynamic = false;
    params.size = 1.0;
    params.mass = 1.0;
    params.position = 1.0;
    params.timestep = 0.005;
    params.insideMargin = 0.0;
    params.outsideMargin = 0.0;
    params.options.dimension = SICONOS_BULLET_2D;

    BounceResult r = bounceTest("disk", "disk", params);

    fprintf(stderr, "\nSize: %g\n", params.size);
    fprintf(stderr, "Final position: %g  (std=%g)\n\n",
            r.final_position, r.final_position_std);
  }
  catch(...)
  {
    siconos::exception::process();
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}
void Contact2dTest::t2()
{
  try
  {
    printf("\n==== t2\n");

    BounceParams params;
    params.trace = true;
    params.dynamic = false;
    params.size = 1.0;
    params.mass = 1.0;
    params.position = 1.0;
    params.timestep = 0.005;
    params.insideMargin = 0.0;
    params.outsideMargin = 0.0;
    params.options.dimension = SICONOS_BULLET_2D;

    BounceResult r = bounceTest("disk", "box", params);

    fprintf(stderr, "\nSize: %g\n", params.size);
    fprintf(stderr, "Final position: %g  (std=%g)\n\n",
            r.final_position, r.final_position_std);
  }
  catch(...)
  {
    siconos::exception::process();
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}
void Contact2dTest::t3()
{
  try
  {
    printf("\n==== t3\n");

    BounceParams params;
    params.trace = true;
    params.dynamic = false;
    params.size = 1.0;
    params.mass = 1.0;
    params.position = 1.0;
    params.timestep = 0.005;
    params.insideMargin = 0.0;
    params.outsideMargin = 0.0;
    params.options.dimension = SICONOS_BULLET_2D;

    BounceResult r = bounceTest("ch2d", "box", params);

    fprintf(stderr, "\nSize: %g\n", params.size);
    fprintf(stderr, "Final position: %g  (std=%g)\n\n",
            r.final_position, r.final_position_std);
  }
  catch(...)
  {
    siconos::exception::process();
    CPPUNIT_ASSERT(0);
  }

  CPPUNIT_ASSERT(1);
}

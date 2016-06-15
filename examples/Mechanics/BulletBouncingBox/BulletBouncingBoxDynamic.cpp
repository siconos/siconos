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

/*!\file BulletBouncingBoxDynamic.cpp
  \brief C++ input file, a Bullet box bouncing on the ground

  A box bouncing on the ground with the use of Bullet collision
  detection.

  This is similar to BulletBouncingBox.cpp, but adds a second box
  part-way through the simulation using
  BulletSpaceFilter::addDynamicObject(), to demonstrate its use.
*/

#include <SiconosBodies.hpp>
#include <SiconosKernel.hpp>

#include <BulletSpaceFilter.hpp>
#include <BulletTimeStepping.hpp>
#include <BulletDS.hpp>
#include <BulletWeightedShape.hpp>

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"
#pragma clang diagnostic ignored "-Woverloaded-virtual"
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
#endif
#include <BulletCollision/CollisionShapes/btConvexHullShape.h>
#include <BulletCollision/CollisionDispatch/btCollisionObject.h>
#include <BulletCollision/CollisionShapes/btBoxShape.h>
#include <BulletCollision/CollisionDispatch/btCollisionWorld.h>
#include <BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h>
#if defined(__clang__)
#pragma clang diagnostic pop
#elif !(__INTEL_COMPILER || __APPLE__ )
#pragma GCC diagnostic pop
#endif

SP::BulletDS makeBox(float g, float pos, std::vector<SP::BulletWeightedShape> shapes)
{
  // note: no rebound with a simple Bullet Box, why ?
  // the distance after the broadphase contact detection is negative
  // and then stay negative.
  // This is ok if we build one with btConveHullShape
  std11::shared_ptr<btConvexHullShape> box(new btConvexHullShape());
  box->addPoint(btVector3(-1.0, 1.0, -1.0));
  box->addPoint(btVector3(-1.0, -1.0, -1.0));
  box->addPoint(btVector3(-1.0, -1.0, 1.0));
  box->addPoint(btVector3(-1.0, 1.0, 1.0));
  box->addPoint(btVector3(1.0, 1.0, 1.0));
  box->addPoint(btVector3(1.0, 1.0, -1.0));
  box->addPoint(btVector3(1.0, -1.0, -1.0));
  box->addPoint(btVector3(1.0, -1.0, 1.0));

  SP::BulletWeightedShape box1(new BulletWeightedShape(box, 1.0));
  shapes.push_back(box1);

  SP::SiconosVector q0(new SiconosVector(7));
  SP::SiconosVector v0(new SiconosVector(6));
  v0->zero();
  q0->zero();

  (*q0)(2) = pos;
  (*q0)(3) = 1.0;
  (*v0)(2) = 0.0;

  // -- The dynamical system --
  // -- the default contactor is the shape given in the constructor
  // -- the contactor id is 0
  SP::BulletDS ds(new BulletDS(box1, q0, v0));

  // -- Set external forces (weight) --
  SP::SiconosVector FExt;
  FExt.reset(new SiconosVector(3)); //
  FExt->zero();
  FExt->setValue(2, - g * box1->mass());
  ds->setFExtPtr(FExt);

  return ds;
};

int main()
{

  // User-defined main parameters
  double t0 = 0;                   // initial computation time
  double T = 20.0;                 // end of computation time
  double h = 0.005;                // time step
  double position_init = 10.0;     // initial position
  double velocity_init = 0.0;      // initial velocity

  double g = 9.81;
  double theta = 0.5;              // theta for MoreauJeanOSI integrator

  // -----------------------------------------
  // --- Dynamical systems && interactions ---
  // -----------------------------------------

  try
  {

    // ------------
    // --- Init ---
    // ------------

    std::cout << "====> Model loading ..." << std::endl << std::endl;


    // -- OneStepIntegrators --
    SP::OneStepIntegrator osi;
    osi.reset(new MoreauJeanOSI(theta));

    // -- Model --
    SP::Model model(new Model(t0, T));

    std::vector<SP::BulletWeightedShape> shapes;

    // -- Create the first dynamical system here, but don't insert it
    //    until after the simulation starts.
    SP::BulletDS body( makeBox(g, position_init, shapes) );

    // -- Add the dynamical system in the non smooth dynamical system
    model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);

    SP::btCollisionObject ground(new btCollisionObject());
    ground->setCollisionFlags(btCollisionObject::CF_STATIC_OBJECT);
    SP::btCollisionShape groundShape(new btBoxShape(btVector3(30, 30, .5)));
    btMatrix3x3 basis;
    basis.setIdentity();
    ground->getWorldTransform().setBasis(basis);
    ground->setCollisionShape(&*groundShape);
    ground->getWorldTransform().getOrigin().setZ(-.50);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    SP::TimeDiscretisation timedisc(new TimeDiscretisation(t0, h));

    // -- OneStepNsProblem --
    SP::FrictionContact osnspb(new FrictionContact(3));

    // -- Some configuration

    osnspb->numericsSolverOptions()->iparam[0] = 1000; // Max number of
    // iterations
    osnspb->numericsSolverOptions()->dparam[0] = 1e-5; // Tolerance


    osnspb->setMaxSize(16384);                        // max number of
    // interactions

    osnspb->setMStorageType(1);                      // Sparse storage

    osnspb->setNumericsVerboseMode(0);               // 0 silent, 1
    // verbose

    osnspb->setKeepLambdaAndYState(true);            // inject
    // previous
    // solution

    // --- Simulation initialization ---

    std::cout << "====> Simulation initialisation ..." << std::endl << std::endl;

    int N = ceil((T - t0) / h); // Number of time steps


    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(0.8, 0., 0.0, 3));

    // -- The space filter performs broadphase collision detection
    SP::BulletSpaceFilter space_filter(new BulletSpaceFilter(model));

    // -- insert a non smooth law for contactors id 0
    space_filter->insert(nslaw, 0, 0);

    // -- add multipoint iterations, this is needed to gather at least
    // -- 3 contact points and avoid objects penetration, see Bullet
    // -- documentation
    space_filter->collisionConfiguration()->setConvexConvexMultipointIterations();
    space_filter->collisionConfiguration()->setPlaneConvexMultipointIterations();

    // -- The ground is a static object
    // -- we give it a group contactor id : 0
    space_filter->addStaticObject(ground, 0);

    // -- MoreauJeanOSI Time Stepping with Bullet Dynamical Systems
    SP::BulletTimeStepping simulation(new BulletTimeStepping(timedisc));

    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);

    model->initialize(simulation);

    std::cout << "====> End of initialisation ..." << std::endl << std::endl;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 4;
    SimpleMatrix dataPlot(N + 1, outputSize);
    dataPlot.zero();

    SP::SiconosVector q = body->q();
    SP::SiconosVector v = body->velocity();

    dataPlot(0, 0) = model->t0();
    dataPlot(0, 1) = (*q)(2);
    dataPlot(0, 2) = (*v)(2);

    // --- Time loop ---

    std::cout << "====> Start computation ... " << std::endl << std::endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (simulation->hasNextEvent())
    {
	  // --- Add a dynamic object at step 100 of the simulation ---
      if (k==100)
      {
	    space_filter->addDynamicObject( makeBox(g, 3, shapes), simulation );
      }

      space_filter->buildInteractions(model->currentTime());

      simulation->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  simulation->nextTime();
      dataPlot(k, 1) = (*q)(2);
      dataPlot(k, 2) = (*v)(2);

      // If broadphase collision detection shows some contacts then we may
      // display contact forces.
      if (space_filter->collisionWorld()->getDispatcher()->getNumManifolds() > 0)
      {

        // we *must* have an indexSet0, filled by Bullet broadphase
        // collision detection and an indexSet1, filled by
        // TimeStepping::updateIndexSet with the help of Bullet
        // getDistance() function
        if (model->nonSmoothDynamicalSystem()->topology()->numberOfIndexSet() == 2)
        {
          SP::InteractionsGraph index1 = simulation->indexSet(1);

          // This is the narrow phase contact detection : if
          // TimeStepping::updateIndexSet has filled indexSet1 then we
          // have some contact forces to display
          if (index1->size() > 0)
          {

            // Four contact points for a cube with a side facing the
            // ground. Note : changing Bullet margin for collision
            // detection may lead this assertion to be false.
            if (index1->size() == 4)
            {
              InteractionsGraph::VIterator iur = index1->begin();

              // different version of bullet may not gives the same
              // contact points! So we only keep the summation.
              dataPlot(k, 3) =
                index1->bundle(*iur)-> lambda(1)->norm2() +
                index1->bundle(*++iur)->lambda(1)->norm2() +
                index1->bundle(*++iur)->lambda(1)->norm2() +
                index1->bundle(*++iur)->lambda(1)->norm2();
            }
          }
        }
      }

      simulation->nextStep();
      ++show_progress;
      k++;
    }


    std::cout << std::endl << "End of computation - Number of iterations done: " << k - 1 << std::endl;
    std::cout << "Computation Time " << time.elapsed()  << std::endl;

    // --- Output files ---
    std::cout << "====> Output file writing ..." << std::endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result_dynamic.dat", "ascii", dataPlot, "noDim");

    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result_dynamic.ref", "ascii", dataPlotRef);

    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file : "
                << (dataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }



  }

  catch (SiconosException e)
  {
    std::cout << e.report() << std::endl;
    exit(1);
  }
  catch (...)
  {
    std::cout << "Exception caught in BulletBouncingBox" << std::endl;
    exit(1);
  }

  return 0;
}

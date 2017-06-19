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
  part-way through the simulation to demonstrate dynamic changes to
  the graph.
*/

#include <SiconosBodies.hpp>
#include <SiconosKernel.hpp>

#include <SiconosBulletCollisionManager.hpp>
#include <BodyDS.hpp>

SP::BodyDS makeBox(float g, float pos, float vel)
{
  // -- Shape: cube with all dimensions=1.0
  SP::SiconosBox box1(std11::make_shared<SiconosBox>(1.0, 1.0, 1.0));

  // -- Initial position and velocity
  SP::SiconosVector q0(std11::make_shared<SiconosVector>(7));
  SP::SiconosVector v0(std11::make_shared<SiconosVector>(6));
  v0->zero();
  q0->zero();

  (*q0)(2) = pos;
  (*q0)(3) = 1.0;
  (*v0)(2) = vel;

  // -- The dynamical system --
  SP::BodyDS body(std11::make_shared<BodyDS>(q0, v0, 1.0));

  // -- add the box to the body's set of contactactors
  // -- by default, the contactor id is 0 with no position offset,
  //    see SiconosContactor.hpp for how to change these.
  body->contactors()->push_back(std11::make_shared<SiconosContactor>(box1));

  // -- Set external forces (weight) --
  SP::SiconosVector FExt(std11::make_shared<SiconosVector>(3));
  FExt->zero();
  FExt->setValue(2, - g * body->scalarMass());
  body->setFExtPtr(FExt);

  return body;
}

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

    // -- Moving object --
    SP::BodyDS body(makeBox(g, position_init, velocity_init));

    // -- Add the dynamical system in the non smooth dynamical system
    model->nonSmoothDynamicalSystem()->insertDynamicalSystem(body);

    SP::SiconosPlane ground(std11::make_shared<SiconosPlane>());

    // -- Create a Z-offset of -0.5 for the ground so that contact is
    //    at zero.
    SP::SiconosVector groundOffset(std11::make_shared<SiconosVector>(7));
    (*groundOffset)(2) = -0.5;  // translation 0,0,-0.5
    (*groundOffset)(3) = 1;     // orientation 1,0,0,0

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

    // some options for the Bullet collision manager:
    // -- defaults are okay, see SiconosBulletCollisionManager.hpp
    // -- in particular we want to leave multipoint iterations enabled
    //    to allow Bullet to collect more points for plane-plane
    //    collisions.
    SiconosBulletOptions options;

    // -- The collision manager performs broadphase collision
    //    detection, we use the Bullet implementation here.
    SP::SiconosBulletCollisionManager collision_manager(
      std11::make_shared<SiconosBulletCollisionManager>(options));

    // -- insert a non smooth law for contactors id 0
    collision_manager->insertNonSmoothLaw(nslaw, 0, 0);

    // -- The ground is a static object.  The collision manager
    //    maintains a list of contact sets for static objects, so we
    //    add one.
    // -- We give it a group contactor id : 0
    // -- We apply the groundOffset to the SiconosContactorSet, but
    //    equivalently it could be applied to the SiconosContactor
    //    inside the set, this is a design choice allowing for re-use
    //    for more complex compound contactor sets.
    SP::SiconosContactorSet staticCtrSet(std11::make_shared<SiconosContactorSet>());
    staticCtrSet->push_back(std11::make_shared<SiconosContactor>(ground));
    collision_manager->insertStaticContactorSet(staticCtrSet, groundOffset);

    // -- MoreauJeanOSI Time Stepping with Bullet collision manager as
    // -- the interaction manager.
    SP::TimeStepping simulation(new TimeStepping(timedisc));
    simulation->insertInteractionManager(collision_manager);

    simulation->insertIntegrator(osi);
    simulation->insertNonSmoothProblem(osnspb);
    model->setSimulation(simulation);

    model->initialize();

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
        SP::BodyDS ds(makeBox(g, 3.0, 0));
        simulation->nonSmoothDynamicalSystem()->insertDynamicalSystem(ds);
        simulation->prepareIntegratorForDS(osi, ds, model, simulation->nextTime());
      }

      collision_manager->resetStatistics();
      simulation->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  simulation->nextTime();
      dataPlot(k, 1) = (*q)(2);
      dataPlot(k, 2) = (*v)(2);

      // If broadphase collision detection shows some contacts then we may
      // display contact forces.
      if ((collision_manager->statistics().new_interactions_created
           + collision_manager->statistics().existing_interactions_processed) > 0)
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

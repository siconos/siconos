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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 3.0;
    double h = 0.001;               // time step
    double position_init = 0.1;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.1; // Ball radius
    double m1 = 1; // Ball mass
    double m2 = 2.0; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<  endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m1;
    (*Mass)(1, 1) = m1;
    (*Mass)(2, 2) = 2. / 5 * m1 * R * R;

    SP::SiconosMatrix Mass2(new SimpleMatrix(nDof, nDof));
    (*Mass2)(0, 0) = m2;
    (*Mass2)(1, 1) = m2;
    (*Mass2)(2, 2) = 2. / 5 * m2 * R * R;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;
    
    SP::SiconosVector q0_2(new SiconosVector(nDof));
    SP::SiconosVector v0_2(new SiconosVector(nDof));
    (*q0_2)(0) = position_init  + 2*R + 0.001;
    (*v0_2)(0) = velocity_init;

    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball1(new LagrangianLinearTIDS(q0, v0, Mass));
    SP::LagrangianLinearTIDS ball2(new LagrangianLinearTIDS(q0_2, v0_2, Mass2));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m1 * g;
    ball1->setFExtPtr(weight);

    SP::SiconosVector weight2(new SiconosVector(nDof));
    (*weight2)(0) = -m2 * g;
    ball2->setFExtPtr(weight2);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.95;

    // Interaction ball-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(nDof, nDof));
    (*H)(0, 0) = 1.0;
    (*H)(1, 1) = 1.0;
    (*H)(2, 2) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(e, e, 0.6, 3));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));
    
    // Interaction ball-ball
    //
    SP::SimpleMatrix H_bb(new SimpleMatrix(nDof, 2*nDof));
    (*H_bb)(0, 0) = -1.0;
    (*H_bb)(1, 1) = -1.0;
    (*H_bb)(2, 2) = -1.0;
    (*H_bb)(0, 3) = 1.0;
    (*H_bb)(1, 4) = 1.0;
    (*H_bb)(2, 5) = 1.0;
    
    SP::SiconosVector b_bb(new SiconosVector(3));
    (*b_bb)(0) = -2 * R;
    (*b_bb)(1) = 0.0;
    (*b_bb)(2) = 0.0;

    SP::Relation relation_bb(new LagrangianLinearTIR(H_bb,b_bb));

    SP::Interaction inter_bb(new Interaction(nslaw, relation_bb));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem bouncingBall(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->insertDynamicalSystem(ball1);
    bouncingBall->insertDynamicalSystem(ball2);

    // link the interaction and the dynamical system
    bouncingBall->link(inter, ball1);
    bouncingBall->link(inter_bb, ball1, ball2);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanGOSI OSI(new MoreauJeanGOSI(theta));


    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new GlobalFrictionContact(3,SICONOS_GLOBAL_FRICTION_3D_NSGS_WR));
    //SP::OneStepNSProblem osnspb(new GlobalFrictionContact(3,SICONOS_GLOBAL_FRICTION_3D_NSN_AC));
    //SP::OneStepNSProblem osnspb(new GlobalFrictionContact(3));
    assert(osnspb->numericsSolverOptions());
    SolverOptions * options = osnspb->numericsSolverOptions().get();
    //solver_options_print(options);
    options->internalSolvers->dparam[0] = 1e-10;
    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(bouncingBall, t, OSI, osnspb));
 
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q1 = ball1->q();
    SP::SiconosVector v1 = ball1->velocity();
    SP::SiconosVector p1 = ball1->p(1);
    SP::SiconosVector lambda = inter->lambda(1);
    SP::SiconosVector q2 = ball2->q();
    SP::SiconosVector v2 = ball2->velocity();
    SP::SiconosVector p2 = ball2->p(1);
    // SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q1)(0);
    dataPlot(0, 2) = (*v1)(0);
    dataPlot(0, 3) = (*p1)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*q2)(0);
    dataPlot(0, 6) = (*v2)(0);
    dataPlot(0, 7) = (*p2)(0);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (s->hasNextEvent())
    {
      osnspb->setNumericsVerboseMode(0);

      //std::cout << "############################  time step = " << k <<std::endl;
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*v1)(0);
      dataPlot(k, 3) = (*p1)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q2)(0);
      dataPlot(k, 6) = (*v2)(0);
      dataPlot(k, 7) = (*p2)(0);
      //osnspb->display();
      s->nextStep();
      ++show_progress;
      k++;
    }
    cout  << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

    double error=0.0, eps=1e-12;
    if ((error=ioMatrix::compareRefFile(dataPlot,
                                        "TwoBouncingBallTS-MoreauJeanGOSI.ref",
                                        eps)) >= 0.0
        && error > eps)
      return 1;

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingBallTS.cpp" << endl;
  }



}

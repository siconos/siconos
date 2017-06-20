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
    double T = 10;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 2. / 5 * m * R * R;

    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    (*q0)(0) = position_init;
    (*v0)(0) = velocity_init;

    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));

    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
    (*H)(0, 0) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::Model bouncingBall(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->insertDynamicalSystem(ball);

    // link the interaction and the dynamical system
    bouncingBall->nonSmoothDynamicalSystem()->link(inter, ball);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanCombinedProjectionOSI OSI(new MoreauJeanCombinedProjectionOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());
    SP::OneStepNSProblem osnspb_pos(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeSteppingCombinedProjection s(new TimeSteppingCombinedProjection(t, OSI, osnspb, osnspb_pos));
    s->setProjectionMaxIteration(4);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);
    bouncingBall->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize();
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 7;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();

    SP::SiconosVector p1 = ball->p(1);
    SP::SiconosVector lambda1 = inter->lambda(1);
    SP::SiconosVector lambda0 = inter->lambda(0);
    SP::SiconosVector p0 = ball->p(0);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p1)(0);
    dataPlot(0, 4) = (*lambda1)(0);
    dataPlot(0, 5) = (*p0)(0);
    dataPlot(0, 6) = (*lambda0)(0);

    double maxviolation = std::max(0.0,  -(*y0)(0));
    //dataPlot(0, 6) = (*lambda2)(0);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    //while (s->hasNextEvent() && k < 95)
    while (s->hasNextEvent())
    {

      // for( int toto=0; toto<3;toto++)
      //  std::cout <<"============> Step Number : " << k << "======================"  <<std::endl;



      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p1)(0);
      dataPlot(k, 4) = (*lambda1)(0);
      dataPlot(k, 5) = (*p0)(0);
      dataPlot(k, 6) = (*lambda0)(0);

      maxviolation = s->maxViolationUnilateral();

      if (maxviolation >= 1e-08)
      {
        std::cout << std::endl << "maxviolation = " << maxviolation << std::endl;
        std::cout << "(*lambda1)(0) "  << (*lambda1)(0) << std::endl;
        std::cout << "(*lambda0)(0) "  << (*lambda0)(0) << std::endl;
        //RuntimeException::selfThrow("maxviolation > 0   ");
      }

      //dataPlot(k, 6) = (*lambda2)(0);
      s->nextStep();
      //osnspb->display();
      //osnspb_pos->display();
      //std::cout <<" (*lambda1)(0) "  <<(*lambda1)(0)<< std::endl;
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;
    cout << "Max violation " << maxviolation << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    cout << "====> Comparison with a reference file ..." << endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result-CombinedProjection.ref", "ascii", dataPlotRef);
    std::cout << "Error w.r.t reference file = " << (dataPlot - dataPlotRef).normInf() << std::endl;

    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
      std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
      (dataPlot - dataPlotRef).display();
      return 1;
    }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingBallTS.cpp" << endl;
  }



}

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

/*!\file
  C++ input file, D1MinusLinearOSI-Time-Stepping version
  T. Schindler, V. Acary

  A Ball bouncing on the ground.
  Direct description of the model
  Simulation with a D1MinusLinearOSI-Time-Stepping scheme.
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
    double T = 10.;                  // final computation time
    double h = 5e-4;                 // time step
    double hplot = 0.005;            // plot step size (larger than time step)
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double R = 0.1;                  // Ball radius
    double m = 1;                    // Ball mass
    double g = 9.81;                 // Gravity

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
    SP::D1MinusLinearOSI OSI(new D1MinusLinearOSI());

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) One step non smooth problem
    SP::OneStepNSProblem impact(new LCP()); // impulse right limit right side
    SP::OneStepNSProblem force(new LCP()); // contact force right limit left side

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeSteppingD1Minus s(new TimeSteppingD1Minus(t, 2));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, SICONOS_OSNSP_TS_VELOCITY + 1);
    bouncingBall->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize();
    int N = ceil((T - t0) / h); // Number of time steps
    int Nplot = (int)((T - t0) / hplot); // Number of plot steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 7;
    SimpleMatrix dataPlot(N + 10, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector lambda = inter->lambda(1);
    SP::SiconosVector p2 = ball->p(2);
    SP::SiconosVector lambda2 = inter->lambda(2);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = (*p2)(0);
    dataPlot(0, 6) = (*lambda2)(0);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (s->hasNextEvent())
    {
      s->advanceToEvent();
      // ball->display();
      // --- Get values to be plotted ---
      //  if (fmod(s->nextTime(), hplot) < h)
      {
 
        //std::cout << "k=" << k <<std::endl;
        dataPlot(k, 0) =  s->nextTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*v)(0);
        dataPlot(k, 3) = (*p)(0);
        dataPlot(k, 4) = (*lambda)(0);  
        dataPlot(k, 5) = (*p2)(0);
        dataPlot(k, 6) = (*lambda2)(0);
        k++;
      }

      s->processEvents();
      ++show_progress;
    }

    cout << endl << "End of computation - Number of iterations done: " << N - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result_tdg.dat", "ascii", dataPlot, "noDim");
    cout << "====> Comparison with a reference file ..." << endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("result_tdg.ref", "ascii", dataPlotRef);
    double error = (dataPlot - dataPlotRef).normInf()/ dataPlotRef.normInf();
    std::cout << "Error = "<< error << std::endl;
    if (error > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
      std::cout <<  "error  = " << (dataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingBallD1MinusLinearOSI.cpp" << endl;
  }
}

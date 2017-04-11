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
  C++ input file, MoreauJeanOSI-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a MoreauJeanOSI-Time-Stepping scheme

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody systems with
  translational clearance joints based on the non-smooth dynamics approach
  */

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3; // degrees of freedom for robot arm
    double t0 = 0;         // initial computation time
    double T = 0.2;       // final computation time
    double h = 1e-4;       // time step : do not decrease, because of strong penetrations

    // geometrical characteristics
    double l1 = 0.1530;
    double l2 = 0.3060;
    double a = 0.05;
    double b = 0.025;
    double c = 0.001;

    // contact parameters
    double eN1 = 0.4;
    double eN2 = 0.4;
    double eN3 = 0.4;
    double eN4 = 0.4;
    double eT1 = 0.;
    double eT2 = 0.;
    double eT3 = 0.;
    double eT4 = 0.;
    double mu1 = 0.01;
    double mu2 = 0.01;
    double mu3 = 0.01;
    double mu4 = 0.01;

    // initial conditions
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();
    (*v0)(0) = 150.;
    (*v0)(1) = -75.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    SP::LagrangianDS slider(new LagrangianDS(q0, v0, "SliderCrankPlugin:mass"));
    slider->setComputeFGyrFunction("SliderCrankPlugin", "FGyr");
    slider->setComputeJacobianFGyrqFunction("SliderCrankPlugin", "jacobianFGyrq");
    slider->setComputeJacobianFGyrqDotFunction("SliderCrankPlugin", "jacobianFGyrqDot");
    slider->setComputeFIntFunction("SliderCrankPlugin", "FInt");
    slider->setComputeJacobianFIntqFunction("SliderCrankPlugin", "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("SliderCrankPlugin", "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eN1, eT1, mu1, 2));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1", "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(eN2, eT2, mu2, 2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2", "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactFrictionNSL(eN3, eT3, mu3, 2));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3", "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction(nslaw3, relation3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactFrictionNSL(eN4, eT4, mu4, 2));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4", "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction(nslaw4, relation4));

    // -------------
    // --- Model ---
    // -------------
    SP::Model sliderWithClearance(new Model(t0, T));
    sliderWithClearance->nonSmoothDynamicalSystem()->insertDynamicalSystem(slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter1, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter2, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter3, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter4, slider);

    // ----------------
    // --- Simulation ---
    // ----------------
    SP::MoreauJeanCombinedProjectionOSI OSI(new MoreauJeanCombinedProjectionOSI(0.5));
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::OneStepNSProblem impact(new FrictionContact(2, SICONOS_FRICTION_2D_ENUM));
    impact->numericsSolverOptions()->dparam[0] = 1e-08;
    impact->numericsSolverOptions()->iparam[0] = 100;
    impact->numericsSolverOptions()->iparam[2] = 1; // random
    SP::OneStepNSProblem position(new MLCPProjectOnConstraints());

    SP::TimeSteppingCombinedProjection s(new TimeSteppingCombinedProjection(t, OSI, impact, position, 2));
    s->setProjectionMaxIteration(500);
    s->setConstraintTolUnilateral(1e-10);
    s->setConstraintTol(1e-10);
    sliderWithClearance->setSimulation(s);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;
    sliderWithClearance->initialize();
    int N = ceil((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 35;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = slider->q();
    SP::SiconosVector v = slider->velocity();

    dataPlot(0, 0) = sliderWithClearance->t0();
    dataPlot(0, 1) = (*q)(0) / (2.*M_PI); // crank revolution
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);
    dataPlot(0, 7) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 1 (normalized)
    dataPlot(0, 8) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 2 (normalized)
    dataPlot(0, 9) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (-c); // y corner 3 (normalized)
    dataPlot(0, 10) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (-c); // y corner 4 (normalized)
    dataPlot(0, 11) = (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1; // x slider (normalized)
    dataPlot(0, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c; // y slider (normalized
    dataPlot(0, 13) = (*inter1->y(0))(0) ; // g1
    dataPlot(0, 14) = (*inter2->y(0))(0) ; // g2
    dataPlot(0, 15) = (*inter3->y(0))(0) ; // g3
    dataPlot(0, 16) = (*inter4->y(0))(0) ; // g4
    dataPlot(0, 17) = (*inter1->y(1))(0) ; // dot g1
    dataPlot(0, 18) = (*inter2->y(1))(0) ; // dot g2
    dataPlot(0, 19) = (*inter3->y(1))(0) ; // dot g3
    dataPlot(0, 20) = (*inter4->y(1))(0) ; // dot g4
    dataPlot(0, 21) = (*inter1->lambda(1))(0) ; // lambda1
    dataPlot(0, 22) = (*inter2->lambda(1))(0) ; // lambda2
    dataPlot(0, 23) = (*inter3->lambda(1))(0) ; // lambda3
    dataPlot(0, 24) = (*inter4->lambda(1))(0) ; // lambda4
    dataPlot(0, 25) = (*inter1->lambda(0))(0) ; // lambda1 projection
    dataPlot(0, 26) = (*inter2->lambda(0))(0) ; // lambda2
    dataPlot(0, 27) = (*inter3->lambda(0))(0) ; // lambda3
    dataPlot(0, 28) = (*inter4->lambda(0))(0) ; // lambda4
    dataPlot(0, 29) = 1;
    dataPlot(0, 30) = 0;
    dataPlot(0, 31) = 0;
    dataPlot(0, 32) = 0;
    dataPlot(0, 33) = 0;
    dataPlot(0, 34) = 0;
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

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) / (2.*M_PI); // crank revolution
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 1 (normalized)
      dataPlot(k, 8) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 2 (normalized)
      dataPlot(k, 9) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (-c); // y corner 3 (normalized)
      dataPlot(k, 10) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (-c); // y corner 4 (normalized)
      dataPlot(k, 11) = (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1; // x slider (normalized)
      dataPlot(k, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c; // y slider (normalized)
      dataPlot(k, 13) = (*inter1->y(0))(0) ; // g1
      dataPlot(k, 14) = (*inter2->y(0))(0) ; // g2
      dataPlot(k, 15) = (*inter3->y(0))(0) ; // g3
      dataPlot(k, 16) = (*inter4->y(0))(0) ; // g4
      dataPlot(k, 17) = (*inter1->y(1))(0) ; // dot g1
      dataPlot(k, 18) = (*inter2->y(1))(0) ; // dot g2
      dataPlot(k, 19) = (*inter3->y(1))(0) ; // dot g3
      dataPlot(k, 20) = (*inter4->y(1))(0) ; // dot g4
      dataPlot(k, 21) = (*inter1->lambda(1))(0) ; // lambda1
      dataPlot(k, 22) = (*inter2->lambda(1))(0) ; // lambda2
      dataPlot(k, 23) = (*inter3->lambda(1))(0) ; // lambda3
      dataPlot(k, 24) = (*inter4->lambda(1))(0) ; // lambda4
      dataPlot(k, 25) = (*inter1->lambda(0))(0) ; // lambda1 projection
      dataPlot(k, 26) = (*inter2->lambda(0))(0) ; // lambda2
      dataPlot(k, 27) = (*inter3->lambda(0))(0) ; // lambda3
      dataPlot(k, 28) = (*inter4->lambda(0))(0) ; // lambda4
      dataPlot(k, 29) = s->getNewtonNbIterations();
      dataPlot(k, 30) = s->nbProjectionIteration();
      dataPlot(k, 31) = s->maxViolationUnilateral();
      dataPlot(k, 32) = s->nbIndexSetsIteration();
      dataPlot(k, 33) = s->cumulatedNewtonNbIterations();
      dataPlot(k, 34) = s->nbCumulatedProjectionIteration();
      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << endl << "Max violation unilateral = " << s->maxViolationUnilateral() << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    cout << "====> Comparison with a reference file ..." << endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("SliderCrankMoreauJeanCombinedProjectionOSI.ref", "ascii", dataPlotRef);
    double error = (dataPlot - dataPlotRef).normInf()/ dataPlotRef.normInf();
    std::cout << "Error = "<< error << std::endl;
    if (error > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
      std::cout <<  "Absolute error  = " << (dataPlot - dataPlotRef).normInf() << std::endl;
      SP::SiconosVector err(new SiconosVector(dataPlot.size(1)));
      (dataPlot - dataPlotRef).normInfByColumn(err);
      err->display();
      return 1;
    }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in SliderCrankD1MinusLinearOSI.cpp" << endl;
  }
}

/* Siconos-sample , Copyright INRIA 2005-2011.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */


/*!\file
  C++ input file, Moreau-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a Moreau-Time-Stepping scheme

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
    double T = 0.15;       // final computation time
    double h = 1e-5;       // time step : do not decrease, because of strong penetrations

    // geometrical characteristics
    double l1 = 0.1530;
    double l2 = 0.3060;
    double a = 0.05;
    double b = 0.025;
    double c = 0.001;

    // contact parameters
    double e1 = 0.4;
    double e2 = 0.4;
    double e3 = 0.4;
    double e4 = 0.4;
    //double mu1 = 0.01;
    //double mu2 = 0.01;
    //double mu3 = 0.01;
    //double mu4 = 0.01;

    // initial conditions
    SP::SimpleVector q0(new SimpleVector(nDof));
    SP::SimpleVector v0(new SimpleVector(nDof));
    q0->zero();
    v0->zero();
    (*v0)(0) = 150.;
    (*v0)(1) = -75.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    SP::LagrangianDS slider(new LagrangianDS(q0, v0, "SliderCrankPlugin:mass"));
    slider->setComputeNNLFunction("SliderCrankPlugin.so", "NNL");
    slider->setComputeJacobianNNLqFunction("SliderCrankPlugin.so", "jacobianNNLq");
    slider->setComputeJacobianNNLqDotFunction("SliderCrankPlugin.so", "jacobianNNLqDot");
    slider->setComputeFIntFunction("SliderCrankPlugin.so", "FInt");
    slider->setComputeJacobianFIntqFunction("SliderCrankPlugin.so", "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("SliderCrankPlugin.so", "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
    SP::NonSmoothLaw nslaw1(new NewtonImpactNSL(e1));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1", "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction(1, nslaw1, relation1, 1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(e2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2", "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction(1, nslaw2, relation2, 2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactNSL(e3));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3", "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction(1, nslaw3, relation3, 3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactNSL(e4));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4", "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction(1, nslaw4, relation4, 4));

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
    SP::Moreau OSI(new Moreau(slider, 0.5));
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem force(new LCP());

    SP::TimeStepping s(new TimeStepping(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;
    sliderWithClearance->initialize(s);
    int N = (int)((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 13;
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

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (s->nextTime() < T)
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

      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    dataPlot.resize(k, outputSize);
    io.write(dataPlot, "noDim");
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in SliderCrankD1MinusLinear.cpp" << endl;
  }
}

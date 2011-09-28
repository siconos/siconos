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
  C++ input file, D1MinusLinear-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a D1MinusLinear-Time-Stepping scheme

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
    double T = 1.0;        // final computation time
    double h = 0.0005;     // time step

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

    DynamicalSystemsSet allDS;
    SP::LagrangianDS slider(new LagrangianDS(q0, v0, "SliderCrankPlugin:mass"));
    slider->setComputeNNLFunction("SliderCrankPlugin.so", "NNL");
    slider->setComputeJacobianNNLqFunction("SliderCrankPlugin.so", "jacobianNNLq");
    slider->setComputeJacobianNNLqDotFunction("SliderCrankPlugin.so", "jacobianNNLqDot");
    slider->setComputeFIntFunction("SliderCrankPlugin.so", "FInt");
    slider->setComputeJacobianFIntqFunction("SliderCrankPlugin.so", "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("SliderCrankPlugin.so", "jacobianFIntqDot");
    allDS.insert(slider);

    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;

    // -- corner 1 --
    SP::NonSmoothLaw nslaw1(new NewtonImpactNSL(e1));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1", "SliderCrankPlugin:W1"));
    SP::Interaction inter1(new Interaction("corner1", allDS, 1, 1, nslaw1, relation1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(e2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2", "SliderCrankPlugin:W2"));
    SP::Interaction inter2(new Interaction("corner2", allDS, 2, 1, nslaw2, relation2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactNSL(e3));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3", "SliderCrankPlugin:W3"));
    SP::Interaction inter3(new Interaction("corner3", allDS, 3, 1, nslaw3, relation3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactNSL(e4));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4", "SliderCrankPlugin:W4"));
    SP::Interaction inter4(new Interaction("corner4", allDS, 4, 1, nslaw4, relation4));

    allInteractions.insert(inter1);
    allInteractions.insert(inter2);
    allInteractions.insert(inter3);
    allInteractions.insert(inter4);

    // -------------
    // --- Model ---
    // -------------
    SP::Model SliderWithClearance(new Model(t0, T, allDS, allInteractions));

    // ----------------
    // --- Simulation ---
    // ----------------
    SP::D1MinusLinear OSI(new D1MinusLinear(slider));
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem force(new LCP());

    SP::TimeSteppingD1Minus s(new TimeSteppingD1Minus(t, 2));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, SICONOS_OSNSP_TS_VELOCITY + 1);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;
    SliderWithClearance->initialize(s);
    int N = (int)((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 7;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = slider->q();
    SP::SiconosVector v = slider->velocity();

    dataPlot(0, 0) =  SliderWithClearance->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*q)(1);
    dataPlot(0, 3) = (*q)(2);
    dataPlot(0, 4) = (*v)(0);
    dataPlot(0, 5) = (*v)(1);
    dataPlot(0, 6) = (*v)(2);

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
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);

      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("SliderCrankResult.dat", "ascii");
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

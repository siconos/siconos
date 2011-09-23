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

/*!\file BouncingBallTS.cpp
  \brief \ref EMBouncingBall - C++ input file, Time-Stepping version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model without XML input.
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"
#include "KneeJointR.hpp"
#include "PrismaticJointR.hpp"
#include <math.h>
using namespace std;

int main(int argc, char* argv[])
{
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;
    unsigned int qDim = 7;
    unsigned int nDim = 6;
    double t0 = 0;                   // initial computation time
    double h = 0.001;                // time step
    double T = 10;
    double L1 = 1.0;
    double L2 = 2.0;
    double L3 = 1.0;
    double theta = 1.0;              // theta for Moreau integrator
    double g = 9.81; // Gravity
    double m = 1.;
    double wx = 0.0;
    double wz = 0.0;
    double wy = 0.0;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    DynamicalSystemsSet allDS1; // the list of DS

    // -- Initial positions and velocities --
    SP::SimpleVector q10(new SimpleVector(qDim));
    SP::SimpleVector v10(new SimpleVector(nDim));
    SP::SimpleMatrix I1(new SimpleMatrix(3, 3));
    v10->zero();
    I1->eye();
    q10->zero();
    q10->setValue(0, 1);
    q10->setValue(1, 1);
    q10->setValue(2, 1);

    double angle = M_PI / 5;
    SimpleVector V1(3);
    V1.zero();
    V1.setValue(2, 1);
    V1.setValue(0, 3);
    V1.setValue(1, 2);
    double Vnorm = V1.norm2();
    V1.setValue(1, V1.getValue(1) / Vnorm);
    V1.setValue(0, V1.getValue(0) / Vnorm);
    V1.setValue(2, V1.getValue(2) / Vnorm);
    q10->setValue(3, cos(angle));
    q10->setValue(4, V1.getValue(0)*sin(angle));
    q10->setValue(5, V1.getValue(1)*sin(angle));
    q10->setValue(6, V1.getValue(2)*sin(angle));


    // -- The dynamical system --
    SP::NewtonEulerDS beam1(new NewtonEulerDS(q10, v10, m, I1));
    allDS1.insert(beam1);
    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof));
    (*weight)(2) = -m * g;
    beam1->setFExtPtr(weight);





    // --------------------
    // --- Interactions ---
    // --------------------

    InteractionsSet allInteractions;


    // Interaction ball-floor
    //
    SP::SimpleMatrix H1(new SimpleMatrix(PrismaticJointR::_sNbEqualities, qDim));
    H1->zero();
    SP::NonSmoothLaw nslaw1(new EqualityConditionNSL(PrismaticJointR::_sNbEqualities));
    SP::SimpleVector axe1(new SimpleVector(3));
    axe1->zero();
    axe1->setValue(2, 1);

    SP::NewtonEulerR relation1(new PrismaticJointR(beam1, axe1));
    relation1->setJachq(H1);
    SP::Interaction inter1(new Interaction("axis-beam1", allDS1, 0, PrismaticJointR::_sNbEqualities, nslaw1, relation1));
    allInteractions.insert(inter1);
    // -------------
    // --- Model ---
    // -------------
    SP::Model bouncingBall(new Model(t0, T, allDS1, allInteractions));

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));
    //    s->setComputeResiduY(true);
    //  s->setUseRelativeConvergenceCriteron(false);

    // -- OneStepIntegrators --
    SP::Moreau OSI1(new Moreau(beam1, theta));
    s->insertIntegrator(OSI1);


    string solverName = "LinearSystem" ;
    // -- OneStepNsProblem --
    SP::OneStepNSProblem osnspb(new Equality());
    s->insertNonSmoothProblem(osnspb);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize(s);
    int N = 2000; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector q1 = beam1->q();
    std::cout << "computeh1\n";
    relation1->computeh(0.);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    int NewtonIt = 0;
    Index dimIndex(2);
    Index startIndex(4);
    int cmp = 0;
    for (cmp = 0; cmp < N; cmp++)
    {
      // solve ...
      s->newtonSolve(1e-4, 50);

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*q1)(1);
      dataPlot(k, 3) = (*q1)(2);
      dataPlot(k, 4) = (*q1)(3);
      dataPlot(k, 5) = (*q1)(4);
      dataPlot(k, 6) = (*q1)(5);
      dataPlot(k, 7) = (*q1)(6);

      s->nextStep();
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix ref("prismatic.ref", "ascii");
    ref.read(dataPlotRef);
    if ((dataPlot - dataPlotRef).normInf() > 1e-7)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
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

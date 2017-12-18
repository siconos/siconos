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

/*!\file PrismaticTest.cpp
  \brief \ref 

  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"
#include "KneeJointR.hpp"
#include "PrismaticJointR.hpp"
#include <boost/math/quaternion.hpp>
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
    double theta = 1.0;              // theta for MoreauJeanOSI integrator
    double g = 9.81; // Gravity
    double m = 1.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    // -- Initial positions and velocities --
    SP::SiconosVector q10(new SiconosVector(qDim));
    SP::SiconosVector v10(new SiconosVector(nDim));
    SP::SimpleMatrix I1(new SimpleMatrix(3, 3));
    v10->zero();
    I1->eye();
    q10->zero();
    q10->setValue(0, 1);
    q10->setValue(1, 1);
    q10->setValue(2, 1);

    double angle = M_PI / 5;
    SiconosVector V1(3);
    V1.zero();
    V1.setValue(0, 3);
    V1.setValue(1, 2);
    V1.setValue(2, 1);
    double Vnorm = V1.norm2();
    V1.setValue(0, V1.getValue(0) / Vnorm);
    V1.setValue(1, V1.getValue(1) / Vnorm);
    V1.setValue(2, V1.getValue(2) / Vnorm);
    q10->setValue(3, cos(angle));
    q10->setValue(4, V1.getValue(0)*sin(angle));
    q10->setValue(5, V1.getValue(1)*sin(angle));
    q10->setValue(6, V1.getValue(2)*sin(angle));

    // -- The dynamical system --
    SP::NewtonEulerDS beam1(new NewtonEulerDS(q10, v10, m, I1));
    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(2) = -m * g;
    beam1->setFExtPtr(weight);





    // --------------------
    // --- Interactions ---
    // --------------------

    // Interaction ball-floor
    // -- prismatic axis 0,0,1 in absolute frame: ball can only move in Z
    SP::SiconosVector axis1(new SiconosVector(3));
    axis1->setValue(0, 0);
    axis1->setValue(1, 0);
    axis1->setValue(2, 1);

    SP::PrismaticJointR relation1(new PrismaticJointR(axis1, true, beam1));

    SP::SimpleMatrix H1(new SimpleMatrix(relation1->numberOfConstraints(), qDim));
    H1->zero();
    relation1->setJachq(H1);

    SP::NonSmoothLaw nslaw1(new EqualityConditionNSL(relation1->numberOfConstraints()));

    SP::Interaction inter1(new Interaction(nslaw1, relation1));
    // -------------
    // --- Model ---
    // -------------
    SP::Model bouncingBall(new Model(t0, T));
    bouncingBall->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam1);
    bouncingBall->nonSmoothDynamicalSystem()->link(inter1, beam1);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));
    //    s->setComputeResiduY(true);
    //  s->setUseRelativeConvergenceCriteron(false);

    // -- OneStepIntegrators --
    SP::MoreauJeanOSI OSI1(new MoreauJeanOSI(theta));
    s->insertIntegrator(OSI1);


    // -- OneStepNsProblem --
    SP::OneStepNSProblem osnspb(new Equality());
    s->insertNonSmoothProblem(osnspb);
    bouncingBall->setSimulation(s);
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize();
    int N = 2000; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector q1 = beam1->q();
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
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
    ioMatrix::write("PrismaticTest.dat", "ascii", dataPlot, "noDim");

    cout << "====> Comparison with a reference file ..." << endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();

    ioMatrix::read("PrismaticTest.ref", "ascii", dataPlotRef);
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
    cout << "Exception caught in BouncingBallTS.cpp" << endl;
  }

}

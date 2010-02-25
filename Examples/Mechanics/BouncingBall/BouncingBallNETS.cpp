/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
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

using namespace std;

int main(int argc, char* argv[])
{
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    unsigned int qDim = 7;           // degrees of freedom for the ball
    unsigned int nDim = 6;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 2.0;      // initial velocity for lowest bead.
    double omega_init = 1.0;      // initial velocity for lowest bead.
    double theta = 1.0;              // theta for Moreau integrator
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    DynamicalSystemsSet allDS; // the list of DS

    // -- Initial positions and velocities --
    SP::SimpleVector q0(new SimpleVector(qDim));
    SP::SimpleVector v0(new SimpleVector(nDim));
    SP::SimpleMatrix I(new SimpleMatrix(3, 3));
    I->eye();
    (*q0)(0) = position_init;
    /*initial quaternion equal to (1,0,0,0)*/
    (*q0)(3) = 1.0;

    (*v0)(0) = velocity_init;
    (*v0)(3) = omega_init;
    // -- The dynamical system --
    SP::NewtonEulerDS ball(new NewtonEulerDS(q0, v0, m, I));
    allDS.insert(ball);
    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof));
    (*weight)(0) = -m * g;
    ball->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    InteractionsSet allInteractions;

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    SP::SiconosMatrix H(new SimpleMatrix(1, qDim));
    (*H)(0, 0) = 1.0;
    //     vector<SP::SiconosMatrix> vecMatrix1;
    //     vecMatrix1.push_back(H);
    //     SP::BlockMatrix H_block(new BlockMatrix(vecMatrix1,1,1));

    //     SP::SiconosMatrix HT(new SimpleMatrix(1,nDim));
    //     vector<SP::SiconosMatrix> vecMatrix2;
    //     vecMatrix2.push_back(HT);
    //     SP::BlockMatrix HT_block(new BlockMatrix(vecMatrix2,1,1));


    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
    SP::NewtonEulerR relation0(new NewtonEulerR());
    relation0->setJachq(H);
    //    relation0->setJacQH(H_block);
    //    relation0->setJacQHT(HT_block);
    cout << "main jacQH" << endl;
    relation0->jachq()->display();
    SP::Interaction inter(new Interaction("floor-ball", allDS, 0, 1, nslaw0, relation0));
    allInteractions.insert(inter);
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
    SP::Moreau OSI(new Moreau(ball, theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI, osnspb));


    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize(s);
    int N = (int)((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 6;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(2);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0, 5) = acos((*q)(3));
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (s->nextTime() < T)
    {
      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k, 5) = (*q)(3);

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

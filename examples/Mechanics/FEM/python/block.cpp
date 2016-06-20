
#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {


    // User-defined main parameters
    unsigned int ndof = 24;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10;                  // final computation time
    double h = 0.005;                // time step
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity

    cout << "====> Model loading ..." << endl << endl;

    SP::SiconosMatrix Mass(new SimpleMatrix(ndof, ndof));
    for (int i = 0; i < ndof; ++i)
      (*Mass)(i, i) = m;

    SP::SiconosMatrix K(new SimpleMatrix(*Mass));
    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(ndof));
    SP::SiconosVector v0(new SiconosVector(ndof));
    //(*q0)(0) = position_init;
    //(*v0)(0) = velocity_init;

    // -- The dynamical system --
    SP::LagrangianLinearTIDS ball(new LagrangianLinearTIDS(q0, v0, Mass));


    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(ndof));
    (*weight)(0) = -m * g;
    for (int i = 2; i < ndof; i = i + 3)
      (*weight)(i) = -m * g;

    ball->setFExtPtr(weight);
    ball->setKPtr(K);

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //

    int diminter = 4;

    SP::SiconosMatrix H(new SimpleMatrix(4, ndof));
    (*H)(0, 2) = 1.0;
    (*H)(1, 5) = 1.0;
    (*H)(2, 8) = 1.0;
    (*H)(3, 11) = 1.0;
    SP::SiconosVector b(new SiconosVector(4));
    (*b)(0) = 3.0;
    (*b)(1) = 3.0;
    (*b)(2) = 3.0;
    (*b)(3) = 3.0;



    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H, b));

    SP::Interaction inter(new Interaction(diminter, nslaw, relation));

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
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(ball, theta));

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
    bouncingBall->initialize();

    inter->y(0)->display();

    return 0;
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (s->hasNextEvent())
    {
      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);
      s->nextStep();
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");



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

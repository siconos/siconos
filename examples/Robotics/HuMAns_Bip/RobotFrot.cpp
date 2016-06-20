
// =============================== Robot sample (Bip) ===============================
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 21;           // degrees of freedom for robot
    double t0 = 0;                   // initial computation time
    double T = 1.3;//0.005;                   // final computation time
    double h = 0.001;                // time step
    double criterion = 0.5;
    unsigned int maxIter = 5000;
    double en = 0.0, et = 0.0, mu = 0.1;      // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;

    // --- DS: robot bip ---

    // The dof are angles between differents parts of the robot.

    // Initial position (angles in radian)
    SiconosVector q0(nDof), v0(nDof);
    q0(1) = -0.1;
    q0(2) = 0.2;
    q0(3) = -0.1;
    q0(5) = -0.1;
    q0(6) = 0.2;
    q0(7) = -0.1;

    SP::LagrangianDS bip(new LagrangianDS(q0, v0));

    // external plug-in
    bip->setComputeMassFunction("RobotFrotPlugin", "mass");
    bip->setComputeFGyrFunction("RobotFrotPlugin", "FGyr");
    bip->setComputeJacobianFGyrFunction(0, "RobotFrotPlugin", "jacobianFGyrq");
    bip->setComputeJacobianFGyrFunction(1, "RobotFrotPlugin", "jacobianVFGyr");

    // -------------------
    // --- Interactions---
    // -------------------

    // Two interactions:
    //  - one with Lagrangian non linear relation to define contact with ground
    //  - the other to define angles limitations (articular stops), with lagrangian linear relation
    //  Both with newton impact nslaw.

    // -- relations --

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(en, et, mu, 3));
    string G = "RobotFrotPlugin:G0";
    SP::Relation relation(new LagrangianScleronomousR("RobotFrotPlugin:h0", G));
    SP::Interaction inter(new Interaction(69, nslaw, relation));

    //The linear contraint corresponding to joints limits (hq+b>0)
    SP::NonSmoothLaw nslaw2(new NewtonImpactFrictionNSL(en, et, 0, 3));
    SimpleMatrix H(90, 21);
    SiconosVector b(90);
    H.zero();
    b.zero();
    for (unsigned int i = 0; i < 15; ++i)
    {
      for (unsigned int j = 0; j < 3; ++j)
      {
        H(6 * i + j, i) = -1;
        H(6 * i + j + 3, i) = 1;
      }
    }

    unsigned int i = 0;
    b(3 * i + 0) = 0.21;
    b(3 * i + 1) = 0.21;
    i++;
    b(3 * i + 2) = 0.64;
    b(3 * i + 3) = 0.64;
    i++;
    b(3 * i + 4) = 1.41;
    b(3 * i + 5) = -0.11;
    i++;
    b(3 * i + 6) = 0.20;
    b(3 * i + 7) = 1.17;
    i++;
    b(3 * i + 8) = 0.21;
    b(3 * i + 9) = 0.21;
    i++;
    b(3 * i + 10) = 0.64;
    b(3 * i + 11) = 0.64;
    i++;
    b(3 * i + 12) = 1.41;
    b(3 * i + 13) = -0.11;
    i++;
    b(3 * i + 14) = 0.2;
    b(3 * i + 15) = 1.17;
    i++;
    b(3 * i + 16) = 0.17;
    b(3 * i + 17) = 0.17;
    i++;
    b(3 * i + 18) = 0.14;
    b(3 * i + 19) = 0.14;
    i++;
    b(3 * i + 20) = 0.17;
    b(3 * i + 21) = 0.17;
    i++;
    b(3 * i + 22) = 0.14;
    b(3 * i + 23) = 0.14;
    i++;
    b(3 * i + 24) = 0.17;
    b(3 * i + 25) = 0.17;
    i++;
    b(3 * i + 26) = 0.08;
    b(3 * i + 27) = 0.08;
    i++;
    b(3 * i + 28) = 0.08;
    b(3 * i + 29) = 0.21;


    SP::Relation relation2(new LagrangianLinearTIR(H, b));
    SP::Interaction inter2(new Interaction(90, nslaw2, relation2));

    // -------------
    // --- Model ---
    // -------------

    SP::Model Robot(new Model(t0, T));
    Robot->nonSmoothDynamicalSystem()->insertDynamicalSystem(bip);
    Robot->nonSmoothDynamicalSystem()->link(inter1, bip);
    Robot->nonSmoothDynamicalSystem()->link(inter2, bip);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator OSI(new MoreauJeanOSI(bip, 0.500001));
    s->insertIntegrator(OSI);

    // -- OneStepNsProblem --
    string solverName = "NSGS";      // solver algorithm used for non-smooth problem
    IntParameters iparam(5);
    iparam[0] = 10100; // Max number of iteration
    // Solver/formulation
    // 0: projection, 1: Newton/AlartCurnier, 2: Newton/Fischer-Burmeister, 3: Path/Glocker
    iparam[4] = 4;
    DoubleParameters dparam(5);
    dparam[0] = 1e-6; // Tolerance
    dparam[2] = 1e-8; // Local Tolerance
    SP::NonSmoothSolver Mysolver(new NonSmoothSolver(solverName, iparam, dparam));
    SP::NonSmoothSolver mySolver(new NonSmoothSolver(solverName, iparam, dparam));
    SP::FrictionContact osnspb(new FrictionContact(3, mySolver));
    s->insertNonSmoothProblem(osnspb);
    cout << "=== End of model loading === " << endl;
    Robot->setSimulation(s);
    // =========================== End of model definition ===========================


    // ================================= Computation =================================

    // --- Simulation initialization ---
    Robot->initialize();
    cout << "End of model initialisation" << endl;

    int k = 0; // Current step
    unsigned int N = (unsigned int)((T - t0) / h);

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 22;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    dataPlot(k, 0) = Robot->t0();

    for (unsigned int i = 1; i < 22; i++)
      dataPlot(k, i) = (*bip->q())(i - 1);

    // --- Compute elapsed time ---
    boost::timer tt;
    tt.restart();

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (s->hasNextEvent())
    {
      // get current time step
      k++;
      cout << k << endl;
      s->newtonSolve(criterion, maxIter);
      s->nextStep();
      dataPlot(k, 0) = s->startingTime();

      for (unsigned int i = 1; i < 22; i++)
        dataPlot(k, i) = (*bip->q())(i - 1);
    }
    cout << "time = " << tt.elapsed() << endl;
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
}

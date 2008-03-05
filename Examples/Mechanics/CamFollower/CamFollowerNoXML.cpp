/* Siconos-Examples version 3.0.0, Copyright INRIA 2005-2008.
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

/*!\file CamFollowerNoXML.cpp
\brief \ref EMCamFollower - C++ input file version - M. di Bernardo, G. Osorio, S. Santini.
*/

#include "SiconosKernel.h"
#include "CamState.h"

using namespace std;

int main(int argc, char* argv[])
{
  double rpm = 358;
  try
  {

    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int dsNumber = 1;       // the Follower and the ground
    unsigned int nDof = 1;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 1;                   // final computation time
    double h = 0.0001;                // time step
    double position_init = 0.;//;40;      // initial position for lowest bead.
    double velocity_init = 0.;//4;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for Moreau integrator
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    SiconosMatrix *Mass, *K, *C;        // mass/rigidity/viscosity
    Mass = new SimpleMatrix(nDof, nDof);
    //    for(i=0;i<nDof;i++)
    (*Mass)(0, 0) = 1.221;
    K = new SimpleMatrix(nDof, nDof);
    (*K)(0, 0) = 1430.8;
    C = new SimpleMatrix(nDof, nDof);
    (*C)(0, 0) = 0;

    // -- Initial positions and velocities --
    vector<SimpleVector *> q0;
    vector<SimpleVector *> velocity0;
    q0.resize(dsNumber, NULL);
    velocity0.resize(dsNumber, NULL);
    q0[0] = new SimpleVector(nDof);
    velocity0[0] = new SimpleVector(nDof);
    (*(q0[0]))(0) = position_init;
    (*(velocity0[0]))(0) = velocity_init;
    LagrangianLinearTIDS * lds = new LagrangianLinearTIDS(0, *(q0[0]), *(velocity0[0]), *Mass, *K, *C);
    lds->setComputeFExtFunction("FollowerPlugin.so", "FollowerFExtR");

    // Example to set a list of parameters in FExt function.
    // 1 - Create a simple vector that contains the required parameters.
    SimpleVector * param = new SimpleVector(1); // Here we only set one parameter, the DS number.
    //    (*param)(0) = vectorDS[0]->getNumber();
    (*param)(0) = rpm;
    // 2 - Assign this param to the function FExt
    lds->setZPtr(param);
    // 2 corresponds to the position of FExt in the stl vector of possible parameters. 0 is mass, 1 FInt and so on.
    // Now the DS number will be available in FExt plugin.

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.8;

    // Interaction Follower-floor
    //
    SiconosMatrix *H = new SimpleMatrix(1, nDof);
    (*H)(0, 0) = 1.0;
    NonSmoothLaw * nslaw0 = new NewtonImpactNSL(e);
    Relation * relation0 = new LagrangianLinearR(*H);
    DynamicalSystemsSet dsConcerned;
    dsConcerned.insert(lds);

    Interaction * inter =  new Interaction("Follower-Ground", dsConcerned, 0, 1, nslaw0, relation0);
    InteractionsSet allInteractions;
    allInteractions.insert(inter);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(dsConcerned, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * Follower = new Model(t0, T);
    Follower->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, Follower);

    TimeStepping* S = new TimeStepping(t);


    // -- OneStepIntegrator --
    OneStepIntegrator*  OSI = new Moreau(lds, theta, S);

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 101; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-5; // Tolerance
    string solverName = "QP" ;
    NonSmoothSolver * mySolver = new NonSmoothSolver(solverName, iparam, dparam);
    OneStepNSProblem * osnspb = new LCP(S, mySolver);

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    S->initialize();
    cout << "End of simulation initialisation" << endl;


    int k = 0;
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    SimpleMatrix DataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    DataPlot(k, 0) = k * t->getH();
    DataPlot(k, 1) = lds->getQ()(0);
    DataPlot(k, 2) = lds->getVelocity()(0);
    DataPlot(k, 3) = (Follower->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
    DataPlot(k, 4) = lds->getFExt()(0);

    // State of the Cam
    //    double rpm=358;
    double CamEqForce, CamPosition, CamVelocity, CamAcceleration;

    CamEqForce = CamState(k * t->getH(), rpm, CamPosition, CamVelocity, CamAcceleration);
    // Position of the Cam
    DataPlot(k, 5) = CamPosition;
    // Velocity of the Cam
    DataPlot(k, 6) = CamVelocity;
    // Acceleration of the Cam
    DataPlot(k, 7) = CamPosition + lds->getQ()(0);

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (k < N)
    {
      // get current time step
      k++;
      // solve ...
      S->computeOneStep();

      // --- Get values to be plotted ---

      DataPlot(k, 0) = k * t->getH();
      DataPlot(k, 1) = lds->getQ()(0);
      DataPlot(k, 2) = lds->getVelocity()(0);
      DataPlot(k, 3) = (Follower->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      DataPlot(k, 4) = lds->getFExt()(0);

      CamEqForce = CamState(k * t->getH(), rpm, CamPosition, CamVelocity, CamAcceleration);
      DataPlot(k, 5) = CamPosition;
      DataPlot(k, 6) = CamVelocity;
      DataPlot(k, 7) = CamPosition + lds->getQ()(0);
      // transfer of state i+1 into state i and time incrementation
      S->nextStep();
    }
    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(DataPlot, "noDim");

    // --- Free memory ---
    delete osnspb;
    delete t;
    delete S;
    delete OSI;
    delete Follower;
    delete nsds;
    delete relation0;
    delete nslaw0;
    delete H;
    delete inter;
    delete lds;
    delete q0[0];
    delete velocity0[0];
    delete C;
    delete K;
    delete Mass;

    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;
    cout << "End of computation - Number of iterations done: " << k << endl;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/CamFollower\'" << endl;
  }
}

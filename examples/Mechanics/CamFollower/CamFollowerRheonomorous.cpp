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

// =============================== Cam Follower (1DOF Impact System) ===============================
//
// The Cam Follower system is modelled as a Generalised Langrangian System impacting against a fixed wall
// the moving constraint (i.e. a rotational cam) is modelled as an input force
//
// Direct description of the model.
//
// Keywords: LagrangianLinearDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// ======================================================================================================

#include "SiconosKernel.hpp"
#include "CamState.h"

using namespace std;

int main(int argc, char* argv[])
{
  double rpm = 358;
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int dsNumber = 1;       // the Follower and the ground
    unsigned int nDof = 1;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 1;                   // final computation time
    double h = 0.0001;                // time step
    double position_init = 0.40;      // initial position for lowest bead.
    double velocity_init = 0.4;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    SP::SimpleMatrix Mass(new SimpleMatrix(nDof, nDof));
    SP::SimpleMatrix K(new SimpleMatrix(nDof, nDof));
    SP::SimpleMatrix C(new SimpleMatrix(nDof, nDof));       // mass/rigidity/viscosity
    (*Mass)(0, 0) = 1.221;
    (*K)(0, 0) = 1430.8;

    // -- Initial positions and velocities --
    vector<SP::SiconosVector> q0;
    vector<SP::SiconosVector> velocity0;
    q0.resize(dsNumber);
    velocity0.resize(dsNumber);
    q0[0].reset(new SiconosVector(nDof));
    velocity0[0].reset(new SiconosVector(nDof));
    (*(q0[0]))(0) = position_init;
    (*(velocity0[0]))(0) = velocity_init;
    SP::LagrangianLinearTIDS lds(new LagrangianLinearTIDS(q0[0], velocity0[0], Mass, K, C));
    lds->setComputeFExtFunction("FollowerPlugin", "FollowerFExtR");

    // Example to set a list of parameters in FExt function.
    // 1 - Create a simple vector that contains the required parameters.
    SP::SiconosVector param(new SiconosVector(1)); // Here we only set one parameter, the DS number.
    //    (*param)(0) = vectorDS[0]->getNumber();
    (*param)(0) = rpm;
    // 2 - Assign this param to the function FExt
    lds->setzPtr(param);
    // 2 corresponds to the position of FExt in the stl vector of possible parameters. 0 is mass, 1 FInt and so on.
    // Now the DS number will be available in FExt plugin.
    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.8;

    // Interaction Follower-floor
    //
    SimpleMatrix H(1, nDof);
    H(0, 0) = 1.0;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
    vector<string> listofG2(1);
    listofG2[0] = "FollowerPlugin:FollowerComputeG0";
    SP::Relation relation0(new LagrangianRheonomousR("FollowerPlugin:FollowerComputeH1", "FollowerPlugin:FollowerComputeG10", "FollowerPlugin:FollowerComputeG11"));

    SP::SiconosVector param2(new SiconosVector(1)); // Here we only set one parameter, the DS number.
    (*param2)(0) = rpm;

    SP::Interaction inter(new Interaction(nslaw0, relation0));
    // -------------
    // --- Model ---
    // -------------

    SP::Model Follower(new Model(t0, T));
    Follower->nonSmoothDynamicalSystem()->insertDynamicalSystem(lds);
    Follower->nonSmoothDynamicalSystem()->link(inter,lds);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping S(new TimeStepping(t));


    // -- OneStepIntegrator --
    SP::OneStepIntegrator OSI(new MoreauJeanOSI(theta));
    S->insertIntegrator(OSI);

    // -- OneStepNsProblem --

    SP::OneStepNSProblem osnspb(new LCP(SICONOS_LCP_QP));

    // solver
    //osnspb->numericsSolverOptions()->solverId=SICONOS_LCP_QP;

    // max number of iterations
    osnspb->numericsSolverOptions()->iparam[0] = 101;

    // tolerance
    osnspb->numericsSolverOptions()->dparam[0] = 1e-6;

    S->insertNonSmoothProblem(osnspb);
    Follower->setSimulation(S);

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    Follower->initialize();
    cout << "End of model initialisation" << endl;


    int k = 0;
    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 8;
    SimpleMatrix DataPlot(N + 1, outputSize);
    // For the initial time step:
    // time
    DataPlot(k, 0) = t0;
    DataPlot(k, 1) = (*lds->q())(0);
    DataPlot(k, 2) = (*lds->velocity())(0);
    DataPlot(k, 3) = (*inter->lambda(1))(0);
    DataPlot(k, 4) = (*lds->fExt())(0);

    // State of the Cam
    //    double rpm=358;
    double CamEqForce, CamPosition, CamVelocity, CamAcceleration;

    CamEqForce = CamState(t0, rpm, CamPosition, CamVelocity, CamAcceleration);
    // Position of the Cam
    DataPlot(k, 5) = CamPosition;
    // Velocity of the Cam
    DataPlot(k, 6) = CamVelocity;
    // Acceleration of the Cam
    DataPlot(k, 7) = CamPosition + (*lds->q())(0);
    boost::timer tt;
    tt.restart();
    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (k < N)
    {
      // get current time step
      k++;
      // solve ...
      S->computeOneStep();

      // --- Get values to be plotted ---

      DataPlot(k, 0) = S->nextTime();
      DataPlot(k, 1) = (*lds->q())(0);
      DataPlot(k, 2) = (*lds->velocity())(0);
      DataPlot(k, 3) = (*inter->lambda(1))(0);
      DataPlot(k, 4) = (*lds->fExt())(0);

      CamEqForce = CamState(S->nextTime(), rpm, CamPosition, CamVelocity, CamAcceleration);
      DataPlot(k, 5) = CamPosition;
      DataPlot(k, 6) = CamVelocity;
      DataPlot(k, 7) = CamPosition + (*lds->q())(0);
      // transfer of state i+1 into state i and time incrementation
      S->nextStep();
    }
    // --- Output files ---
    ioMatrix::write("result.dat", "ascii", DataPlot, "noDim");
    cout << "time = " << tt.elapsed() << endl;
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

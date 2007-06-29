/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2006.
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
 *
 * Multi-beads 3D frictionl contact problem in presence of a rigid foundations
 * 14/05/2007- Authors: houari khenous

*/
// =============================== Multi bouncing beads couple simulation ===============================
//
// Keywords: LagrangianLinearDS, LagrangianDS relation, Moreau TimeStepping, newton method.
//
// ======================================================================================================


#include "SiconosKernel.h"

using namespace std;


int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters

    unsigned int COUPLE = 1;       // the number of dynamical systems

    unsigned int nDof = 6;            // degrees of freedom for beads

    double t0 = 0;                    // initial computation time
    double T = 10.;                    // final computation time
    double h = 0.005;                 // time step
    double PI = 3.14;

    string solverName = "NEWTONFUNCTION";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGSNEWTON";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    //string solverName = "Lemke";      // solver algorithm used for non-smooth problem

    double e  = 0.9;                  // nslaw
    double mu = 0.1;


    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    unsigned int i;

    DynamicalSystemsSet allDS; // the list of DS
    LagrangianDS *GLOB_tabLDS; // table of Lagrangian DS


    // -- Initial positions and velocities --

    vector<SimpleVector *> q0;
    vector<SimpleVector *> v0;
    q0.resize(COUPLE, NULL);
    v0.resize(COUPLE, NULL);

    // Memory allocation for q0[i] and v0[i]
    for (i = 0; i < COUPLE; i++)
    {
      q0[i] = new SimpleVector(nDof);
      v0[i] = new SimpleVector(nDof);
    }

    // set values

    // (*(q0[0]))(0) =  0.0;  (*(q0[0]))(1) =  0.;  (*(q0[0]))(2) =  0.3; (*(q0[0]))(3) =  PI/2;

    // multi couple without contact between coulpes

    for (i = 0; i < COUPLE; i++)
    {
      (*(q0[i]))(0) = 0.;
      (*(q0[i]))(1) = 0.5 * (i + 1.);
      (*(q0[i]))(2) =  0.3;
      (*(q0[i]))(3) =  PI / 3;
      //  (*(q0[i]))(1) = 0.5*i;
    }

    GLOB_tabLDS = new LagrangianDS(0, *(q0[0]), *(v0[0]));

    GLOB_tabLDS->setComputeFExtFunction("NLMPlugin.so", "gravity");
    GLOB_tabLDS->setComputeMassFunction("NLMPlugin.so", "mass");
    GLOB_tabLDS->setComputeNNLFunction("NLMPlugin.so", "NNL");
    GLOB_tabLDS->setComputeJacobianNNLFunction(0, "NLMPlugin.so", "jacobianQNNL");
    GLOB_tabLDS->setComputeJacobianNNLFunction(1, "NLMPlugin.so", "jacobianVNNL");

    allDS.insert(GLOB_tabLDS);

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;
    DynamicalSystemsSet dsConcernedi;

    NonSmoothLaw * nslaw1 = new NewtonImpactFrictionNSL(e, e, mu, 3);


    Relation * relation1 = new LagrangianScleronomousR("NLMPlugin:h1", "NLMPlugin:G1");
    Interaction * inter1 = new Interaction("bead1", allDS, 0, 3, nslaw1, relation1);

    Relation * relation2 = new LagrangianScleronomousR("NLMPlugin:h2", "NLMPlugin:G2");
    Interaction * inter2 = new Interaction("bead2", allDS, 0, 3, nslaw1, relation2);

    allInteractions.insert(inter1);
    allInteractions.insert(inter2);


    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * BeadsCOUPLE = new Model(t0, T);
    BeadsCOUPLE->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    TimeDiscretisation * GLOB_T = new TimeDiscretisation(h, BeadsCOUPLE);
    TimeStepping* GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --

    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 100, 0.001);

    // OneStepNSProblem * osnspb = new LCP(GLOB_SIM,"name",solverName,100, 0.001,1);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition
    // ================================= Computation

    // ================================= Computation =================================

    // --- Simulation initialization ---

    GLOB_SIM->initialize();
    cout << "End of simulation initialisation" << endl;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot

    int k = 0; // index for output.
    int N = GLOB_T->getNSteps(); // Number of time steps

    unsigned int outputSize = 1 + 2 * COUPLE;
    SimpleMatrix dataPlot(N + 1, outputSize);

    dataPlot(k, 0) = k * GLOB_T->getH();
    dataPlot(k, 1) = (BeadsCOUPLE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getY(0))(0); //GLOB_tabLDS->getQ()(2);
    dataPlot(k, 2) = GLOB_tabLDS->getVelocity()(2);

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    EventsManager * eventsManager = GLOB_SIM->getEventsManagerPtr();

    while (eventsManager->hasNextEvent())
    {
      GLOB_SIM->computeOneStep();
      GLOB_SIM->advanceToEvent();
      GLOB_SIM->processEvents();
      // --- Get values to be plotted ---
      k++;
      dataPlot(k, 0) = k * GLOB_T->getH();
      dataPlot(k, 1) = (BeadsCOUPLE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getY(0))(0); //GLOB_tabLDS->getQ()(2);
      dataPlot(k, 2) = GLOB_tabLDS->getVelocity()(2);

      GLOB_SIM->nextStep();
    }
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    delete OSI;
    delete osnspb;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/NLM Init\'" << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl << endl;
}






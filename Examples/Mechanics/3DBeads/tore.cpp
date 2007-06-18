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
 * Torus 3D frictionl contact problem in presence of a rigid foundations
 * 30/01/2007- Authors: houari khenous & Roger Pissard

*/
// =============================== Torus simulation ===============================
//
// Keywords: LagrangianLinearTIDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
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

    //    unsigned int DSNUMBER = 1;       // the number of dynamical systems

    unsigned int FEM = 48;

    // unsigned int nDof = 3;            // degrees of freedom for beads

    double t0 = 0;                    // initial computation time
    double T = 1;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "NLGSNEWTON";      // solver algorithm used for non-smooth problem

    double e = 0.8;                  // nslaw
    double mu = 10.;


    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    unsigned int i;
    unsigned int j;
    unsigned int k;

    DynamicalSystemsSet allDS; // the list of DS
    LagrangianLinearTIDS *GLOB_tabLDS; // table of Lagrangian DS


    SiconosMatrix *M;
    SiconosMatrix *K;
    SiconosMatrix *C;
    SiconosMatrix *Position;

    // mass matrix
    M = new SimpleMatrix("mass", 1);

    // rigid matrix
    K = new SimpleMatrix("rigid", 1);

    // maillage matrix
    Position = new SimpleMatrix("position", 1);


    // amortissement matrix
    C = new SimpleMatrix(FEM, FEM);

    // -- Initial positions and velocities --
    // q0 and v0.

    SimpleVector* q0 = new SimpleVector(FEM);
    SimpleVector* v0 = new SimpleVector(FEM);

    // Memory allocation for q0[i] and v0[i]
    for (i = 0; i < FEM; i++)
      (*q0)(i) = 1.;

    // set values

    GLOB_tabLDS = new LagrangianLinearTIDS(0, *q0, *v0, *M, *K, *C);


    allDS.insert(GLOB_tabLDS);



    // ==> at this point, all the required dynamical systems are saved in allDS.


    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;

    cp = []

         SiconosMatrix * H1 = new SimpleMatrix(cp, FEM);
    for (i = 0; i < FEM; i++)
    {
      (*H1)(i, 0) = 1.;
      (*H1)(i, 1) = 1.;
      (*H1)(i, 2) = 1.;
    }

    SimpleVector *b1 = new SimpleVector(10);
    j = 0.;
    for (dal::bv_visitor i(cn); !i.finished(); ++i)
      if (i % 3 == 0)
      {
        (*H1)(j, i + 2) = 1.;
        (*H1)(2 * j, i) = 1.;
        (*H1)(2 * j + 1, i + 1) = 1.;
        ++j;
        (*b1)(j) = (*Position)(i, j);
      }


    //  vector<SimpleVector *> Pos;
    //     Pos.resize(FEM,NULL);

    //     for (i=0;i<FEM;i++){
    //       Pos[i] = new SimpleVector(3);
    //       for (j=0;j<3;j++)
    //  (*(Pos[i]))(j) = ;
    //}



    for (i = 0; i < FEM; i++)
    {
      (*b1)(i) = 0.;
      (*b1)(i + 1) = 0.;
      (*b1)(i + 2) = -10.;
    }
    NonSmoothLaw* nslaw1 = new NewtonImpactFrictionNSL(e, e, mu, 3);
    Relation* relation1 = new LagrangianLinearR(*H1, *b1);
    Interaction * inter1 = new Interaction("bead1", allDS, 0, 3, nslaw1, relation1);


    allInteractions.insert(inter1);

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * TORE = new Model(t0, T);
    TORE->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    TimeDiscretisation * GLOB_T = new TimeDiscretisation(h, TORE);
    TimeStepping* GLOB_SIM = new TimeStepping(GLOB_T);

    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , GLOB_SIM);

    // -- OneStepNsProblem --
    //OneStepNSProblem * osnspb = new LCP(GLOB_SIM,"FrictionContact3D",solverName,101,0.001);

    OneStepNSProblem * osnspb = new FrictionContact3D(GLOB_SIM , "FrictionContact3D", solverName, 100, 0.001);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition
    // ================================= Computation

    // ================================= Computation =================================

    // --- Simulation initialization ---

    GLOB_SIM->initialize();
    cout << "End of simulation initialisation" << endl;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot

    //  unsigned int outputSize = 1+3*DSNUMBER;
    //     SimpleMatrix dataPlot(N+1,outputSize);
    int k_iter = 0; // index for output.
    //     dataPlot(k_iter,0) = k_iter*GLOB_T->getH();
    //     dataPlot(k_iter,1) = GLOB_tabLDS[0]->getQ()(2);
    //     dataPlot(k_iter,2) = GLOB_tabLDS[0]->getVelocity()(2);
    //     dataPlot(k_iter,3) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);

    //   dataPlot(k_iter,4) = GLOB_tabLDS[1]->getQ()(2);
    //    dataPlot(k_iter,5) = GLOB_tabLDS[1]->getVelocity()(2);
    //    dataPlot(k_iter,6) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    EventsManager * eventsManager = GLOB_SIM->getEventsManagerPtr();

    while (eventsManager->hasNextEvent())
    {
      GLOB_SIM->computeOneStep();
      //  GLOB_SIM->advanceToEvent();
      //  GLOB_SIM->processEvents();
      // --- Get values to be plotted ---
      k_iter++;
      //  dataPlot(k_iter,0) = k_iter*GLOB_T->getH();
      //  dataPlot(k_iter,1) = GLOB_tabLDS[0]->getQ()(2);
      //  dataPlot(k_iter,2) = GLOB_tabLDS[0]->getVelocity()(2);
      //  dataPlot(k_iter,3) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      //  dataPlot(k_iter,4) = GLOB_tabLDS[1]->getQ()(2);
      //  dataPlot(k_iter,5) = GLOB_tabLDS[1]->getVelocity()(2);
      //  dataPlot(k_iter,6) = (TORE->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);
      GLOB_SIM->nextStep();
    }
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    //    ioMatrix io("result.dat", "ascii");
    //     io.write(dataPlot,"noDim");
    //    cout<<"End of computation - Number of iterations done: "<<k<<endl;

    delete OSI;
    delete osnspb;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/TORE Init\'" << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl << endl;
}






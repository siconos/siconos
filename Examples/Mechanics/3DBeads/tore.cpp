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
 * 21/06/2007- Authors: houari khenous

*/
// =============================== Torus simulation ===============================
//
// Keywords: LagrangianLinearTIDS relation, Moreau TimeStepping, NLGS, NLGSNEWTON.
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

    unsigned int FEM = 324;

    // unsigned int nDof = 3;            // degrees of freedom for beads

    double t0 = 0;                    // initial computation time
    double T = 1;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "NLGSNEWTON";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem

    double e = 0.8;                  // nslaw
    double mu = 10.;
    double rho = 6e-6;
    double g = 9.81;

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

    double gap = 1.;

    // Memory allocation for q0[i] and v0[i]
    for (i = 1; i < FEM; i++)
    {
      if (i % 3 == 0)
      {
        (*q0)(i - 1) = (*Position)(i - 1, 0) + gap;
      }
    }

    GLOB_tabLDS = new LagrangianLinearTIDS(0, *q0, *v0, *M, *K, *C);


    allDS.insert(GLOB_tabLDS);


    // -- Set external forces (weight) --
    SiconosVector * weight = new SimpleVector(FEM);
    for (i = 1; i < FEM; i++)
      if (i % 3 == 0)
        (*weight)(i - 1) = -rho * g;
    GLOB_tabLDS->setFExtPtr(weight);

    // ==> at this point, all the required dynamical systems are saved in allDS.


    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;
    int cp = 72;

    //cp = [0 1 2 6 7 8 12 13 14 18 19 20 24 25 26 30 31 32 54 55 56 57 58 59 66 67 68 72 73 74 75 76 77 78 79 80 81 82 83 87 88 89 99 100 101 105 106 107 198 199 200 201 202 203 204 205 206 213 214 215 288 289 290 291 292 293 294 295 296 303 304 305]

    SiconosMatrix *H = new SimpleMatrix(cp, FEM);
    SimpleVector *b = new SimpleVector(cp);


    std::vector<int> v;
    v.resize(cp / 3);

    v[0] = 0;
    v[1] = 6;
    v[2] = 12;
    v[3] = 18;
    v[4] = 24;
    v[5] = 30;
    v[6] = 54;
    v[7] = 57;
    v[8] = 66;
    v[9] = 72;
    v[10] = 75;
    v[11] = 78;
    v[12] = 81;
    v[13] = 87;
    v[14] = 99;
    v[15] = 105;
    v[16] = 198;
    v[17] = 201;
    v[18] = 204;
    v[19] = 213;
    v[20] = 288;
    v[21] = 291;
    v[22] = 294;
    v[23] = 303;


    j = 0;
    for (size_t m = 0, size = v.size(); m < size; ++m)
    {
      (*H)(3 * j + 2, v[m] + 2) = 1.;
      (*H)(3 * j, v[m]) = 1.;
      (*H)(3 * j + 1, v[m] + 1) = 1.;
      (*b)(3 * j + 2) = (*Position)(v[m] + 2, 0) + gap;
      ++j;
    }

    NonSmoothLaw* nslaw = new NewtonImpactFrictionNSL(e, e, mu, 3);
    Relation* relation = new LagrangianLinearR(*H, *b);
    Interaction * inter = new Interaction("bead1", allDS, 0, cp, nslaw, relation);


    allInteractions.insert(inter);

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

    // --- Free memory ---
    delete osnspb;
    delete GLOB_T;
    delete GLOB_SIM;
    delete OSI;
    delete GLOB_tabLDS;
    delete TORE;
    delete nsds;
    delete inter;
    delete relation;
    delete nslaw;
    delete H;
    delete Position;
    delete b;
    delete q0;
    delete v0;
    delete K;
    delete M;
    delete C;


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






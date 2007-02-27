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
 * Two beads 3D frictionless contact problem in presence of a rigid foundation
 * 26/12/2006 - houari.khenous@inrialpes.fr

*/
// =============================== Multi bouncing beads column simulation ===============================
//  N beads between a floor and a ceiling ...
// Keywords: LagrangianLinearDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// ======================================================================================================
#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
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
    unsigned int dsNumber = 2;        // the number of dynamical systems
    unsigned int nDof = 6;            // degrees of freedom for beads

    double m1 = 1.;                   // mass of ball 1
    double R1 = 0.1;                   // radius of ball 1

    double m2 = 1.;                   // mass of ball 2
    double R2 = 0.1;                   // radius of ball 2

    double t0 = 0;                    // initial computation time
    double T = 10;                    // final computation time
    double h = 0.005;                 // time step
    string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    double criterion = 0.001;
    unsigned int maxIter = 5000;
    double e = 0.5;                  // nslaw
    double e2 = 0.9;                  // nslaw2
    double mu = 0.5;
    double mu2 = 0.5;

    // initial position ball 1
    double x1 = 0.;
    double yy1 = 0.;
    double z1 = 0.2;
    double theta1 = 0.;
    double phi1 = 0;
    double psi1 = 0;
    // initial velocity ball 1
    double dot_x1 = 0.;
    double dot_y1 = 0.;
    double dot_z1 = 0.;
    double dot_theta1 = 0.;
    double dot_phi1 = 0.;
    double dot_psi1 = 0.;
    // initial position ball 2
    double x2 = 0.;
    double y2 = 0.;
    double z2 = 0.5;
    double theta2 = 0.0;
    double phi2 = 0;
    double psi2 = 0;
    // initial velocity ball 2
    double dot_x2 = 0.;
    double dot_y2 = 0.;
    double dot_z2 = 0.;
    double dot_theta2 = 0.;
    double dot_phi2 = 0.;
    double dot_psi2 = 0.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    unsigned int i;
    DynamicalSystemsSet allDS; // the list of DS
    CheckInsertDS checkDS;


    // -- Initial positions and velocities --
    vector<SimpleVector *> q0;
    vector<SimpleVector *> velocity0;
    q0.resize(dsNumber, NULL);
    velocity0.resize(dsNumber, NULL);

    q0[0] = new SimpleVector(nDof);
    velocity0[0] = new SimpleVector(nDof);

    (*(q0[0]))(0) = x1;
    (*(velocity0[0]))(0) = dot_x1;
    (*(q0[0]))(1) = yy1;
    (*(velocity0[0]))(1) = dot_y1;
    (*(q0[0]))(2) = z1;
    (*(velocity0[0]))(2) = dot_z1;
    (*(q0[0]))(3) = theta1;
    (*(velocity0[0]))(3) = dot_theta1;
    (*(q0[0]))(4) = phi1;
    (*(velocity0[0]))(4) = dot_phi1;
    (*(q0[0]))(5) = psi1;
    (*(velocity0[0]))(5) = dot_psi1;

    q0[1] = new SimpleVector(nDof);
    velocity0[1] = new SimpleVector(nDof);

    (*(q0[1]))(0) = x2;
    (*(velocity0[1]))(0) = dot_x2;
    (*(q0[1]))(1) = y2;
    (*(velocity0[1]))(1) = dot_y2;
    (*(q0[1]))(2) = z2;
    (*(velocity0[1]))(2) = dot_z2;
    (*(q0[1]))(3) = theta2;
    (*(velocity0[1]))(3) = dot_theta2;
    (*(q0[1]))(4) = phi2;
    (*(velocity0[1]))(4) = dot_phi2;
    (*(q0[1]))(5) = psi2;
    (*(velocity0[1]))(5) = dot_psi2;

    LagrangianDS* lds1 = new LagrangianDS(0, *(q0[0]), *(velocity0[0]));
    LagrangianDS* lds2 = new LagrangianDS(1, *(q0[1]), *(velocity0[1]));

    // weight of beads as internal forces

    SimpleVector * poids1 = new SimpleVector(nDof);
    SimpleVector * poids2 = new SimpleVector(nDof);
    double g1 = 9.81;
    double g2 = 9.81;

    (*poids1)(2) =  -m1 * g1;
    (*poids2)(2) =  -m2 * g2;

    lds1->setFExtPtr(poids1);
    lds2->setFExtPtr(poids2);


    // external forces plug-in

    lds1->setComputeMassFunction("3DBeadsPlugin.so", "Mass1");
    lds2->setComputeMassFunction("3DBeadsPlugin.so", "Mass2");

    allDS.insert(lds1);
    allDS.insert(lds2);

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------

    InteractionsSet allInteractions;
    int interactionNumber = dsNumber + 1;
    vector<string> id;
    id.resize(interactionNumber - 2);
    unsigned int nInter = 3; // number of relations in each interaction
    DynamicalSystemsSet dsConcerned0;
    DynamicalSystemsSet dsConcerned1;

    // Interaction beads and floor
    SiconosVector *b0 = new SimpleVector(3);
    SiconosVector *b1 = new SimpleVector(3);
    (*b0)(0) = -R1;
    (*b1)(0) = -R2;
    SiconosMatrix *H1 = new SimpleMatrix(3, nDof);
    (*H1)(0, 2) = 1.0;
    (*H1)(1, 0) = 1.0;
    (*H1)(1, 4) = -R1;
    (*H1)(2, 1) = 1.0;
    (*H1)(2, 3) =  R1;
    SiconosMatrix *H2 = new SimpleMatrix(3, nDof);
    (*H2)(0, 2) = 1.0;
    (*H2)(1, 0) = 1.0;
    (*H2)(1, 4) = -R2;
    (*H2)(2, 1) = 1.0;
    (*H2)(2, 3) =  R2;


    NonSmoothLaw * nslaw = new NewtonImpactFrictionNSL(e, e, mu, 3);

    Relation * relation0 = new LagrangianLinearR(*H1, *b0);
    Relation * relation1 = new LagrangianLinearR(*H2, *b1);
    dsConcerned0.insert(lds1);
    Interaction * inter0 = new Interaction("floor_bead1", dsConcerned0, 0, 3, nslaw, relation0);
    dsConcerned1.insert(lds2);
    Interaction * inter1 = new Interaction("floor_bead2", dsConcerned1, 1, 3, nslaw, relation1);

    allInteractions.insert(inter0);
    allInteractions.insert(inter1);

    dsConcerned0.clear();
    dsConcerned1.clear();

    // Interaction between beads

    DynamicalSystemsSet dsConcerned2 ;
    CheckInsertInteraction checkInter;
    vector<Relation*> LLR(interactionNumber);

    NonSmoothLaw * nslaw2 = new NewtonImpactFrictionNSL(e2, e2, mu2, 3);// new NewtonImpactNSL(e2);

    Relation * relation2 = new LagrangianScleronomousR("3DBeadsPlugin:h0", "3DBeadsPlugin:G0");

    for (i = 1; (int)i < interactionNumber - 1; i++)
    {
      dsConcerned2.insert(allDS.getDynamicalSystemPtr(i - 1));
      dsConcerned2.insert(allDS.getDynamicalSystemPtr(i));
      ostringstream ostr;
      ostr << i;
      id[i - 1] = ostr.str();
      LLR[i - 1] = new LagrangianScleronomousR("3DBeadsPlugin:h0", "3DBeadsPlugin:G0");
      checkInter = allInteractions.insert(new Interaction(id[i - 1], dsConcerned2, i, nInter, nslaw2, LLR[i - 1]));
      dsConcerned2.clear();
    }



    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    bool isBVP = 0;
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP);

    // -------------
    // --- Model ---
    // -------------

    Model * multiBeads = new Model(t0, T);
    multiBeads->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    TimeDiscretisation * t = new TimeDiscretisation(h, multiBeads);
    TimeStepping* s = new TimeStepping(t);


    // -- OneStepIntegrators --
    OneStepIntegrator * OSI = new Moreau(allDS , 0.5000001 , s);

    // -- OneStepNsProblem --
    OneStepNSProblem * osnspb = new FrictionContact3D(s, "FrictionContact3D", solverName, 101, 0.001);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition
    // ================================= Computation
    // --- Simulation initialization ---
    s->initialize();
    cout << "End of simulation initialisation" << endl;

    int k;
    int N = t->getNSteps(); // Number of time steps

    // Prepare output and save value for the initial time
    unsigned int outputSize = dsNumber * 6 + 1 + 1;
    SimpleMatrix dataPlot(N + 1, outputSize); // Output data matrix

    cout << "Start computation ... " << endl;
    while (k < N)
    {

      // solve ...
      s->newtonSolve(criterion, maxIter);

      //  cout << " position" << endl;
      //  cout << (LLR[0])->getInteractionPtr()->getY(0)(0) << endl;
      //  cout << (lds1)->getQ()(2) << endl;
      //  cout << (lds2)->getQ()(2) << endl;
      k++;
    }
    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");


    // --- Free memory ---
    delete osnspb;
    delete t;
    delete s;
    delete multiBeads;
    delete nsds;
    delete relation0;
    delete relation2;
    delete nslaw;
    delete b0;
    delete H1;
    delete b1;
    delete H2;
    delete OSI;
    for (i = 0; i < dsNumber; i++)
    {
      delete q0[i];
      delete velocity0[i];
    }

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
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
}

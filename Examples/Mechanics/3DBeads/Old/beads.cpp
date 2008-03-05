/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2008.
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
 * 30/01/2007- Authors: houari khenous

*/
// =============================== Multi beads frictional contact simulation ============================

// Keywords: LagrangianLinearDS, LagrangianLinear relation, Moreau TimeStepping, Non-smooth Newton method.
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

    unsigned int DSNUMBER = 5;       // the number of dynamical systems

    unsigned int nDof = 6;            // degrees of freedom for beads

    double m = 1.;                   // mass of balls
    double R = 0.1;                   // radius of balls

    double t0 = 0;                    // initial computation time
    double T = 5;                    // final computation time
    double h = 0.005;                 // time step

    string solverName = "NSGS";      // solver algorithm used for non-smooth problem
    //string solverName = "NLGS";      // solver algorithm used for non-smooth problem
    double e = 0.8;                  // nslaw
    double e2 = 0.8;                  // nslaw2
    double mu = 0.2;


    // 1 to take in account the obstacle and  0 no
    double Top = 2.2;
    double Wall = 1.;
    double Ground = 0.;

    int obst_z_p = 0;                    //  for z --> +
    int obst_z_m = 1;                    //  for z --> -
    int obst_y_p = 1;                    //  for y --> +
    int obst_y_m = 1;                    //  for y --> -
    int obst_x_p = 1;                    //  for x --> +
    int obst_x_m = 1;                    //  for x --> -


    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    unsigned int Fact;
    Fact = (DSNUMBER) * (DSNUMBER - 1) / 2;

    unsigned int i;
    unsigned int j;
    unsigned int l;
    DynamicalSystemsSet allDS; // the list of DS
    LagrangianDS *GLOB_tabLDS[DSNUMBER]; // table of Lagrangian DS
    CheckInsertDS checkDS;

    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    (*Mass)(0, 0) = (*Mass)(1, 1) = (*Mass)(2, 2) = m;    ;
    (*Mass)(3, 3) = (*Mass)(4, 4) = (*Mass)(5, 5) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    // q0[i] and v0[i] correspond to position and velocity of ball i.

    vector<SimpleVector *> q0;
    vector<SimpleVector *> v0;
    q0.resize(DSNUMBER, NULL);
    v0.resize(DSNUMBER, NULL);

    // Memory allocation for q0[i] and v0[i]
    for (i = 0; i < DSNUMBER; i++)
    {
      q0[i] = new SimpleVector(nDof);
      v0[i] = new SimpleVector(nDof);
    }

    // set values

    (*(q0[0]))(0) =  0.0;
    (*(q0[0]))(1) =  0.3;
    (*(q0[0]))(2) =  0.12;
    (*(q0[1]))(0) =  0.2;
    (*(q0[1]))(1) =  0.3;
    (*(q0[1]))(2) =  0.12;


    //   (*(q0[2]))(0) =  0.;    (*(q0[2]))(1) =  0.5;  (*(q0[2]))(2) =  0.1;
    //     (*(q0[3]))(0) =  0.;    (*(q0[3]))(1) =  -0.5;  (*(q0[3]))(2) =  0.1;

    //  (*(q0[0]))(0) =  0.;    (*(q0[0]))(1) =  0.;  (*(q0[0]))(2) =  0.35;
    //     (*(q0[1]))(0) =  0.0;    (*(q0[1]))(1) =  0.;  (*(q0[1]))(2) =  0.12;



    // billard

    //   (*(q0[0]))(0) =  0.;     (*(q0[0]))(1) =  0.;   (*(q0[0]))(2) =  0.1;
    //     (*(q0[1]))(0) =  0.1;    (*(q0[1]))(1) = -0.2;  (*(q0[1]))(2) =  0.1;
    //     (*(q0[2]))(0) = -0.1;    (*(q0[2]))(1) = -0.2;  (*(q0[2]))(2) =  0.1;
    //     (*(q0[3]))(0) =  0.2;    (*(q0[3]))(1) = -0.4;  (*(q0[3]))(2) =  0.1;
    //     (*(q0[4]))(0) = -0.2;    (*(q0[4]))(1) = -0.4;  (*(q0[4]))(2) =  0.1;
    //     (*(q0[5]))(0) =  0.;     (*(q0[5]))(1) = -0.4;  (*(q0[5]))(2) =  0.1;
    //     (*(q0[6]))(0) =  0.1;    (*(q0[6]))(1) = -0.6;  (*(q0[6]))(2) =  0.1;
    //     (*(q0[7]))(0) = -0.1;    (*(q0[7]))(1) = -0.6;  (*(q0[7]))(2) =  0.1;

    //     (*(q0[8]))(0) =  0.;     (*(q0[8]))(1) =  0.8;  (*(q0[8]))(2) =  0.1;

    //     (*(v0[8]))(0) = -1;      (*(v0[8]))(1) = -20.;

    //     (*(q0[9]))(0) =  0.3;    (*(q0[9]))(1) = -0.6;  (*(q0[9]))(2) =  0.1;
    //     (*(q0[10]))(0)= -0.3;    (*(q0[10]))(1)= -0.6;  (*(q0[10]))(2)=  0.1;
    //     (*(q0[11]))(0)=  0.2;    (*(q0[11]))(1)= -0.8;  (*(q0[11]))(2)=  0.1;
    //     (*(q0[12]))(0)= -0.2;    (*(q0[12]))(1)= -0.8;  (*(q0[12]))(2)=  0.1;
    //     (*(q0[13]))(0)=  0.;     (*(q0[13]))(1)= -0.8;  (*(q0[13]))(2)=  0.1;
    //     (*(q0[14]))(0)=  0.4;    (*(q0[14]))(1)= -0.8;  (*(q0[14]))(2)=  0.1;
    //     (*(q0[15]))(0)= -0.4;    (*(q0[15]))(1)= -0.8;  (*(q0[15]))(2)=  0.1;

    // Cube de billes


    //  (*(q0[0]))(0) =  0.2;    (*(q0[0]))(1) = -0.2;  (*(q0[0]))(2) =  0.2;
    //     (*(q0[1]))(0) =  0.2;    (*(q0[1]))(1) =  0.2;  (*(q0[1]))(2) =  0.2;
    //     (*(q0[2]))(0) = -0.2;    (*(q0[2]))(1) =  0.2;  (*(q0[2]))(2) =  0.2;
    //     (*(q0[3]))(0) = -0.2;    (*(q0[3]))(1) = -0.2;  (*(q0[3]))(2) =  0.2;
    //     (*(q0[4]))(0) =  0.25;   (*(q0[4]))(1) = -0.2;  (*(q0[4]))(2) =  0.4;
    //     (*(q0[5]))(0) =  0.25;   (*(q0[5]))(1) =  0.2;  (*(q0[5]))(2) =  0.4;
    //     (*(q0[6]))(0) = -0.25;   (*(q0[6]))(1) =  0.2;  (*(q0[6]))(2) =  0.4;
    //     (*(q0[7]))(0) = -0.25;   (*(q0[7]))(1) = -0.2;  (*(q0[7]))(2) =  0.4;

    //     (*(q0[8]))(0) =  0.;     (*(q0[8]))(1) = 0.3;    (*(q0[8]))(2) =  0.1;

    //     (*(v0[8]))(1) =  -10;

    //     (*(q0[9]))(0) =  0.2;    (*(q0[9]))(1) =  0.2;  (*(q0[9]))(2) =  0.6;
    //     (*(q0[10]))(0)= -0.2;    (*(q0[10]))(1)=  0.2;  (*(q0[10]))(2)=  0.6;
    //     (*(q0[11]))(0)= -0.2;    (*(q0[11]))(1)= -0.2;  (*(q0[11]))(2)=  0.6;
    //     (*(q0[12]))(0)=  0.25;   (*(q0[12]))(1)= -0.2;  (*(q0[12]))(2)=  0.8;
    //     (*(q0[13]))(0)=  0.25;   (*(q0[13]))(1)=  0.2;  (*(q0[13]))(2)=  0.8;
    //     (*(q0[14]))(0)= -0.25;   (*(q0[14]))(1)=  0.2;  (*(q0[14]))(2)=  0.9;
    //     (*(q0[15]))(0)= -0.25;   (*(q0[15]))(1)= -0.2;  (*(q0[15]))(2)=  0.9;
    //     (*(q0[16]))(0)=  0.35;   (*(q0[16]))(1)=-0.35;  (*(q0[16]))(2)=  0.1;
    //     (*(q0[17]))(0)=  0.35;   (*(q0[17]))(1)= 0.35;  (*(q0[17]))(2)=  0.1;
    //     (*(q0[18]))(0)= -0.35;   (*(q0[18]))(1)= 0.35;  (*(q0[18]))(2)=  0.1;
    //     (*(q0[19]))(0)= -0.35;   (*(q0[19]))(1)=-0.35;  (*(q0[19]))(2)=  0.1;

    // 25*etage beads in cube
    unsigned int k;
    unsigned int cote  = 5;
    unsigned int etage  = 4;
    for (j = 0; j < etage; j++)
    {
      for (k = 0; k < cote; k++)
      {
        for (i = 0; i < cote; i++)
        {
          if (j % 2 == 0)
          {
            (*(q0[k * cote + i + j * cote * cote]))(0) = -0.6 + 3 * k * R;
            (*(q0[k * cote + i + j * cote * cote]))(1) = -0.6 + 3 * i * R;
            (*(q0[k * cote + i + j * cote * cote]))(2) = 0.2 + (2 * j) * R;
          }
          else
          {
            (*(q0[k * cote + i + j * cote * cote]))(0) = -0.6 + 3 * k * R;
            (*(q0[k * cote + i + j * cote * cote]))(1) = -0.6 + 3 * i * R + R / 4;
            (*(q0[k * cote + i + j * cote * cote]))(2) = 0.2 + (2 * j) * R;
          }
        }
      }
    }

    //     for (i=0;i<DSNUMBER;i++)
    //       (*(v0[i]))(2) = -1.;

    //     (*(q0[DSNUMBER-2]))(0)= 0.;   (*(q0[DSNUMBER-2]))(1)= 0.;   (*(q0[DSNUMBER-2]))(2)=  1.2*etage/4;    (*(v0[DSNUMBER-2]))(1) =  -1;(*(v0[DSNUMBER-2]))(2) =  -2;
    //     (*(q0[DSNUMBER-1]))(0)= -0.8;  (*(q0[DSNUMBER-1]))(1)= -0.8;  (*(q0[DSNUMBER-1]))(2)=  0.1;  (*(v0[DSNUMBER-1]))(0) =  10;(*(v0[DSNUMBER-1]))(1) =  10;



    for (i = 0; i < DSNUMBER; i++)
    {
      GLOB_tabLDS[i] = new LagrangianDS(i, *(q0[i]), *(v0[i]), *Mass);
      checkDS = allDS.insert(GLOB_tabLDS[i]);
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("3DDrawPlugin.so", "gravity");
    }

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------
    InteractionsSet allInteractions;

    vector<string> id;
    vector<string> id2;
    id.resize(Fact);
    id2.resize(DSNUMBER);

    DynamicalSystemsSet dsConcernedi;
    DynamicalSystemsSet dsConcerned2 ;
    CheckInsertInteraction checkInter;
    vector<Relation*> LLR(Fact);
    vector<Relation*> LLR1(DSNUMBER);
    vector<Relation*> LLR1_(DSNUMBER);
    vector<Relation*> LLR2(DSNUMBER);
    vector<Relation*> LLR2_(DSNUMBER);
    vector<Relation*> LLR3(DSNUMBER);
    vector<Relation*> LLR3_(DSNUMBER);

    // With friction

    NonSmoothLaw * nslaw1 = new NewtonImpactFrictionNSL(e, e, mu, 3);

    // Interaction beads and plan1 (OXY)

    if (obst_z_m)
    {
      SiconosVector *b1 = new SimpleVector(3);
      (*b1)(0) = Ground - R;
      SiconosMatrix *H1 = new SimpleMatrix(3, nDof);
      (*H1)(0, 2) = 1.0;
      (*H1)(1, 0) = 1.0;
      (*H1)(1, 4) = -R;
      (*H1)(2, 1) = 1.0;
      (*H1)(2, 3) =  R;

      for (i = 0; i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR1[i] = new LagrangianLinearR(*H1, *b1);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR1[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction beads and plan1 (-YOX)

    if (obst_z_p)
    {
      SiconosVector *b1_ = new SimpleVector(3);
      (*b1_)(0) = Top - R;
      SiconosMatrix *H1_ = new SimpleMatrix(3, nDof);
      (*H1_)(0, 2) = -1.0;
      (*H1_)(1, 0) = 1.0;
      (*H1_)(1, 4) = -R;
      (*H1_)(2, 1) = 1.0;
      (*H1_)(2, 3) =  R;

      for (i = 0; i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR1_[i] = new LagrangianLinearR(*H1_, *b1_);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR1_[i]));
        dsConcernedi.clear();
      }
    }


    // Interaction beads and plan2 (OXZ)

    if (obst_y_p)
    {
      SiconosVector *b2 = new SimpleVector(3);
      (*b2)(0) = Wall - R;
      SiconosMatrix *H2 = new SimpleMatrix(3, nDof);
      (*H2)(0, 1) = 1.0;
      (*H2)(1, 0) = 1.0;
      (*H2)(1, 5) = -R;
      (*H2)(2, 2) = 1.0;
      (*H2)(2, 3) =  R;

      for (i = 0; i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR2[i] = new LagrangianLinearR(*H2, *b2);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR2[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction beads and plan2 (-ZOX)

    if (obst_y_m)
    {
      SiconosVector *b2_ = new SimpleVector(3);
      (*b2_)(0) = Wall - R;
      SiconosMatrix *H2_ = new SimpleMatrix(3, nDof);
      (*H2_)(0, 1) = -1.0;
      (*H2_)(1, 0) = 1.0;
      (*H2_)(1, 5) = -R;
      (*H2_)(2, 2) = 1.0;
      (*H2_)(2, 3) =  R;

      for (i = 0; i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR2_[i] = new LagrangianLinearR(*H2_, *b2_);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR2_[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction beads and plan3 (OYZ)

    if (obst_x_p)
    {
      SiconosVector *b3 = new SimpleVector(3);
      (*b3)(0) = Wall - R;
      SiconosMatrix *H3 = new SimpleMatrix(3, nDof);
      (*H3)(0, 0) = 1.0;
      (*H3)(1, 1) = 1.0;
      (*H3)(1, 5) = -R;
      (*H3)(2, 2) = 1.0;
      (*H3)(2, 4) =  R;

      for (i = 0; i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR3[i] = new LagrangianLinearR(*H3, *b3);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR3[i]));
        dsConcernedi.clear();
      }
    }
    // Interaction beads and plan3 (-ZOY)

    if (obst_x_m)
    {
      SiconosVector *b3_ = new SimpleVector(3);
      (*b3_)(0) = Wall - R;
      SiconosMatrix *H3_ = new SimpleMatrix(3, nDof);
      (*H3_)(0, 0) = -1.0;
      (*H3_)(1, 1) = 1.0;
      (*H3_)(1, 5) = -R;
      (*H3_)(2, 2) = 1.0;
      (*H3_)(2, 4) =  R;

      for (i = 0; i < DSNUMBER; i++)
      {
        dsConcernedi.insert(GLOB_tabLDS[i]);
        ostringstream ostr;
        ostr << i;
        id2[i] = ostr.str();
        LLR3_[i] = new LagrangianLinearR(*H3_, *b3_);
        checkInter = allInteractions.insert(new Interaction(id2[i], dsConcernedi, i, 3, nslaw1, LLR3_[i]));
        dsConcernedi.clear();
      }
    }

    // Interaction between beads

    // frictional contact condition between beads
    NonSmoothLaw * nslaw2 = new NewtonImpactFrictionNSL(e2, e2, mu, 3);

    l = 0;
    for (i = 0; i < DSNUMBER; i++)
    {
      dsConcerned2.insert(GLOB_tabLDS[i]);
      for (j = 0; j < DSNUMBER; j++)
      {
        if (j > i)
        {
          dsConcerned2.insert(GLOB_tabLDS[j]);
          ostringstream ostr;
          ostr << l;
          id[l] = ostr.str();
          LLR[l] = new LagrangianScleronomousR("3DDrawPlugin:h0", "3DDrawPlugin:G0");
          checkInter = allInteractions.insert(new Interaction(id[l], dsConcerned2, l, 3, nslaw2, LLR[l]));
          dsConcerned2.erase(GLOB_tabLDS[j]);
          l = l + 1;
        }
      }
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

    TimeDiscretisation * GLOB_T = new TimeDiscretisation(h, multiBeads);
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
    int N = GLOB_T->getNSteps(); // Number of time steps

    unsigned int outputSize = 1 + 2 * DSNUMBER;
    SimpleMatrix dataPlot(N + 1, outputSize);
    int iter_k = 0; // index for output.
    dataPlot(iter_k, 0) = multiBeads->getT0();
    for (i = 0; i < DSNUMBER; ++i)
    {
      dataPlot(iter_k, 1 + 2 * i) = GLOB_tabLDS[i]->getQ()(2);
      dataPlot(iter_k, 2 + 2 * i) = GLOB_tabLDS[i]->getVelocity()(2);
      //  dataPlot(iter_k,3+3*i) = (multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);
    }

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    boost::progress_display show_progress(N);
    for (iter_k = 1 ; iter_k < N + 1  ; ++iter_k)
    {
      GLOB_SIM->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(iter_k, 0) =  GLOB_SIM->getNextTime();
      for (i = 0; i < DSNUMBER; ++i)
      {
        dataPlot(iter_k, 1 + 2 * i) = GLOB_tabLDS[i]->getQ()(2);
        dataPlot(iter_k, 2 + 2 * i) = GLOB_tabLDS[i]->getVelocity()(2);
        // dataPlot(iter_k,3+3*i)=(multiBeads->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(1)->getLambda(1))(0);
      }
      GLOB_SIM->nextStep();
      ++show_progress;
    }
    cout << "End of computation - Number of iterations done: " << iter_k - 1 << endl;



    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    //    io.write(dataPlot,"noDim");

    // --- Free memory ---
    delete nsds;
    delete multiBeads;
    delete osnspb;
    delete OSI;
    delete GLOB_SIM;
    delete GLOB_T;
    for (i = 0; i < Fact; ++i)
      delete LLR[i];
    for (i = 0; i < DSNUMBER; ++i)
    {
      delete GLOB_tabLDS[i];
      delete q0[i];
      delete v0[i];
      delete LLR1[i];
      delete LLR1_[i];
      delete LLR2[i];
      delete LLR2_[i];
      delete LLR3[i];
      delete LLR3_[i];
    }
    delete Mass;
    delete nslaw1;
    delete nslaw2;

  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in beads.cpp" << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl << endl;
}






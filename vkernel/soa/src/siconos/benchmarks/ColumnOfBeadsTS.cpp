/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
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

/*!\file ColumnOfBeadsTS.cpp
  \brief \ref EMColumnOfBeads - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include "OneStepNSProblem.hpp"
#include "SolverOptions.h"
#include "SiconosKernel.hpp"
#include <chrono>
#include <stdlib.h> // atoi

using namespace std;

int main(int argc, char* argv[])
{
  unsigned int nBeads = atoi(argv[1]);
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 10.0;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 1.0;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for MoreauJeanOSI integrator
    double R = 0.1; // Ball radius
    double m = 1; // Ball mass
    double g = 9.81; // Gravity
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    // Number of Beads
    double initialGap = 1.0;
    double alert = 0.02;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    std::vector<SP::SiconosVector> q0(nBeads);
    std::vector<SP::SiconosVector> v0(nBeads);

    for(unsigned int i = 0; i < nBeads; i++)
    {
      (q0[i]).reset(new SiconosVector(nDof));
      (v0[i]).reset(new SiconosVector(nDof));
      (q0[i])->setValue(0, position_init + i * initialGap);
      (v0[i])->setValue(0, velocity_init);
    }

    // -- The dynamical system --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(0) = -m * g;


    std::vector<SP::LagrangianLinearTIDS> beads(nBeads);
    for(unsigned int i = 0; i < nBeads; i++)
    {
      beads[i].reset(new LagrangianLinearTIDS(q0[i], v0[i], Mass));
      // -- Set external forces (weight) --
      beads[i]->setFExtPtr(weight);
    }


    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.9;

    // Interaction ball-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(2, nDof));
    (*H)(0, 0) = 1.0;
    (*H)(1, 1) = 1.0;
    (*H)(1, 2) = - R;
    SP::SiconosVector b(new SiconosVector(2));
    (*b)(0) = -R;

    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(e, 0, 0, 2));
    SP::Relation relation(new LagrangianLinearTIR(H, b));

    SP::Interaction inter;


    // beads/beads interactions
    SP::SimpleMatrix HOfBeads(new SimpleMatrix(2, 2 * nDof));
    (*HOfBeads)(0, 0) = -1.0;
    (*HOfBeads)(0, 3) = 1.0;
    (*HOfBeads)(1, 1) = 1.0;
    (*HOfBeads)(1, 3) = 1.0;
    (*HOfBeads)(1, 2) = - R;
    (*HOfBeads)(1, 5) = - R;
    SP::SiconosVector bOfBeads(new SiconosVector(2));
    (*bOfBeads)(0) = -2 * R;

    std::vector<SP::Relation > relationOfBeads(nBeads - 1);
    std::vector<SP::Interaction > interOfBeads(nBeads - 1);


    // --------------------------------------
    // ---      Model and simulation      ---
    // --------------------------------------

    SP::NonSmoothDynamicalSystem columnOfBeads(new NonSmoothDynamicalSystem(t0, T));
    // add the dynamical system in the non smooth dynamical system
    for(unsigned int i = 0; i < nBeads; i++)
    {
      columnOfBeads->insertDynamicalSystem(beads[i]);
    }

    // --  (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    //SP::OneStepNSProblem osnspb(new LCP());
    SP::OneStepNSProblem osnspb(new FrictionContact(2, SICONOS_FRICTION_2D_NSGS));
    std::static_pointer_cast<LinearOSNS>(osnspb)->setMStorageType(NM_SPARSE);
    std::static_pointer_cast<LinearOSNS>(osnspb)->setAssemblyType(REDUCED_DIRECT);
    std::static_pointer_cast<LinearOSNS>(osnspb)->setKeepLambdaAndYState(false);
    std::static_pointer_cast<OneStepNSProblem>(osnspb)->numericsSolverOptions()->iparam[SICONOS_IPARAM_MAX_ITER] = 1;
    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(columnOfBeads, t, OSI, osnspb));


    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 1 + nBeads * 4;
    SimpleMatrix dataPlot(N + 1, outputSize);

    dataPlot(0, 0) = columnOfBeads->t0();

    for(unsigned int i = 0; i < nBeads; i++)
    {
      dataPlot(0, 1 + i * 2) = (beads[i]->q())->getValue(0);
      dataPlot(0, 2 + i * 2) = (beads[i]->velocity())->getValue(0);
      //      dataPlot(0,3+i*4) = (beads[i]->p(1))->getValue(0);
    }

    // for (unsigned int i =1; i< nBeads; i++)
    // {
    // dataPlot(0,4+i*4) = (interOfBeads[i-1]->lambda(1))->getValue(0);
    // }

    // --- Time loop ---
    // cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;

    inter.reset(new Interaction(nslaw, relation));
    columnOfBeads->link(inter, beads[0]);
    for(unsigned int i=0; i<nBeads-1; ++i)
    {
      relationOfBeads[i].reset(new LagrangianLinearTIR(HOfBeads, bOfBeads));
      interOfBeads[i].reset(new Interaction(nslaw, relationOfBeads[i]));
      columnOfBeads->link(interOfBeads[i], beads[i], beads[i+1]);
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    int ncontact = 0 ;
    while(s->hasNextEvent())
    {
      s->computeOneStep();

      // --- Get values to be plotted ---
/*      dataPlot(k, 0) =  s->nextTime();
      for(unsigned int i = 0; i < nBeads; i++)
      {
        dataPlot(k, 1 + i * 2) = (beads[i]->q())->getValue(0);
        dataPlot(k, 2 + i * 2) = (beads[i]->velocity())->getValue(0);
      }
      // for (unsigned int i =1; i< nBeads; i++)
      // {
      //   dataPlot(k,4+i*4) = (interOfBeads[i-1]->lambda(1))->getValue(0);
      // }
      // for (unsigned int i =1; i< nBeads; i++)
      // {
      //   std::cout <<  (interOfBeads[i-1]->y(0))->getValue(0) << std::endl ;
      // }
      */
      s->nextStep();
      cout << k << std::endl;
      k++;
    }
//    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
//    cout << "Computation Time " << endl;;
    end = std::chrono::system_clock::now();
    int elapsed = std::chrono::duration_cast<std::chrono::milliseconds>
                  (end-start).count();
    cout << "RESULTS : -- " << elapsed << " --" << endl;
    // --- Output files ---
/*    cout << "====> Output file writing ..." << endl;
      dataPlot.resize(k, outputSize);
      ioMatrix::write("ColumnOfbeadsTS.dat", "ascii", dataPlot, "noDim");*/

//    double error=0.0, eps=1e-12;
//    if((error=ioMatrix::compareRefFile(dataPlot, "ColumnOfBeadsTS.ref", eps)) >= 0.0
//        && error > eps)
//      return 1;

  }

  catch(...)
  {
    Siconos::exception::process();
    return 1;
  }



}

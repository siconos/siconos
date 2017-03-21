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

/*!\file ColumnOfBeadsTS.cpp
  \brief \ref EMColumnOfBeads - C++ input file, Time-Stepping version -
  V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;           // degrees of freedom for the ball
    double t0 = 0;                   // initial computation time
    double T = 2.0;                  // final computation time
    double h = 0.0005;                // time step
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
    unsigned int nBeads = 10;
    double initialGap = 0.25;
    double alert = 0.02;

    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m;
    (*Mass)(1, 1) = m;
    (*Mass)(2, 2) = 3. / 5 * m * R * R;

    // -- Initial positions and velocities --
    std::vector<SP::SiconosVector> q0(nBeads);
    std::vector<SP::SiconosVector> v0(nBeads);

    for (unsigned int i = 0; i < nBeads; i++)
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
    for (unsigned int i = 0; i < nBeads; i++)
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
    SP::SimpleMatrix H(new SimpleMatrix(1, nDof));
    (*H)(0, 0) = 1.0;
    SP::SiconosVector b(new SiconosVector(1));
    (*b)(0) = -R;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H, b));

    SP::Interaction inter(new Interaction(nslaw, relation));


    // beads/beads interactions
    SP::SimpleMatrix HOfBeads(new SimpleMatrix(1, 2 * nDof));
    (*HOfBeads)(0, 0) = -1.0;
    (*HOfBeads)(0, 3) = 1.0;
    SP::SiconosVector bOfBeads(new SiconosVector(1));
    (*bOfBeads)(0) = -2 * R;

    std::vector<SP::Relation > relationOfBeads(nBeads - 1);
    std::vector<SP::Interaction > interOfBeads(nBeads - 1);
    for (unsigned int i =0; i< nBeads-1; i++)
      {
      relationOfBeads[i].reset(new LagrangianLinearTIR(HOfBeads,bOfBeads));
      interOfBeads[i].reset(new Interaction(nslaw, relationOfBeads[i]));
      }

    // --------------------------------------
    // ---      Model and simulation      ---
    // --------------------------------------
    SP::Model columnOfBeads(new Model(t0, T));
    
    // add the dynamical system in the non smooth dynamical system
    for (unsigned int i = 0; i < nBeads; i++)
    {
      columnOfBeads->nonSmoothDynamicalSystem()->insertDynamicalSystem(beads[i]);
    }

    columnOfBeads->nonSmoothDynamicalSystem()->link(inter, beads[0]);
    // link the interaction and the dynamical systems
    for (unsigned int i =0; i< nBeads-1; i++)
      columnOfBeads->nonSmoothDynamicalSystem()->link(interOfBeads[i],beads[i], beads[i+1]);

    // --  (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::LinearOSNS osnspb(new LCP(SICONOS_LCP_NSGS_SBM));
    osnspb->setMStorageType(1);
    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI, osnspb));

    columnOfBeads->setSimulation(s);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    columnOfBeads->initialize();

    int N = ceil((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 1 + nBeads * 4;
    SimpleMatrix dataPlot(N + 1, outputSize);

    dataPlot(0, 0) = columnOfBeads->t0();

    for (unsigned int i = 0; i < nBeads; i++)
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
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    int ncontact = 0 ;
    bool isOSNSinitialized = false;
    while (s->hasNextEvent())
    {
      // Rough contact detection
      for (unsigned int i = 0; i < nBeads - 1; i++)
      {
        // Between first bead and plane
        if (abs(((beads[i])->q())->getValue(0) - R) < alert)
        {
          if (!inter)
          {
            ncontact++;
            // std::cout << "Number of contact = " << ncontact << std::endl;

            inter.reset(new Interaction(nslaw, relation));
            columnOfBeads->nonSmoothDynamicalSystem()->link(inter, beads[0]);
            s->initializeInteraction(s->nextTime(), inter);

            if (!isOSNSinitialized)
            {
              s->initOSNS();
              isOSNSinitialized = true;
            }

            assert(inter->y(0)->getValue(0) >= 0);
          }
        }

        // Between two beads
        if (abs(((beads[i + 1])->q())->getValue(0) - ((beads[i])->q())->getValue(0) - 2 * R) < alert)
        {
          //std::cout << "Alert distance for declaring contact = ";
          //std::cout << abs(((beads[i])->q())->getValue(0)-((beads[i+1])->q())->getValue(0))   <<std::endl;
          if (!interOfBeads[i].get())
          {
            ncontact++;
            // std::cout << "Number of contact = " << ncontact << std::endl;

            relationOfBeads[i].reset(new LagrangianLinearTIR(HOfBeads, bOfBeads));
            interOfBeads[i].reset(new Interaction(nslaw, relationOfBeads[i]));

            columnOfBeads->nonSmoothDynamicalSystem()->link(interOfBeads[i], beads[i], beads[i+1]);
            s->initializeInteraction(s->nextTime(), interOfBeads[i]);

            if (!isOSNSinitialized)
            {
              s->initOSNS();
              isOSNSinitialized = true;
            }

            assert(interOfBeads[i]->y(0)->getValue(0) >= 0);
          }
        }
      }

      s->computeOneStep();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      for (unsigned int i = 0; i < nBeads; i++)
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

      s->nextStep();
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();

    ioMatrix::read("result-SBM.ref", "ascii", dataPlotRef);

    cout << "====> Comparison with reference file ..." << endl;
    std::cout << "Error w.r.t. reference file : " << (dataPlot - dataPlotRef).normInf() << std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
      return 1;
    }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in ColumnOfBeadsTS.cpp" << endl;
  }



}

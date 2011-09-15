/* Siconos-sample , Copyright INRIA 2005-2011.
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

/*!\file ImpactingBeam.cpp
  V. Acary

  A Beam bouncing on the ground
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
    unsigned int nDof = 100;           // degrees of freedom for the beam
    double t0 = 0;                   // initial computation time
    double T = 0.05;                  // final computation time
    double h = 0.00005;                // time step
    double position_init = 0.001;      // initial position
    double velocity_init =  -.1;      // initial velocity
    double theta = 0.500000000001;              // theta for Moreau integrator
    //theta = .;
    double E = 210e6; // young Modulus
    double S = 0.000314; //  Beam Section 1 cm  for the diameter
    //S=0.1;
    double L = 1.0; // length of the  beam
    double rho = 7800 ; // specific mass
    //rho=1.0;
    double g = 9.81; // Gravity
    g = 0.0;
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    double l = L / nDof; // length of an element

    //    cout << "length" << l <<  endl;

    SP::SiconosMatrix SparseMass(new SimpleMatrix(nDof, nDof, Siconos::SPARSE, nDof));
    SP::SiconosMatrix SparseStiffness(new SimpleMatrix(nDof, nDof, Siconos::SPARSE, 3 * nDof));




    SparseMass->setValue(0, 0, 1.0 / 3.0);
    SparseMass->setValue(0, 1, 1.0 / 6.0);
    SparseStiffness->setValue(0, 0, 1.0);
    SparseStiffness->setValue(0, 1, -1.0);

    for (int i = 1; i < nDof - 1; i++)
    {
      SparseMass->setValue(i, i, 2.0 / 3.0);
      SparseMass->setValue(i, i - 1, 1.0 / 6.0);
      SparseMass->setValue(i, i + 1, 1.0 / 6.0);


      SparseStiffness->setValue(i, i, 2.0);
      SparseStiffness->setValue(i, i - 1, -1.0);
      SparseStiffness->setValue(i, i + 1, -1.0);
    }

    SparseMass->setValue(nDof - 1, nDof - 1, 1.0 / 3.0);
    SparseMass->setValue(nDof - 1, nDof - 2, 1.0 / 6.0);


    SparseStiffness->setValue(nDof - 1, nDof - 2, -1.0);
    SparseStiffness->setValue(nDof - 1, nDof - 1, 1.0);

    *SparseMass  *= rho * S * l;
    *SparseStiffness  *= E * S / l;

    //      SparseMass->display();
    //      SparseStiffness->display();


    // -- Initial positions and velocities --
    SP::SimpleVector q0(new SimpleVector(nDof, position_init));
    SP::SimpleVector v0(new SimpleVector(nDof, velocity_init));

    // -- The dynamical system --
    SP::LagrangianLinearTIDS beam(new LagrangianLinearTIDS(q0, v0, SparseMass));

    // -- Set stiffness matrix (weight) --
    beam->setKPtr(SparseStiffness);



    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof, -g * rho * S / l));


    beam->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction beam-floor
    //
    SP::SiconosMatrix H(new SimpleMatrix(1, nDof));
    (*H)(0, 0) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(1, nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::Model impactingBeam(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    impactingBeam->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam);

    // link the interaction and the dynamical system
    impactingBeam->nonSmoothDynamicalSystem()->link(inter, beam);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::Moreau OSI(new Moreau(beam, theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI, osnspb));

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    impactingBeam->initialize(s);
    int N = floor((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 7;
    SimpleMatrix dataPlot(N, outputSize);

    SP::SiconosVector q = beam->q();
    SP::SiconosVector v = beam->velocity();
    SP::SiconosVector p = beam->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = impactingBeam->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);

    SP::SiconosVector tmp(new SimpleVector(nDof));

    prod(*SparseStiffness, *q, *tmp, true);
    double potentialEnergy = inner_prod(*q,   *tmp);
    prod(*SparseMass, *v, *tmp, true);
    double kineticEnergy = inner_prod(*v, *tmp);

    dataPlot(0, 5) = potentialEnergy;
    dataPlot(0, 6) = kineticEnergy;

    //    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
    //     std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    //    while (s->nextTime() < T)
    while (k < N)
    {
      s->computeOneStep();

      //       std::cout << "position"  << std::endl;
      //       q->display();
      //       std::cout << "velocity"  << std::endl;
      //       v->display();

      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      dataPlot(k, 4) = (*lambda)(0);


      prod(*SparseStiffness, *q, *tmp, true);
      potentialEnergy = inner_prod(*q,   *tmp);
      prod(*SparseMass, *v, *tmp, true);
      kineticEnergy = inner_prod(*v, *tmp);

      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

      //      std::cout << "q" << std::endl;
      //       q->display();


      //       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
      //       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



      s->nextStep();
      ++show_progress;
      k++;
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");



  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in ImpactingBeamTS.cpp" << endl;
  }



}

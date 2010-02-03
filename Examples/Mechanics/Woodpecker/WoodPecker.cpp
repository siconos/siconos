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

#include "SiconosKernel.hpp"
#include "WoodPeckerConsts.h"

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer t;
  t.restart();
  try
  {
    // ================= Model definition =================

    // User-defined main parameters
    unsigned int nDof = 3;            // degrees of freedom
    double t0 = 0;                    // initial computation time
    double T = 0.3;                   // final computation time
    double h = 0.00002;               // time step
    double theta = 0.5;               // theta for Moreau integrator;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    SP::SiconosMatrix Mass, K, C;
    Mass.reset(new SimpleMatrix(nDof, nDof));
    (*Mass)(0, 0) = m_S + m_M;
    (*Mass)(0, 1) = m_S * l_M;
    (*Mass)(0, 2) = m_S * l_G;
    (*Mass)(1, 0) = m_S * l_M;
    (*Mass)(1, 1) = J_M + m_S * l_M * l_M;
    (*Mass)(1, 2) = m_S * l_M * l_G;
    (*Mass)(2, 0) = m_S * l_G;
    (*Mass)(2, 1) = m_S * l_M * l_G;
    (*Mass)(2, 2) = J_S + m_S * l_G * l_G;
    K.reset(new SimpleMatrix(nDof, nDof));
    (*K)(1, 1) = c_phi;
    (*K)(1, 2) = -c_phi;
    (*K)(2, 1) = -c_phi;
    (*K)(2, 2) = c_phi;
    C.reset(new SimpleMatrix(nDof, nDof));

    // -- Initial positions and velocities --
    SP::SimpleVector q0(new SimpleVector(nDof));
    (*q0)(0) = y_0;
    (*q0)(1) = phi_M_0;
    (*q0)(2) = phi_S_0;

    SP::SimpleVector velocity0(new SimpleVector(nDof));
    (*velocity0)(0) = v_0;
    (*velocity0)(1) = omega_M_0;
    (*velocity0)(2) = omega_S_0;

    SP::LagrangianDS dynamicalSystem(new LagrangianLinearTIDS(q0, velocity0, Mass, K, C));
    dynamicalSystem->setComputeFExtFunction("WoodPeckerPlugin.so", "FExt");

    DynamicalSystemsSet allDS;
    allDS.insert(dynamicalSystem);

    // --------------------
    // --- Interactions ---
    // --------------------

    SP::SiconosMatrix H1(new SimpleMatrix(2, nDof));
    (*H1)(0, 0) = 0;
    (*H1)(0, 1) = 0;
    (*H1)(0, 2) = -h_S;
    (*H1)(1, 0) = 1;
    (*H1)(1, 1) = l_M;
    (*H1)(1, 2) = l_G - l_S;
    SP::SimpleVector b1(new SimpleVector(2));
    (*b1)(0) = l_M + l_G - l_S - r_O;
    (*b1)(1) = 0;
    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eps_N_1, eps_T_123, mu_123, 2));
    SP::Relation relation1(new LagrangianLinearTIR(H1, b1));

    SP::SiconosMatrix H23(new SimpleMatrix(4, nDof));
    (*H23)(0, 0) = 0;
    (*H23)(0, 1) = h_M;
    (*H23)(0, 2) = 0;
    (*H23)(1, 0) = 1;
    (*H23)(1, 1) = r_M;
    (*H23)(1, 2) = 0;
    (*H23)(2, 0) = 0;
    (*H23)(2, 1) = -h_M;
    (*H23)(2, 2) = 0;
    (*H23)(3, 0) = 1;
    (*H23)(3, 1) = r_M;
    (*H23)(3, 2) = 0;
    SP::SimpleVector b23(new SimpleVector(4));
    (*b23)(0) = r_M - r_O;
    (*b23)(1) = 0;
    (*b23)(2) = r_M - r_O;
    (*b23)(3) = 0;
    SP::NonSmoothLaw nslaw23(new NewtonImpactFrictionNSL(eps_N_23, eps_T_123, mu_123, 2));
    SP::Relation relation23(new LagrangianLinearTIR(H23, b23));

    SP::Interaction I1(new Interaction("contact1", allDS, 0, 2, nslaw1, relation1));

    SP::Interaction I2(new Interaction("contact23", allDS, 1, 4, nslaw23, relation23));

    InteractionsSet allInteractions;
    allInteractions.insert(I1);
    allInteractions.insert(I2);
    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------

    bool isBVP = 0;
    SP::NonSmoothDynamicalSystem nsds(new NonSmoothDynamicalSystem(allDS, allInteractions, isBVP));

    // -------------
    // --- Model ---
    // -------------

    SP::Model model(new Model(t0, T));
    model->setNonSmoothDynamicalSystemPtr(nsds);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator vOSI(new Moreau(dynamicalSystem, theta));
    s->insertIntegrator(vOSI);


    SP::OneStepNSProblem osnspb(new FrictionContact(2));
    s->insertNonSmoothProblem(osnspb);


    cout << "=== End of model loading === " << endl;

    // ================= Computation =================

    // --- Simulation initialization ---
    model->initialize(s);

    cout << "End of model initialisation" << endl;


    int k = 0;
    int N = (int)((T - t0) / h);

    // --- Get the values to be plotted ---
    unsigned int outputSize = 7;
    SimpleMatrix dataPlot(N + 1, outputSize);
    dataPlot(k, 0) = t0;
    for (int i = 0; i < (int)nDof; i++)
    {
      dataPlot(k, 2 * i + 1) = (*dynamicalSystem->q())(i);
      dataPlot(k, 2 * i + 2) = (*dynamicalSystem->velocity())(i);
    }

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    while (k < N)
    {

      // get current time step
      k++;

      // solve ...
      s->computeOneStep();

      // get values
      dataPlot(k, 0) = s->nextTime();
      for (int i = 0; i < (int)nDof; i++)
      {
        dataPlot(k, 2 * i + 1) = (*dynamicalSystem->q())(i);
        dataPlot(k, 2 * i + 2) = (*dynamicalSystem->velocity())(i);
      }

      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
    }
    cout << "End of computation - Number of iterations done: " << k << endl;

    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'WoodPecker\'" << endl;
  }
  cout << "Computation Time " << t.elapsed()  << endl;
}

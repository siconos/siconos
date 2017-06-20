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
    double theta = 0.5;               // theta for MoreauJeanOSI integrator;

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
    SP::SiconosVector q0(new SiconosVector(nDof));
    (*q0)(0) = y_0;
    (*q0)(1) = phi_M_0;
    (*q0)(2) = phi_S_0;

    SP::SiconosVector velocity0(new SiconosVector(nDof));
    (*velocity0)(0) = v_0;
    (*velocity0)(1) = omega_M_0;
    (*velocity0)(2) = omega_S_0;

    SP::LagrangianDS dynamicalSystem(new LagrangianLinearTIDS(q0, velocity0, Mass, K, C));
    dynamicalSystem->setComputeFExtFunction("WoodPeckerPlugin", "FExt");

    // --------------------
    // --- Interactions ---
    // --------------------

    SP::SimpleMatrix H1(new SimpleMatrix(2, nDof));
    (*H1)(0, 0) = 0;
    (*H1)(0, 1) = 0;
    (*H1)(0, 2) = -h_S;
    (*H1)(1, 0) = 1;
    (*H1)(1, 1) = l_M;
    (*H1)(1, 2) = l_G - l_S;
    SP::SiconosVector b1(new SiconosVector(2));
    (*b1)(0) = l_M + l_G - l_S - r_O;
    (*b1)(1) = 0;
    SP::NonSmoothLaw nslaw1(new NewtonImpactFrictionNSL(eps_N_1, eps_T_123, mu_123, 2));

    SP::Relation relation1(new LagrangianLinearTIR(H1, b1));

    SP::SimpleMatrix H2(new SimpleMatrix(2, nDof));
    (*H2)(0, 0) = 0;
    (*H2)(0, 1) = h_M;
    (*H2)(0, 2) = 0;
    (*H2)(1, 0) = 1;
    (*H2)(1, 1) = r_M;
    (*H2)(1, 2) = 0;
    SP::SimpleMatrix H3(new SimpleMatrix(2, nDof));
    (*H3)(0, 0) = 0;
    (*H3)(0, 1) = -h_M;
    (*H3)(0, 2) = 0;
    (*H3)(1, 0) = 1;
    (*H3)(1, 1) = r_M;
    (*H3)(1, 2) = 0;
    SP::SiconosVector b2(new SiconosVector(2));
    (*b2)(0) = r_M - r_O;
    (*b2)(1) = 0;
    SP::SiconosVector b3(new SiconosVector(2));
    (*b3)(0) = r_M - r_O;
    (*b3)(1) = 0;

    SP::NonSmoothLaw nslaw23(new NewtonImpactFrictionNSL(eps_N_23, eps_T_123, mu_123, 2));

    SP::Relation relation2(new LagrangianLinearTIR(H2, b2));
    SP::Relation relation3(new LagrangianLinearTIR(H3, b3));

    SP::Interaction I1(new Interaction(nslaw1, relation1));

    SP::Interaction I2(new Interaction(nslaw23, relation2));

    SP::Interaction I3(new Interaction(nslaw23, relation3));

    // -------------
    // --- Model ---
    // -------------

    SP::Model model(new Model(t0, T));
    model->nonSmoothDynamicalSystem()->insertDynamicalSystem(dynamicalSystem);
    model->nonSmoothDynamicalSystem()->link(I1, dynamicalSystem);
    model->nonSmoothDynamicalSystem()->link(I2, dynamicalSystem);
    model->nonSmoothDynamicalSystem()->link(I3, dynamicalSystem);

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator vOSI(new MoreauJeanOSI(theta));
    s->insertIntegrator(vOSI);


    SP::OneStepNSProblem osnspb(new FrictionContact(2));
    s->insertNonSmoothProblem(osnspb);

    model->setSimulation(s);

    cout << "=== End of model loading === " << endl;

    // ================= Computation =================

    // --- Simulation initialization ---
    model->initialize();

    cout << "End of model initialisation" << endl;


    int k = 0;
    int N = floor((T - t0) / h);

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
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("Woodpecker.ref", "ascii", dataPlotRef);
    double error = (dataPlot - dataPlotRef).normInf()/  dataPlotRef.normInf();
    std::cout << "Error = "<< error <<std::endl;
    if (error > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
    return 1;
  }
  catch (...)
  {
    cout << "Exception caught in \'WoodPecker\'" << endl;
    return 1;
  }
  cout << "Computation Time " << t.elapsed()  << endl;
  return 0;
}

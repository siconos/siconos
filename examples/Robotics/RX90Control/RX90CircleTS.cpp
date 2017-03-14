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


// =============================== Robot arm sample  ===============================
//
// Keywords: LagrangianDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP.
//
// =============================================================================================

#include "SiconosKernel.hpp"
#include <math.h>

#define PI 3.14159265

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 6;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 6.0;                   // final computation time
    double h = 5e-3;                // time step
    double criterion = 1e-8;
    unsigned int maxIter = 200000;
    double e = 0.0;                  // nslaw

    // -> mind to set the initial conditions below.

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    // unsigned int i;

    // --- DS: manipulator arm ---

    // The dof are angles between ground and arm and between differents parts of the arm. (See corresponding .pdf for more details)

    // Initial position (angles in radian)
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();
    (*q0)(1) = PI / 3;
    (*q0)(2) = -PI / 6;
    (*q0)(4) = PI / 6;
    (*v0)(0) = -0.34;
    (*v0)(3) = 0.59;
    (*v0)(5) = -0.34;


    SP::LagrangianDS arm(new LagrangianDS(q0, v0, "RX90Plugin:mass"));

    // external plug-in
    //    arm->setComputeMassFunction("RX90Plugin","mass");
    arm->setComputeFGyrFunction("RX90Plugin", "FGyr");
    arm->setComputeJacobianFGyrqFunction("RX90Plugin", "jacobianVFGyr");
    arm->setComputeJacobianFGyrqDotFunction("RX90Plugin", "jacobianFGyrq");
    arm->setComputeFIntFunction("RX90Plugin", "U");
    arm->setComputeJacobianFIntqFunction("RX90Plugin", "jacobFintV");
    arm->setComputeJacobianFIntqDotFunction("RX90Plugin", "jacobFintQ");

    // creating Z parameter computed in Actuators and used in FInt
    SP::SiconosVector torques(new SiconosVector(nDof));
    torques->zero();
    arm->setzPtr(torques);

    // -------------------
    // --- Interactions---
    // -------------------

    //  - one with linear relation to define joints limits
    //  Both with newton impact nslaw.

    // -- relations --

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));

    std::vector<SP::SiconosMatrix> Hvector(12);
    std::vector<SP::SiconosVector> bvector(12);


    SP::SiconosVector b(new SiconosVector(12));
    for (unsigned int i = 0; i < nDof; i++)
    {
      Hvector[2 * i].reset(new SimpleMatrix(1, 6));
      Hvector[2 * i + 1].reset(new SimpleMatrix(1, 6));
      Hvector[2 * i]->zero();
      Hvector[2 * i + 1]->zero();
      (*(Hvector[2 * i]))(0, i) = 1;
      (*(Hvector[2 * i + 1]))(0, i) = -1;

    }


    (*b)(0) = PI * 167.0 / 180.0;
    (*b)(1) = (*b)(0);
    (*b)(2) = PI * 137.5 / 180.0;
    (*b)(3) = (*b)(2);
    (*b)(4) = PI * 142.5 / 180.0;
    (*b)(5) = (*b)(4);
    (*b)(6) = PI * 270.0 / 180.0;
    (*b)(7) = (*b)(6);
    (*b)(8) = PI * 112.5 / 180.0;
    (*b)(9) = (*b)(8);
    (*b)(10) = PI * 270.0 / 180.0;
    (*b)(11) = (*b)(10);
    std::vector<SP::Relation> relationVector(12);
    std::vector<SP::Interaction> interactionVector(12);
    // -------------
    // --- Model ---
    // -------------

    SP::Model RX90(new Model(t0, T));
    RX90->nonSmoothDynamicalSystem()->insertDynamicalSystem(arm);

    for (unsigned int i = 0; i < 2 * nDof; i++)
    {
      bvector[i].reset(new SiconosVector(1));
      // (*(bvector[i])) (1) = *b(i);
      (bvector[i])->setValue(0, (*b)(i));
      relationVector[i].reset(new LagrangianLinearTIR(Hvector[i], bvector[i]));
      interactionVector[i].reset(new Interaction(nslaw, relationVector[i]));
      RX90->nonSmoothDynamicalSystem()->link(interactionVector[i], arm);
    }
   

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- Actuation (ici?) --
    //le premier evenement doit etre un evenement doit etre un evenement
    //de calcul, et pas de Actuators, ni sensors, c'est pourquoi le pas
    //de temps de la simu (h = 5e-3) est plus petit que le retard (6e-3)
    SP::TimeDiscretisation Sampling(new TimeDiscretisation(t0, 2 * h));

    unsigned int N = (unsigned int)((T - t0) / h + 1);

    std::vector<double> *tmp = new vector<double>(N);

    (*tmp)[0] = t0;
    for (unsigned int i = 1; i < tmp->size() - 1; i++)
      (*tmp)[i] = t0 + (i - 1) * 2 * h + 6e-3;
    (*tmp)[tmp->size() - 1] = T;

    SP::TimeDiscretisation delay(new TimeDiscretisation(*tmp));

    //Creation du control et ajout du Sensor et Actuator
    SP::ControlManager control(new ControlManager(RX90));
    control->addSensor(2, Sampling);
    control->addActuator(2, delay);
    (*(control->getActuators().begin()))->addSensorPtr(*((control->getSensors()).begin()));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    SP::OneStepIntegrator OSI(new MoreauJeanOSI(arm, 0.5));
    s->insertIntegrator(OSI);

    // -- OneStepNsProblem --
    SP::OneStepNSProblem osnsp(new LCP(SICONOS_LCP_LEMKE));
    s->insertNonSmoothProblem(osnsp);

    cout << "=== End of model loading === " << endl;

    // =========================== End of model definition ===========================


    // ================================= Computation
    // --- Simulation initialization ---


    RX90->initialize();
    cout << "End of simulation initialisation" << endl;

    int k = 0;
    //    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 13;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    SP::SiconosVector q = arm->q();
    SP::SiconosVector v = arm->velocity();
    SP::EventsManager eventsManager = s->eventsManager();

    dataPlot(k, 0) =  RX90->t0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (*q)(2);
    dataPlot(k, 4) = (*q)(3);
    dataPlot(k, 5) = (*q)(4);
    dataPlot(k, 6) = (*q)(5);
    dataPlot(k, 7) = (*v)(0);
    dataPlot(k, 8) = (*v)(1);
    dataPlot(k, 9) = (*v)(2);
    dataPlot(k, 10) = (*v)(3);
    dataPlot(k, 11) = (*v)(4);
    dataPlot(k, 12) = (*v)(5);


    while (s->hasNextEvent())
    {
      s->advanceToEvent();
      s->processEvents();
      // get current time step
      std::cout << k << std::endl;
      if (abs(s->startingTime() - (k + 1)*h) < 1e-12)
      {
        k++;
        dataPlot(k, 0) =  s->startingTime();
        dataPlot(k, 1) = (*q)(0);
        dataPlot(k, 2) = (*q)(1);
        dataPlot(k, 3) = (*q)(2);
        dataPlot(k, 4) = (*q)(3);
        dataPlot(k, 5) = (*q)(4);
        dataPlot(k, 6) = (*q)(5);
        dataPlot(k, 7) = (*v)(0);
        dataPlot(k, 8) = (*v)(1);
        dataPlot(k, 9) = (*v)(2);
        dataPlot(k, 10) = (*v)(3);
        dataPlot(k, 11) = (*v)(4);
        dataPlot(k, 12) = (*v)(5);
      }
    }

    // --- Output files ---
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in RX90Control" << endl;
  }
  cout << "Computation Time " << time.elapsed()  << endl;
}

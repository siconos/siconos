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

/*!\file PID.cpp
  \brief \ref EMPID - C++ input file -
  O. Huber.

  Simple PID example.
  The controlled plant is a double integrator
*/

#include "SiconosKernel.hpp"
#include "PID.hpp"
#include "LinearSensor.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for the system
    double t0 = 0;                   // initial computation time
    double T = 100;                  // final computation time
    double h = 0.05;                // time step
    double position_init = 10;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for Moreau integrator
    double xFinal = 0;              // final value
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;

    SP::SiconosMatrix A(new SimpleMatrix(nDof, nDof));
    A->zero();
    (*A)(0, 1) = 1;
    (*A)(1, 0) = -1;

    SP::SiconosVector B(new SiconosVector(nDof));
    B->zero();


    // -- Initial positions and velocities --
    SP::SiconosVector x0(new SiconosVector(nDof));
    (*x0)(0) = position_init;
    (*x0)(1) = velocity_init;

    // -- The dynamical system --
    SP::FirstOrderLinearTIDS doubleIntegrator(new FirstOrderLinearTIDS(x0, A, B));

    // -------------
    // --- Model ---
    // -------------
    SP::Model process(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    process->nonSmoothDynamicalSystem()->insertDynamicalSystem(doubleIntegrator);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::ZeroOrderHoldOSI OSI(new ZeroOrderHoldOSI(doubleIntegrator));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::TimeDiscretisation tSensor(new TimeDiscretisation(t0, 2*h));
    SP::TimeDiscretisation tActuator(new TimeDiscretisation(t0, 2*h));
    SP::TimeDiscretisation tObserver(new TimeDiscretisation(t0, 2*h));

    // -- (3) Simulation setup with (1) (2)
    SP::TimeStepping processSimulation(new TimeStepping(t, 0));
    processSimulation->insertIntegrator(OSI);

    // define the control manager
    SP::ControlManager control(new ControlManager(processSimulation));

    // use a controlSensor
    SP::SimpleMatrix C(new SimpleMatrix(1, 2, 0));
    (*C)(0, 0) = 1;
    SP::LinearSensor sens(new LinearSensor(doubleIntegrator, C));
    control->addSensorPtr(sens, tSensor);

    // add the Observer
    SP::SiconosMatrix L(new SimpleMatrix(2, 1));
    (*L)(0, 0) = -8.71455866;
    (*L)(1, 0) = -24.02084375;
    SP::SiconosVector xHat0(new SiconosVector(2));
    (*xHat0)(0) = -10;
    (*xHat0)(1) = -5;
    SP::Observer obs(new LuenbergerObserver(sens, *xHat0, C, L));
    control->addObserverPtr(obs, tObserver);
    // add the PID controller
    SP::SiconosVector K(new SiconosVector(3, 0));
    (*K)(0) = .25;
    (*K)(1) = .125;
    (*K)(2) = 2;
    SP::PID act = std11::static_pointer_cast<PID>
                                 (control->addActuator(PID_, tActuator, sens));

    // To store the nextEvent
    SP::Event currentEvent;
    SP::Event nextEvent;

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    // Initialize the model and the controlManager
    process->initialize(processSimulation);
    control->initialize(*process);
    act->setRef(xFinal);
    act->setK(*K);

    SP::EventsManager eventsManager = processSimulation->eventsManager();
    unsigned int N = ceil((T - t0) / h + 10); // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 7;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SiconosVector& xProc = *doubleIntegrator->x();
    const SiconosVector& u = act->u();
    SiconosVector& e = *obs->e();
    SiconosVector& xHat = *obs->xHat();
    dataPlot(0, 0) = process->t0();
    dataPlot(0, 1) = (xProc)(0);
    dataPlot(0, 2) = (xProc)(1);
    dataPlot(0, 3) = (u)(0);
    dataPlot(0, 4) = (e)(0);
    dataPlot(0, 5) = (xHat)(0);
    dataPlot(0, 6) = (xHat)(1);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (processSimulation->hasNextEvent())
    {
      nextEvent = eventsManager->nextEvent();
      // --- Get values to be plotted ---
      // the following check prevents saving the same data multiple times
      // XXX what happends if we have NS Events ?
      if (nextEvent->getType() == TD_EVENT)
        processSimulation->computeOneStep();

      processSimulation->nextStep();
      currentEvent = eventsManager->currentEvent();
      if (currentEvent->getType() == OBSERVER_EVENT)
      {
        dataPlot(k, 0) =  processSimulation->startingTime();
        dataPlot(k, 1) = (xProc)(0);
        dataPlot(k, 2) = (xProc)(1);
        dataPlot(k, 3) = (u)(0);
        dataPlot(k, 4) = (e)(0);
        dataPlot(k, 5) = (xHat)(0);
        dataPlot(k, 6) = (xHat)(1);
        ++show_progress;
        k++;
      }
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("LuenbergerObserver.dat", "ascii", dataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("LuenbergerObserver.ref", "ascii", dataPlotRef);

    std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;

    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
      std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }
    else
    {
      return 0;
    }

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in PID.cpp" << endl;
  }



}

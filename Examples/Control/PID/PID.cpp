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
#include "SampledPIDActuator.hpp"
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

    SP::SiconosVector B(new SimpleVector(nDof));
    B->zero();


    // -- Initial positions and velocities --
    SP::SimpleVector x0(new SimpleVector(nDof));
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
    SP::Moreau OSI(new Moreau(doubleIntegrator, theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::TimeDiscretisation tSensor(new TimeDiscretisation(t0, h));
    SP::TimeDiscretisation tActuator(new TimeDiscretisation(t0, h));

    // -- (3) Simulation setup with (1) (2)
    SP::TimeStepping s(new TimeStepping(t, 0));
    s->insertIntegrator(OSI);

    // define the control manager
    SP::ControlManager control(new ControlManager(process));

    // use a controlSensor
    SP::SimpleMatrix C(new SimpleMatrix(1, 2, 0));
    (*C)(0, 0) = 1;
    SP::SimpleMatrix D(new SimpleMatrix(1, 2, 0));
    SP::LinearSensor sens(new LinearSensor(tSensor, doubleIntegrator, C, D));
    control->addSensorPtr(sens);
    // add the PID controller
    SP::SimpleVector K(new SimpleVector(3, 0));
    (*K)(0) = .25;
    (*K)(1) = .125;
    (*K)(2) = 2;
    SP::SampledPIDActuator act = static_pointer_cast<SampledPIDActuator>(control->addActuator(100, tActuator));
    act->addSensorPtr(sens);

    // To store the nextEvent
    SP::Event nextEvent;

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    // Initialize the model and the controlManager
    process->initialize(s);
    control->initialize();
    act->setRef(xFinal);
    act->setK(*K);

    SP::EventsManager eventsManager = s->eventsManager();
    unsigned int N = ceil((T - t0) / h + 10); // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 5;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector xProc = doubleIntegrator->x();
    SP::SiconosVector u = doubleIntegrator->b();
    dataPlot(0, 0) = process->t0();
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    dataPlot(0, 3) = (*u)(0);
    dataPlot(0, 4) = (*u)(1);
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

    while (s->nextTime() < T)
    {
      s->computeOneStep();
      nextEvent = eventsManager->followingEvent(eventsManager->currentEvent());
      // --- Get values to be plotted ---
      // the following check prevents saving the same data multiple times
      // XXX what happends if we have NS Events ?
      if (nextEvent->getType() == 1)
      {
        dataPlot(k, 0) =  s->nextTime();
        dataPlot(k, 1) = (*xProc)(0);
        dataPlot(k, 2) = (*xProc)(1);
        dataPlot(k, 3) = (*u)(0);
        dataPlot(k, 4) = (*u)(1);
        ++show_progress;
        k++;
      }
      s->nextStep();
    }
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    dataPlot.resize(k, outputSize);
    io.write(dataPlot, "noDim");
    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix ref("result.ref", "ascii");
    ref.read(dataPlotRef);

    if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    {
      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
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

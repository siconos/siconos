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
#include "sampledPIDActuator.hpp"
#include "linearSensor.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for the system
    double t0 = 0;                   // initial computation time
    double T = 10;                  // final computation time
    double h = 0.005;                // time step
    double position_init = 10;      // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double theta = 0.5;              // theta for Moreau integrator

    SP::SimpleVector xFinal(new SimpleVector(2, 0));
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
    SP::SimpleMatrix C(new SimpleMatrix(2, 2));
    C->eye();
    SP::SimpleMatrix D(new SimpleMatrix(2, 2, 0));
    SP::linearSensor sens(new linearSensor(100, tSensor, process, C, D));
    control->addSensorPtr(sens);
    // add the PID controller
    SP::SimpleVector K(new SimpleVector(3, 0));
    (*K)(0) = 1.1;
    (*K)(1) = .5;
    SP::sampledPIDActuator act = static_pointer_cast<sampledPIDActuator>(control->addActuator(100, tActuator));
    act->addSensorPtr(sens);

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    // Initialize the model and the controlManager
    process->initialize(s);
    control->initialize();
    act->setRef(*xFinal);
    act->setK(*K);

    SP::EventsManager eventsManager = s->eventsManager();
    //    eventsManager->display();
    int N = (int)(3 * (T - t0) / h); // Number of time steps
    //    int N = 10000;
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 3;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector xProc = doubleIntegrator->x();

    dataPlot(0, 0) = process->t0();
    dataPlot(0, 1) = (*xProc)(0);
    dataPlot(0, 2) = (*xProc)(1);
    //    dataPlot(0, 3) = (*yProc)(0);
    //    dataPlot(0, 4) = (*lambda)(0);
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
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      s->nextStep();
      ++show_progress;
      k++;
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

/* Siconos-sample , Copyright INRIA 2005-2012.
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

/* !\file SMCExampleImplicit.cpp
  \brief Two independent systems of dimension one controlled to slide
  on \f$x = 0\f$. An implicit scheme is used
  O. Huber
  */

#include "SiconosKernel.hpp"
using namespace std;

// main program
int main(int argc, char* argv[])
{

double t0 = 0.0;
double T = 10;
double h = 1.0e-4;
double hControl = 1.0e-2;

unsigned int N = ceil((T-t0)/h + 10);

unsigned int n = 4;
unsigned int outputSize = 1+n+2;

/*
A = [[-0.1452, 0, 0, 0.1611],
     [-0.0274, 0, 0, 1.879],
     [1, 0, 0, 0],
     [0, 1, 0, 0]]

Brel = [[1.452], [0.274], [0], [0]]
Csurface = np.array(Brel).T
*/
  SP::SiconosMatrix A(new SimpleMatrix(n, n, 0));
  (*A)(0, 0) = -0.1452;
  (*A)(0, 3) = 0.1611;
  (*A)(1, 0) = -0.0274;
  (*A)(1, 3) = 1.879;
  (*A)(2, 0) = 1;
  (*A)(3, 1) = 1;

  SP::SiconosVector x0(new SiconosVector(n));
  (*x0)(0) = -4;
  (*x0)(1) = 2.21;
  (*x0)(2) = 3;
  (*x0)(3) = -2.5;

  // For the Sensor
  SP::SimpleMatrix sensorC(new SimpleMatrix(n, n));
  sensorC->eye();
  // For the Actuator
  SP::SimpleMatrix Csurface(new SimpleMatrix(1, n));
  (*Csurface)(0, 0) = 1.452;
  (*Csurface)(0, 1) = 0.274;
  (*Csurface)(0, 2) = 0;
  (*Csurface)(0, 3) = 0;
  SP::SimpleMatrix Brel(new SimpleMatrix(n, 1));
  (*Brel)(0, 0) = 1.452;
  (*Brel)(1, 0) = 0.274;
  (*Brel)(2, 0) = 0;
  (*Brel)(3, 0) = 0;

  SP::SimpleMatrix Drel(new SimpleMatrix(1, 1, 0));
  // Dynamical Systems
  SP::FirstOrderLinearDS processDS(new FirstOrderLinearDS(x0, A));
  processDS->setComputebFunction("RelayPluginUnperturbed.so", "computeB");

  // -------------
  // --- Model process ---
  // -------------
  SP::Model process(new Model(t0, T));
  process->nonSmoothDynamicalSystem()->insertDynamicalSystem(processDS);

  // ------------------
  // --- Simulation ---
  // ------------------
  // TimeDiscretisation
  SP::TimeDiscretisation processTD(new TimeDiscretisation(t0, h));
  SP::TimeDiscretisation tSensor(new TimeDiscretisation(t0, hControl));
  SP::TimeDiscretisation tActuator(new TimeDiscretisation(t0, hControl));
  // == Creation of the Simulation ==
  SP::TimeStepping processSimulation(new TimeStepping(processTD, 0));
  processSimulation->setName("plant simulation");
  // -- OneStepIntegrators --
//  SP::OneStepIntegrator processIntegrator(new ZeroOrderHold(processDS, theta));
  SP::OneStepIntegrator processIntegrator(new ZeroOrderHold(processDS));
  processSimulation->insertIntegrator(processIntegrator);

  // Control stuff
  SP::ControlManager control(new ControlManager(process));
  // use a controlSensor
  SP::LinearSensor sens(new LinearSensor(tSensor, processDS, sensorC));
  control->addSensorPtr(sens);
  // add the sliding mode controller
  LinearSMC& act = *std11::static_pointer_cast<LinearSMC>(control->addActuator
                   (LINEAR_SMC, tActuator));
  act.setCsurfacePtr(Csurface);
  act.setBPtr(Brel);
  act.setDPtr(Drel);
  act.addSensorPtr(sens);
  act.setTheta(0);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ..." << endl << endl;
  // initialise the process and the ControlManager
  process->initialize(processSimulation);
  control->initialize();

  // --- Get the values to be plotted ---
  SiconosVector& xProc = *processDS->x();
  // Save data in a matrix dataPlot
  SimpleMatrix dataPlot(N, outputSize);

  dataPlot(0, 0) = process->t0(); // Initial time of the model
  for(unsigned int i = 0; i < n; i++)
  {
    dataPlot(0, i+1) = xProc(i);
  }

  EventsManager& eventsManager = *processSimulation->eventsManager();

  // ==== Simulation loop =====
  cout << "====> Start computation ... " << endl << endl;
  unsigned int k = 0; // Current step
  // Simulation loop
  boost::progress_display show_progress(N);
  boost::timer time;
  time.restart();
  while (processSimulation->hasNextEvent())
  {
    Event& nextEvent = *eventsManager.nextEvent();
    if (nextEvent.getType() == TD_EVENT)
    {
      k++;
      processSimulation->computeOneStep();
      dataPlot(k, 0) = processSimulation->nextTime();
      for(unsigned int i = 0; i < n; i++)
      {
          dataPlot(k, i+1) = xProc(i);
      }
      dataPlot(k, n+1) = act.ueq()(0);
      dataPlot(k, n+2) = act.us()(0);
      ++show_progress;
    }
    processSimulation->nextStep();
  }
  cout << endl << "Computation Time " << time.elapsed()  << endl;
  eventsManager.display();

  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  dataPlot.resize(k, outputSize);
  ioMatrix::write("SMCExampleExplicitImplicit.dat", "ascii", dataPlot, "noDim");

#if 0
  // Comparison with a reference file
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  ioMatrix::read("SMCExampleImplicit.ref", "ascii", dataPlotRef);
  std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;
  if ((dataPlot - dataPlotRef).normInf() > 1e-12)
  {
    std::cout << "Warning. The results is rather different from the reference file." << std::endl;
    (dataPlot - dataPlotRef).display();
    return 1;
  }
#endif
}


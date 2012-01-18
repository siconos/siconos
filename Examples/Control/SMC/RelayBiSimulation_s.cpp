
#include "SiconosKernel.hpp"
using namespace std;

// main program
int main(int argc, char* argv[])
{
  // User-defined parameters
  unsigned int ndof = 2;          // Number of degrees of freedom of your system
  double t0 = 0.0;                // Starting time
  double T = 1;                   // Total simulation time
  double h = 1.0e-4;              // Time step for simulation
  double hControl = 1.0e-2;       // Time step for control
  double Xinit = 1.0;
  double theta = 0.5;

  if (h > hControl)
  {
    RuntimeException::selfThrow("hControl must be bigger than h");
  }

  // ================= Creation of the model =======================
  // Steps:
  // - create a Dynamical System
  // - add a Simulation to the model

  // -------------------------
  // --- Dynamical systems ---
  // -------------------------

  // First System:
  // dx/dt = Ax + u(t) + r
  // x(0) = x0
  // Note: r = Blambda, B defines in relation below.

  // Matrix declaration
  SP::SiconosMatrix A(new SimpleMatrix(ndof, ndof, 0));
  SP::SiconosVector x0(new SimpleVector(ndof));
  (*x0)(0) = Xinit;
  (*x0)(1) = -Xinit;
  SP::SimpleMatrix sensorC(new SimpleMatrix(2, 2));
  sensorC->eye();
  SP::SimpleMatrix sensorD(new SimpleMatrix(2, 2, 0));
  SP::SimpleMatrix Csurface(new SimpleMatrix(2, 2));
  Csurface->eye();

  // Dynamical Systems
  SP::ControlFirstOrderLinearDS controlProcess(new ControlFirstOrderLinearDS(t0, T, h, x0, A));
  SP::FirstOrderLinearDS processDS = controlProcess->processDS();
  processDS->setComputebFunction("RelayPlugin.so", "computeB");
  controlProcess->initialize();
  SP::Model controlModel = controlProcess->model();

  // TimeDiscretisation
  SP::TimeDiscretisation tSensor(new TimeDiscretisation(t0, hControl));
  SP::TimeDiscretisation tActuator(new TimeDiscretisation(t0, hControl));

  // Control stuff
  SP::ControlManager control(new ControlManager(controlModel));
  // use a controlSensor
  SP::linearSensor sens(new linearSensor(100, tSensor, controlModel, sensorC, sensorD));
  control->addSensorPtr(sens);
  // add the sliding mode controller
  SP::linearSMC act = static_pointer_cast<linearSMC>(control->addActuator(101, tActuator));
  act->addSensorPtr(sens);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ..." << endl << endl;
  // initialise the process and the ControlManager
  control->initialize();
  // Get the simulation
  SP::TimeStepping processSimulation = controlProcess->simulation();
  // Only now we can add the surface
  act->setCsurfacePtr(Csurface);

  // --- Get the values to be plotted ---
  unsigned int outputSize = 3; // number of required data
  int N = (int)((T - t0) / h) + 10; // Number of time steps

  SP::SiconosVector xProc = processDS->x();
  // Save data in a matrix dataPlot
  SimpleMatrix dataPlot(N, outputSize);
  dataPlot(0, 0) = controlModel->t0(); // Initial time of the model
  dataPlot(0, 1) = (*xProc)(0);
  dataPlot(0, 2) = (*xProc)(1);

  SP::EventsManager eventsManager = processSimulation->eventsManager();

  // ==== Simulation loop =====
  cout << "====> Start computation ... " << endl << endl;
  int k = 0; // Current step
  // Simulation loop
  boost::progress_display show_progress(N);
  boost::timer time;
  time.restart();
  SP::Event nextEvent;
  while (processSimulation->nextTime() < T)
  {
    processSimulation->computeOneStep();
    nextEvent = eventsManager->followingEvent(eventsManager->currentEvent());
    if (nextEvent->getType() == 1)
    {
      k++;
      dataPlot(k, 0) = processSimulation->nextTime();
      dataPlot(k, 1) = (*xProc)(0);
      dataPlot(k, 2) = (*xProc)(1);
      ++show_progress;
    }
    processSimulation->nextStep();
  }
  cout << endl << "Computation Time " << time.elapsed()  << endl;

  // --- Output files ---
  cout << "====> Output file writing ..." << endl;
  ioMatrix io("RelayBiSimulation_s.dat", "ascii");
  dataPlot.resize(k, outputSize);
  io.write(dataPlot, "noDim");

  // Comparison with a reference file
  SimpleMatrix dataPlotRef(dataPlot);
  dataPlotRef.zero();
  ioMatrix ref("RelayBiSimulation.ref", "ascii");
  ref.read(dataPlotRef);
  std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;

  if ((dataPlot - dataPlotRef).normInf() > 1e-12)
  {
    std::cout << "Warning. The results is rather different from the reference file." << std::endl;
    return 1;
  }

}


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
  // For the DynamicalSystem
  SP::SiconosMatrix A(new SimpleMatrix(ndof, ndof, 0));
  SP::SiconosVector x0(new SimpleVector(ndof));
  (*x0)(0) = Xinit;
  (*x0)(1) = -Xinit;
  // For the Sensor
  SP::SimpleMatrix sensorC(new SimpleMatrix(2, 2));
  sensorC->eye();
  SP::SimpleMatrix sensorD(new SimpleMatrix(2, 2, 0));
  // For the Actuator
  SP::SimpleMatrix Csurface(new SimpleMatrix(2, 2));
  Csurface->eye();
  SP::SimpleMatrix Brel(new SimpleMatrix(2, 2));
  Brel->eye();
  *(Brel) *= 2;
  SP::SimpleMatrix Drel(new SimpleMatrix(2, 2, 0));

  // Dynamical Systems
  SP::ControlFirstOrderLinearDS controlProcess(new ControlFirstOrderLinearDS(t0, T, h, x0, A));
  SP::FirstOrderLinearDS processDS = controlProcess->processDS();
  processDS->setComputebFunction("RelayPlugin.so", "computeB");
  controlProcess->initialize();

  // TimeDiscretisation
  SP::TimeDiscretisation tSensor(new TimeDiscretisation(t0, hControl));
  SP::TimeDiscretisation tActuator(new TimeDiscretisation(t0, hControl));

  // Control stuff
  SP::ControlManager control = controlProcess->CM();
  // use a controlSensor
  SP::LinearSensor sens(new LinearSensor(tSensor, processDS, sensorC, sensorD));
  control->addSensorPtr(sens);
  // add the sliding mode controller
  SP::LinearSMC act = static_pointer_cast<LinearSMC>(control->addActuator(101, tActuator));
  act->setCsurfacePtr(Csurface);
  act->setBPtr(Brel);
  act->setDPtr(Drel);
  act->addSensorPtr(sens);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  cout << "====> Simulation initialisation ..." << endl << endl;
  // initialise the ControlManager
  control->initialize();

  // ==== Simulation loop =====
  cout << "====> Start computation ... " << endl << endl;
  controlProcess->run();
  // Simulation loop

  cout << endl << "Computation Time " << controlProcess->elapsedTime()  << endl;
  SP::SimpleMatrix dataPlot = controlProcess->data();
  // We are only interested in the state
  dataPlot->resize(dataPlot->size(0), 3);

  // ==== Output files ====
  cout << "====> Output file writing ..." << endl;
  ioMatrix io("RelayBiSimulation_s.dat", "ascii");
  io.write(*dataPlot, "noDim");

  // Comparison with a reference file
  SimpleMatrix dataPlotRef(*dataPlot);
  dataPlotRef.zero();
  ioMatrix ref("RelayBiSimulation-SMC.ref", "ascii");
  ref.read(dataPlotRef);
  std::cout << (*dataPlot - dataPlotRef).normInf() << std::endl;

  if ((*dataPlot - dataPlotRef).normInf() > 1e-12)
  {
    std::cout << "Warning. The results is rather different from the reference file." << std::endl;
    return 1;
  }

}

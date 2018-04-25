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

/* !\file ZeroOrderHold.cpp
  */

#include "SiconosKernel.hpp"
#include "SiconosControl.hpp"
using namespace std;
class MyDS : public FirstOrderLinearDS
{
public:
  MyDS(SP::SiconosVector x0) : FirstOrderLinearDS(x0)
  {
    _b.reset(new SiconosVector(x0->size()));
    _A.reset(new SimpleMatrix(x0->size(), x0->size(), 0));
  };
  MyDS(const MyDS & myds) : FirstOrderLinearDS(myds)
  {
    // copy constructor
    std::cout << "copy constructor" << std::endl;
  };

  

  
  void computeA(double time)
  {
    printf("computeA\n");
    double t = sin(50 * time);
    _A->eye();
    *_A  = t * *_A;
    _A->display();
  };
  
  void computeb(double time)
  {
    printf("computeB\n");
    double t = sin(50 * time);
    _b->setValue(0,t);
    _b->setValue(1,-t);
    printf("b[0] = %g, b[1] = %g\n", _b->getValue(0), _b->getValue(1));
  };
};
// main program
int main(int argc, char* argv[])
{
  // User-defined parameters
  unsigned int ndof = 2;          // Number of degrees of freedom of your system
  double t0 = 0.0;                // Starting time
  double T = 1;                   // Total simulation time
  double h = 1.0e-3;              // Time step for simulation
  double Xinit = 1.0;

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
  SP::SiconosVector x0(new SiconosVector(ndof));
  (*x0)(0) = Xinit;
  (*x0)(1) = -Xinit;
  SP::SimpleMatrix sensorC(new SimpleMatrix(2, 2));
  sensorC->eye();
  SP::SimpleMatrix sensorD(new SimpleMatrix(2, 2, 0));
  SP::SimpleMatrix Csurface(new SimpleMatrix(1, 2, 0));
  (*Csurface)(0, 1) = 1;
  SP::SimpleMatrix Brel(new SimpleMatrix(2, 1, 0));
  (*Brel)(1, 0) = 2;

  // Dynamical Systems
  SP::FirstOrderLinearDS processDS(new MyDS(x0));


  SP::NonSmoothDynamicalSystem nsds(new NonSmoothDynamicalSystem(t0, T));
  nsds->insertDynamicalSystem(processDS);

  // TimeDiscretisation
  SP::TimeDiscretisation td(new TimeDiscretisation(t0, h));
  // == Creation of the Simulation ==
  SP::TimeStepping s(new TimeStepping(nsds, td));

  // -- OneStepIntegrators --
  SP::ZeroOrderHoldOSI myIntegrator(new ZeroOrderHoldOSI());
  s->insertIntegrator(myIntegrator);

  // =========================== End of model definition ===========================

  // ================================= Computation =================================

  // --- Simulation initialization ---

  // --- Get the values to be plotted ---
  unsigned int outputSize = 3; // number of required data
  unsigned int N = ceil((T - t0) / h); // Number of time steps

  SimpleMatrix dataPlot(N, outputSize);

  SP::SiconosVector xProc = processDS->x();


  // -> saved in a matrix dataPlot
  dataPlot(0, 0) = nsds->t0(); // Initial time of the model
  dataPlot(0, 1) = (*xProc)(0);
  dataPlot(0, 2) = (*xProc)(1);

  // ==== Simulation loop =====
  cout << "====> Start computation ... " << endl << endl;

  // *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
  unsigned int k = 0; // Current step

  // Simulation loop
  boost::timer time;
  time.restart();
  while (k < N - 1)
  {
    k++;
    //  *z = *(myProcessInteraction->y(0)->getVectorPtr(0));
    s->computeOneStep();
    dataPlot(k, 0) = s->nextTime();
    dataPlot(k, 1) = (*xProc)(0);
    dataPlot(k, 2) = (*xProc)(1);
    s->nextStep();
  }
  cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
  cout << "Computation Time " << time.elapsed()  << endl;


  // --- Output files ---
  cout << "====> Output file writing ..." << endl;

  ioMatrix::write("ZeroOrderHold.dat", "ascii", dataPlot, "noDim");

  // // Comparison with a reference file
  // SimpleMatrix dataPlotRef(dataPlot);
  // dataPlotRef.zero();
  // ioMatrix::read("SMCExampleImplicit.ref", "ascii", dataPlotRef);
  // std::cout << (dataPlot - dataPlotRef).normInf() << std::endl;

  // if ((dataPlot - dataPlotRef).normInf() > 1e-12)
  // {
  //   std::cout << "Warning. The results is rather different from the reference file." << std::endl;
  //   return 1;
  // }

}

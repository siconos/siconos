/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2008.
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

/*!\file DryFriction.cpp
\brief \ref EMDryFriction - C++/XML input file, Time-Stepping - F. Perignon.

Description of the model with XML input.
Simulation with a Time-Stepping scheme.
*/
#include "SiconosKernel.h"

using namespace std;
int main(int argc, char* argv[])
{
  try
  {

    // --- Model loading from xml file ---
    Model oscillator("./DryFriction.xml");
    cout << "\n *** DryFriction.xml file loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(oscillator.getSimulationPtr());
    cout << "simulation initialization" << endl;
    s->initialize();
    cout << "\n **** the simulation is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = 0;
    int N = t->getNSteps(); // Number of time steps

    //t->display();

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 5);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) = k * t->getH();
    // state q for the first dynamical system (ball)
    LagrangianDS* oscillo = static_cast<LagrangianDS*>(oscillator.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    dataPlot(k, 1) = (oscillo->getQ())(0);
    // velocity for the oscillo
    dataPlot(k, 2) = (oscillo->getVelocity())(0);
    dataPlot(k, 3) = (oscillator.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
    dataPlot(k, 4) = (oscillator.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(1);

    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    cout << "Computation ... " << endl;
    // --- Time loop  ---
    while (k < N)
    {
      // get current time step
      k++;
      //  cout << " Pas " << k;
      // solve ...
      s->computeOneStep();
      // --- Get values to be plotted ---
      //time
      dataPlot(k, 0) = k * t->getH();
      // Oscillo: state q
      dataPlot(k, 1) = (oscillo->getQ())(0);
      // Oscillo: velocity
      dataPlot(k, 2) = (oscillo->getVelocity())(0);
      dataPlot(k, 3) = (oscillator.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      dataPlot(k, 4) = (oscillator.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(1);
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
    }

    // --- elapsed time computing ---
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    // Xml output
    //  oscillator.saveToXMLFile("./BouncingOscillo_TIDS.xml.output");
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/DryFriction\'" << endl;
  }
}

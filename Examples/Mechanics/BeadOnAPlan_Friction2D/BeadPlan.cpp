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
#include "SiconosKernel.hpp"
using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // --- Model loading from xml file ---
    SP::Model oscillator(new Model("./BeadPlan.xml"));
    cout << "\n *** BeadPlan.xml file loaded ***" << endl;
    oscillator->initialize();

    // --- Get and initialize the simulation ---
    SP::TimeStepping s = boost::static_pointer_cast<TimeStepping>(oscillator->getSimulationPtr());

    // --- Get the time discretisation scheme ---
    SP::TimeDiscretisation t = s->getTimeDiscretisationPtr();
    int k = 0;
    double t0 = oscillator->getT0();
    double T = oscillator->getFinalT();
    double h = s->getTimeStep();
    int N = (int)((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 3);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    SP::LagrangianDS oscillo = boost::static_pointer_cast<LagrangianDS> (oscillator->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtrNumber(1));
    SP::SiconosVector q = oscillo->getQPtr();
    SP::SiconosVector v = oscillo->getVelocityPtr();
    SP::SiconosVector p = oscillo->getPPtr(2);

    dataPlot(k, 0) = t0;
    dataPlot(0, 1) = (*q)(0);
    // velocity for the oscillo
    dataPlot(0, 2) = (*v)(0);

    cout << "Computation ... " << endl;
    // --- Time loop  ---
    while (s->getNextTime() <= oscillator->getFinalT())
    {
      // solve ...
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->getNextTime();
      dataPlot(k, 1) = (oscillo->getQ())(0);
      dataPlot(k, 2) = (oscillo->getVelocity())(0);
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
      k++;
    }

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
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

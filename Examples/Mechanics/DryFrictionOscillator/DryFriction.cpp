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

/*!\file DryFriction.cpp
\brief \ref EMDryFriction - C++/XML input file, Time-Stepping - F. Perignon.

Description of the model with XML input.
Simulation with a Time-Stepping scheme.
*/
#include "SiconosKernel.hpp"

using namespace std;
int main(int argc, char* argv[])
{
  try
  {

    // --- Model loading from xml file ---
    SP::Model oscillator(new Model("./DryFriction.xml"));
    cout << "\n *** DryFriction.xml file loaded ***" << endl;
    oscillator->initialize();

    // --- Get the simulation ---
    SP::TimeStepping s = boost::static_pointer_cast<TimeStepping>(oscillator->simulation());
    // --- Get the time discretisation scheme ---
    SP::TimeDiscretisation t = s->timeDiscretisation();
    int k = 0;
    double T = oscillator->finalT();
    double t0 = oscillator->t0();
    double h = s->timeStep();
    int N = ceil((T - t0) / h);
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 5);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) = t0;
    // state q for the first dynamical system (ball)
    SP::LagrangianDS oscillo = boost::static_pointer_cast<LagrangianDS> (oscillator->nonSmoothDynamicalSystem()->dynamicalSystemNumber(1));
    dataPlot(k, 1) = ((*oscillo->q()))(0);
    // velocity for the oscillo
    dataPlot(k, 2) = ((*oscillo->velocity()))(0);
    dataPlot(k, 3) = (*oscillator->nonSmoothDynamicalSystem()->topology()->interactions()->getPtr(1)->lambda(1))(0);
    dataPlot(k, 4) = (*oscillator->nonSmoothDynamicalSystem()->topology()->interactions()->getPtr(1)->lambda(1))(1);

    // --- Compute elapsed time ---
    boost::timer tt;
    tt.restart();
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
      dataPlot(k, 0) = s->nextTime();
      // Oscillo: state q
      dataPlot(k, 1) = ((*oscillo->q()))(0);
      // Oscillo: velocity
      dataPlot(k, 2) = ((*oscillo->velocity()))(0);

      dataPlot(k, 3) = (*oscillator->nonSmoothDynamicalSystem()->topology()->interactions()->getPtr(1)->lambda(1))(0);
      dataPlot(k, 4) = (*oscillator->nonSmoothDynamicalSystem()->topology()->interactions()->getPtr(1)->lambda(1))(1);
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
    }

    // --- elapsed time computing ---
    cout << "time = " << tt.elapsed() << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

    // Xml output
    //  oscillator->saveToXMLFile("./BouncingOscillo_TIDS.xml.output");
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

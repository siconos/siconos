/* Siconos-sample version 2.0.1, Copyright INRIA 2005-2006.
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
#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // --- Model loading from xml file ---
    Model uBeads("./ThreeBeadsColumn.xml");
    cout << "\n *** ThreeBeadsColumn.xml loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(uBeads.getSimulationPtr());
    s->initialize();
    cout << "\n **** the simulation is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = 0;
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 7);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) = k * t->getH();

    // state q and velocity for the first dynamical system
    LagrangianLinearTIDS* bead = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    dataPlot(k, 1) = (bead->getQ())(0);
    dataPlot(k, 2) = (bead->getVelocity())(0);

    // state q and velocity for the second dynamical system
    LagrangianLinearTIDS* bead2 = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(1));
    dataPlot(k, 3) = (bead2->getQ())(0);
    dataPlot(k, 4) = (bead2->getVelocity())(0);

    // state q and velocity for the third dynamical system
    LagrangianLinearTIDS* bead3 = static_cast<LagrangianLinearTIDS*>(uBeads.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(2));
    dataPlot(k, 5) = (bead3->getQ())(0);
    dataPlot(k, 6) = (bead3->getVelocity())(0);

    cout << " Computation ... " << endl;
    while (k < N)
    {
      // get current time step
      k++;
      // solve ...
      s->computeOneStep();

      //cout<<"Iteration: "<<k<<endl;
      // --- Get values to be plotted ---

      dataPlot(k, 0) = k * t->getH();
      dataPlot(k, 1) = (bead->getQ())(0);
      dataPlot(k, 2) = (bead->getVelocity())(0);
      dataPlot(k, 3) = (bead2->getQ())(0);
      dataPlot(k, 4) = (bead2->getVelocity())(0);
      dataPlot(k, 5) = (bead3->getQ())(0);
      dataPlot(k, 6) = (bead3->getVelocity())(0);
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();

    }
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
    cout << "End of computation - Number of iterations  done: " << k << endl;

    //    uBeads.saveToXMLFile("./ThreeBeadsColumn.xml.output");

  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/ThreeBeadsColumn\'" << endl;
  }
}

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
#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // --- Model loading from xml file ---
    Model * BeadsModel = new Model("./ThreeBeadsColumn.xml");

    cout << "\n *** ThreeBeadsColumn.xml loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping * s = static_cast<TimeStepping*>(BeadsModel->getSimulationPtr());
    s->initialize();

    cout << "\n **** the simulation is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = 0;
    double t0 = BeadsModel->getT0();
    double T = BeadsModel->getFinalT();
    double h = s->getTimeStep();
    int N = (int)((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 7);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) = BeadsModel->getT0();

    // state q and velocity for the first dynamical system
    LagrangianLinearTIDS* bead = static_cast<LagrangianLinearTIDS*>(BeadsModel->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    SiconosVector * q1 = bead->getQPtr();
    SiconosVector * v1 = bead->getVelocityPtr();
    dataPlot(k, 1) = (*q1)(0);
    dataPlot(k, 2) = (*v1)(0);

    // state q and velocity for the second dynamical system
    LagrangianLinearTIDS* bead2 = static_cast<LagrangianLinearTIDS*>(BeadsModel->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(1));
    SiconosVector * q2 = bead2->getQPtr();
    SiconosVector * v2 = bead2->getVelocityPtr();
    dataPlot(k, 3) = (*q2)(0);
    dataPlot(k, 4) = (*v2)(0);

    // state q and velocity for the third dynamical system
    LagrangianLinearTIDS* bead3 = static_cast<LagrangianLinearTIDS*>(BeadsModel->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(2));
    SiconosVector * q3 = bead3->getQPtr();
    SiconosVector * v3 = bead3->getVelocityPtr();
    dataPlot(k, 5) = (*q3)(0);
    dataPlot(k, 6) = (*v3)(0);

    cout << " Computation ... " << endl;
    while (s->getNextTime() <= BeadsModel->getFinalT())
    {
      // solve ...
      s->computeOneStep();

      //cout<<"Iteration: "<<k<<endl;
      // --- Get values to be plotted ---

      dataPlot(k, 0) = s->getStartingTime();
      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*v1)(0);
      dataPlot(k, 3) = (*q2)(0);
      dataPlot(k, 4) = (*v2)(0);
      dataPlot(k, 5) = (*q3)(0);
      dataPlot(k, 6) = (*v3)(0);
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
      k++;

    }
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");
    cout << "End of computation - Number of iterations  done: " << k << endl;

    //    BeadsModel.saveToXMLFile("./ThreeBeadsColumn.xml.output");

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

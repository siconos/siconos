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
    SP::Model BeadsModel(new Model("./ThreeBeadsColumn.xml"));

    cout << "\n *** ThreeBeadsColumn.xml loaded ***" << endl;
    BeadsModel->initialize();

    // --- Get and initialize the simulation ---
    SP::TimeStepping s = boost::static_pointer_cast<TimeStepping>(BeadsModel->simulation());

    cout << "\n **** the simulation is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    SP::TimeDiscretisation t = s->timeDiscretisation();
    int k = 0;
    double t0 = BeadsModel->t0();
    double T = BeadsModel->finalT();
    double h = s->timeStep();
    int N = (int)((T - t0) / h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 7);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) = BeadsModel->t0();

    // state q and velocity for the first dynamical system
    SP::LagrangianLinearTIDS bead = boost::static_pointer_cast<LagrangianLinearTIDS> (BeadsModel->nonSmoothDynamicalSystem()->dynamicalSystemNumber(1));
    SP::SiconosVector q1 = bead->q();
    SP::SiconosVector v1 = bead->velocity();
    dataPlot(k, 1) = (*q1)(0);
    dataPlot(k, 2) = (*v1)(0);

    // state q and velocity for the second dynamical system
    SP::LagrangianLinearTIDS bead2 = boost::static_pointer_cast<LagrangianLinearTIDS> (BeadsModel->nonSmoothDynamicalSystem()->dynamicalSystemNumber(2));
    SP::SiconosVector q2 = bead2->q();
    SP::SiconosVector v2 = bead2->velocity();
    dataPlot(k, 3) = (*q2)(0);
    dataPlot(k, 4) = (*v2)(0);

    // state q and velocity for the third dynamical system
    SP::LagrangianLinearTIDS bead3 = boost::static_pointer_cast<LagrangianLinearTIDS> (BeadsModel->nonSmoothDynamicalSystem()->dynamicalSystemNumber(3));
    SP::SiconosVector q3 = bead3->q();
    SP::SiconosVector v3 = bead3->velocity();
    dataPlot(k, 5) = (*q3)(0);
    dataPlot(k, 6) = (*v3)(0);

    cout << " Computation ... " << endl;
    while (s->nextTime() <= BeadsModel->finalT())
    {
      // solve ...
      s->computeOneStep();

      //cout<<"Iteration: "<<k<<endl;
      // --- Get values to be plotted ---

      dataPlot(k, 0) = s->startingTime();
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

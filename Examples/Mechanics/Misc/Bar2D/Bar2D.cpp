/* Siconos-sample version 2.1.1, Copyright INRIA 2005-2007.
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
    Model bar2D("./Bar2D_TIDS.xml");
    cout << "\n *** Bar2D.xml file loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(bar2D.getSimulationPtr());
    cout << "simulation initialization" << endl;
    s->initialize();
    cout << "\n **** the simulation is ready ****" << endl;

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int N = t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 6);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    int k = 0;
    dataPlot(k, 0) = k * t->getH();
    // state q for the first dynamical system (bar2D)
    LagrangianDS* bar = static_cast<LagrangianDS*>(bar2D.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    dataPlot(k, 1) = (bar->getQ())(0);
    // velocity for the bar2D
    dataPlot(k, 2) = (bar->getVelocity())(0);
    // 2nd DS (ground), state q
    LagrangianDS* ground = static_cast<LagrangianDS*>(bar2D.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(1));
    dataPlot(k, 3) = (ground->getQ())(0);
    // velocity
    dataPlot(k, 4) = (ground->getVelocity())(0);
    // reaction (lambda)
    //dataPlot(k, 5) =    (bouncingBar2D.getNonSmoothDynamicalSystem()->getInteractionOnNumber(0)->getLambda())(0);
    dataPlot(k, 5) = 0.0;

    // --- Get interaction between DSs ---
    //cout << "Get interaction ... " << endl;
    //Interaction* I =  bouncingBar2D.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0);
    //I->display();

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
      // solve ...
      s->computeOneStep();
      // --- Get values to be plotted ---
      //time
      dataPlot(k, 0) = k * t->getH();
      // Bar2D: state q
      dataPlot(k, 1) = (bar->getQ())(0);
      // Bar2D: velocity
      dataPlot(k, 2) = (bar->getVelocity())(0);
      // Ground: q
      dataPlot(k, 3) = (ground->getQ())(0);
      // Ground: velocity
      dataPlot(k, 4) = (ground->getVelocity())(0);
      // Reaction
      dataPlot(k, 5) = (bar2D.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
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
    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot);

    // Xml output
    //  Bar2D.saveToXMLFile("./Bar2D_TIDS.xml.output");
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/Bar2D\'" << endl;
  }
}

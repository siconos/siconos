/* Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
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
    Model bouncingBall("./BallTS.xml");
    cout << "\n *** BallTS.xml file loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(bouncingBall.getSimulationPtr());
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
    SimpleMatrix dataPlot(N + 1, 6);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) =  bouncingBall.getT0();
    // state q for the first dynamical system (ball)
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    dataPlot(k, 1) = (ball->getQ())(0);
    dataPlot(k, 2) = (ball->getVelocity())(0);
    dataPlot(k, 3) = (ball->getQ())(1);
    dataPlot(k, 4) = (ball->getVelocity())(1);
    dataPlot(k, 5) = (bouncingBall.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);

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
      dataPlot(k, 0) = bouncingBall.getCurrentT();;
      // Ball: state q
      dataPlot(k, 1) = (ball->getQ())(0);
      // Ball: velocity
      dataPlot(k, 2) = (ball->getVelocity())(0);
      // Ground: q
      dataPlot(k, 3) = (ball->getQ())(1);
      // Ground: velocity
      dataPlot(k, 4) = (ball->getVelocity())(1);
      // Reaction
      dataPlot(k, 5) = (bouncingBall.getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
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
    //  bouncingBall.saveToXMLFile("./BouncingBall_TIDS.xml.output");
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/BouncingBall\'" << endl;
  }
}

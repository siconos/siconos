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


/*!\file BouncingBallTSXml.cpp
\brief \ref EMBouncingBall - C++/XML input file, Time-Stepping version - V. Acary, F. Perignon.

A Ball bouncing on the ground.
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
    Model * bouncingBall = new Model("./BallTS.xml");
    cout << "\n *** BallTS.xml file loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(bouncingBall->getSimulationPtr());
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    cout << "simulation initialization ..." << endl;
    s->initialize();

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();

    int N = t->getNSteps(); // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 4);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time

    int k = 0;
    dataPlot(k, 0) = bouncingBall->getT0();
    // state q for the ball
    dataPlot(k, 1) = (ball->getQ())(0);
    // velocity for the ball
    dataPlot(k, 2) = (ball->getVelocity())(0);
    // Reaction
    dataPlot(k, 3) = (bouncingBall->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);

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
    for (k = 1 ; k < N + 1 ; ++k)
    {
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  bouncingBall->getCurrentT();
      dataPlot(k, 1) = ball->getQ()(0);
      dataPlot(k, 2) = ball->getVelocity()(0);
      dataPlot(k, 3) = (bouncingBall->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);
      s->nextStep();
    }
    cout << "End of computation - Number of iterations done: " << k - 1 << endl;

    // --- elapsed time computing ---
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");


    // Xml output
    //  bouncingBall->saveToXMLFile("./BouncingBall_TIDS.xml.output");
    delete bouncingBall;
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/BouncingBallXml\'" << endl;
  }
}

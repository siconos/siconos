/* Siconos-sample version 1.2.0, Copyright INRIA 2005-2006.
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
#include "Model.h"

#include "LagrangianLinearTIDS.h"

#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "LCP.h"
#include <sys/time.h>
#include <iostream>
#include <math.h>

using namespace std;

bool BallBowl()
{
  bool res = false;

  try
  {

    cout << " **************************************" << endl;
    cout << " ******** Start Ball in a Bowl *********" << endl << endl << endl;
    // --- Model loading from xml file ---
    Model bouncingBall("./BallBowl.xml");
    // --- Get and initialize the simulation ---
    Simulation* s = bouncingBall.getSimulationPtr();
    s->initialize();
    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = t->getK(); // Current step
    int N = t->getNSteps(); // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 3);

    // For the initial time step:
    // time
    dataPlot(k, 0) = k * t->getH();
    // state q for the first dynamical system (ball)
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    dataPlot(k, 1) = (ball->getQ())(0);
    // velocity for the ball
    dataPlot(k, 2) = (ball->getVelocity())(0);
    // --- Time loop  ---
    while (k < N)
    {
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
      // get current time step
      k = t->getK();
      // solve ...
      s->computeOneStep();

      // --- Get values to be plotted ---
      //time
      dataPlot(k, 0) = k * t->getH();
      // Ball: state q
      dataPlot(k, 1) = (ball->getQ())(0);
      // Ball: velocity
      dataPlot(k, 2) = (ball->getVelocity())(0);
    }
    SiconosMatrix * dataRef = new SimpleMatrix("refBallBowl.dat", true);
    double tol = 1e-7;
    double norm = (dataPlot - (*dataRef)).normInf() ;
    cout << endl << endl;
    dataPlot.rawWrite("result.dat", "ascii");
    if (norm < tol)
    {
      cout << " ******** BallBowl global test ended with success ********" << endl;
      res = true;
    }
    else
    {
      cout << " ******** BallBowl global test failed, results differ from those of reference file. ********" << endl;
      res = false;
    }

    cout << endl << endl;

    delete dataRef;
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

  return res;
}

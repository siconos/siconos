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
    TimeStepping* s = static_cast<TimeStepping*>(bouncingBall.getSimulationPtr());
    s->initialize();
    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();
    int k = 0; // Current step
    int N = t->getNSteps(); // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot


    SimpleMatrix dataPlot(N + 1, 6);

    // For the initial time step:
    // time
    dataPlot(k, 0) =  bouncingBall.getT0();
    // state q for the first dynamical system (ball)
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall.getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    SiconosVector * q = ball->getQPtr();
    SiconosVector * v = ball->getVelocityPtr();
    SiconosVector * p = ball->getPPtr(2);

    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*v)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*v)(1);
    dataPlot(k, 5) = (*p)(0);

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
      dataPlot(k, 0) = s->getStartingTime();;
      // Ball: state q
      dataPlot(k, 1) = (*q)(0);
      // Ball: velocity
      dataPlot(k, 2) = (*v)(0);
      // Ground: q
      dataPlot(k, 3) = (*q)(1);
      // Ground: velocity
      dataPlot(k, 4) = (*v)(1);
      // Reaction
      dataPlot(k, 5) = (*p)(0);

      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
    }
    SiconosMatrix * dataRef = new SimpleMatrix("refBallBowl.dat", true);
    double tol = 1e-7;
    double norm = (dataPlot - (*dataRef)).normInf() ;
    cout << endl << endl;
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

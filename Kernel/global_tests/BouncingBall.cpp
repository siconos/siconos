/* Siconos-sample version 2.1.1, Copyright INRIA 2005-2006.
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

// =============================== Ball bouncing on the ground ===============================
//
// Direct description of the model with XML input.
//
// Keywords: LagrangianLinearDS, LagrangianLinear relation, Moreau TimeStepping, LCP, XML.
//
// ======================================================================================================

#include "GlobalTests.h"
#include "SiconosKernel.h"
using namespace std;

bool BouncingBall()
{

  bool res = false;
  try
  {
    cout << " **************************************" << endl;
    cout << " ******** Start Bouncing Ball *********" << endl << endl << endl;
    // --- Model loading from xml file ---
    Model * bouncingBall = new Model("./Ball.xml");
    cout << "\n *** BallTS.xml file loaded ***" << endl;

    // --- Get and initialize the simulation ---
    TimeStepping* s = static_cast<TimeStepping*>(bouncingBall->getSimulationPtr());
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    s->initialize();

    // --- Get the time discretisation scheme ---
    TimeDiscretisation* t = s->getTimeDiscretisationPtr();

    int N = t->getNSteps(); // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 4);
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

    SiconosMatrix * dataRef = new SimpleMatrix("refBouncingBall.dat", true);
    double tol = 1e-9;
    double norm = (dataPlot - (*dataRef)).normInf() ;// diff->normInf();
    if (norm < tol)
    {
      cout << " ******** Bouncing Ball global test ended with success ********" << endl;
      res = true;
    }
    else
    {
      cout << " ******** Bouncing Ball global test failed, results differ from those of reference file. ********" << endl;
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

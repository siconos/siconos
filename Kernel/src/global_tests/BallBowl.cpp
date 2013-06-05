/* Siconos-sample version 3.0.0, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/
#include "SiconosKernel.hpp"


bool BallBowl()
{
  bool res = false;

  try
  {

    std::cout << " **************************************" <<std::endl;
    std::cout << " ******** Start Ball in a Bowl *********" <<std::endl <<std::endl <<std::endl;
    // --- Model loading from xml file ---
    SP::Model bouncingBall(new Model("./BallBowl.xml"));
    bouncingBall->initialize();

    // --- Get and initialize the simulation ---
    SP::TimeStepping s = std11::static_pointer_cast<TimeStepping>(bouncingBall->simulation());
    // --- Get the time discretisation scheme ---
    SP::TimeDiscretisation t = s->timeDiscretisation();
    int k = 0; // Current step
    SP::SiconosMatrix dataRef(new SimpleMatrix("refBallBowl.dat", true));
    int N = dataRef->size(0);
    SimpleMatrix dataPlot(N, 6);

    // For the initial time step:
    // time
    dataPlot(k, 0) =  bouncingBall->t0();
    // state q for the first dynamical system (ball)
    SP::LagrangianDS ball = std11::static_pointer_cast<LagrangianDS> (bouncingBall->nonSmoothDynamicalSystem()->dynamicalSystemNumber(1));
    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(2);

    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*v)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*v)(1);
    dataPlot(k, 5) = (*p)(0);

    std::cout << "Computation ... " <<std::endl;
    // --- Time loop  ---
    while (s->hasNextEvent())
    {
      // solve ...
      s->computeOneStep();

      // --- Get values to be plotted ---
      //time
      dataPlot(k, 0) = s->nextTime();;
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
      k++;
    }
    double tol = 1e-7;
    double norm = (dataPlot - (*dataRef)).normInf() ;
    std::cout <<std::endl <<std::endl;
    if (norm < tol)
    {
      std::cout << " ******** BallBowl global test ended with success ********" <<std::endl;
      res = true;
    }
    else
    {
      std::cout << " ******** BallBowl global test failed, results differ from those of reference file. ********" <<std::endl;
      res = false;
    }

    std::cout <<std::endl <<std::endl;
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    std::cout << e.report() <<std::endl;
  }
  catch (...)
  {
    std::cout << "Exception caught in \'sample/BouncingBall\'" <<std::endl;
  }

  return res;
}

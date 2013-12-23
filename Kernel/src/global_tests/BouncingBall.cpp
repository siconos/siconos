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

// =============================== Ball bouncing on the ground ===============================
//
// Direct description of the model with XML input.
//
// Keywords: LagrangianLinearDS, LagrangianLinear relation, MoreauJeanOSI TimeStepping, LCP, XML.
//
// ======================================================================================================

#include "GlobalTests.hpp"
#include "SiconosKernel.hpp"


bool BouncingBall()
{

  bool res = false;
  try
  {
    std::cout << " **************************************" <<std::endl;
    std::cout << " ******** Start Bouncing Ball *********" <<std::endl <<std::endl <<std::endl;
    // --- Model loading from xml file ---
    SP::Model bouncingBall(new Model("./Ball.xml"));
    // --- Get and initialize the simulation ---
    SP::TimeStepping s = std11::static_pointer_cast<TimeStepping>(bouncingBall->simulation());
    SP::LagrangianDS ball = std11::static_pointer_cast<LagrangianDS> (bouncingBall->nonSmoothDynamicalSystem()->dynamicalSystem(1));
    bouncingBall->initialize();

    int N = 2000; // Number of time steps
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 4);
    // For the initial time step:
    // time

    SP::SiconosVector q = ball->q();
    SP::SiconosVector v = ball->velocity();
    SP::SiconosVector p = ball->p(2);

    dataPlot(0, 0) = bouncingBall->t0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);

    int k = 1;
    while (s->hasNextEvent())
    {
      s->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0);
      s->nextStep();
      k++;
    }

    SP::SiconosMatrix dataRef(new SimpleMatrix("refBouncingBall.dat", true));
    double tol = 1e-9;
    double norm = (dataPlot - (*dataRef)).normInf() ;// diff->normInf();
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");

    if (norm < tol)
    {
      std::cout << " ******** Bouncing Ball global test ended with success ********" <<std::endl;
      res = true;
    }
    else
    {
      std::cout << " ******** Bouncing Ball global test failed, results differ from those of reference file. ********" <<std::endl;
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

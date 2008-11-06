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
    SP::Model bouncingBall(new Model("./BallTS.xml"));
    cout << "\n *** BallTS.xml file loaded ***" << endl;

    cout << "====> Initialisation ..." << endl << endl;
    bouncingBall->initialize();

    // --- Get the simulation ---
    SP::TimeStepping s = boost::static_pointer_cast<TimeStepping>(bouncingBall->getSimulationPtr());
    // --- Get the time discretisation scheme ---
    SP::TimeDiscretisation t = s->getTimeDiscretisationPtr();
    int k = 0;
    double T = bouncingBall->getFinalT();
    double t0 = bouncingBall->getT0();
    double h = s->getTimeStep();
    int N = (int)((T - t0) / h);

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 6);

    cout << "Prepare data for plotting ... " << endl;
    // For the initial time step:
    // time
    dataPlot(k, 0) =  bouncingBall->getT0();
    // state q for the first dynamical system (ball)
    SP::LagrangianDS ball = boost::static_pointer_cast<LagrangianDS> (bouncingBall->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtrNumber(1));
    SP::SiconosVector q = ball->getQPtr();
    SP::SiconosVector v = ball->getVelocityPtr();
    SP::SiconosVector p = ball->getPPtr(2);

    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*v)(0);
    dataPlot(k, 3) = (*q)(1);
    dataPlot(k, 4) = (*v)(1);
    dataPlot(k, 5) = (*p)(0);

    // --- Compute elapsed time ---
    cout << "Computation ... " << endl;
    // --- Time loop  ---
    while (s->getNextTime() <= bouncingBall->getFinalT())
    {
      // solve ...
      s->computeOneStep();

      // --- Get values to be plotted ---
      //time
      dataPlot(k, 0) = s->getNextTime();;
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
      //  dataPlot(k, 6) = osi->computeResidu();
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
      k++;
    }

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix io("result.dat", "ascii");
    //io.write(dataPlot,"noDim");
    io.write(dataPlot);
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

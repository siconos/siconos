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

/*!\file BouncingBallEDXml.cpp
\brief \ref EMBouncingBall - C++/XML input file, Event-Driven version - V. Acary, F. Perignon.

A Ball bouncing on the ground.
Description of the model with XML input.
Simulation with an Event-Driven scheme.
*/

#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {

    // --- Model loading from xml file ---
    cout << "====> Model loading (XML) ..." << endl << endl;
    Model * bouncingBall = new Model("./BallED.xml");
    cout << "\n *** BallED.xml file loaded ***" << endl << endl;

    // --- Get and initialize the simulation ---
    EventDriven* s = static_cast<EventDriven*>(bouncingBall->getSimulationPtr());
    LagrangianDS* ball = static_cast<LagrangianDS*>(bouncingBall->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    cout << "====> Simulation initialisation ..." << endl << endl;
    s->initialize();

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot

    int N = 12368; // Number of saved points: depends on the number of events ...
    unsigned int outputSize = 4;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SiconosVector * q = ball->getQPtr();
    SiconosVector * v = ball->getVelocityPtr();
    SiconosVector * lambda = bouncingBall->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambdaPtr(1);

    dataPlot(0, 0) = bouncingBall->getT0();
    dataPlot(0, 1) = (*q)(0);
    dataPlot(0, 2) = (*v)(0);
    dataPlot(0, 3) = (*lambda)(0);

    cout << "====> Start computation ... " << endl << endl;
    // --- Time loop  ---
    EventsManager * eventsManager = s->getEventsManagerPtr();
    unsigned int numberOfEvent = 0 ;
    int k = 0;
    while (s->hasNextEvent())
    {
      k++;
      s->advanceToEvent();

      s->processEvents();
      // If the treated event is non smooth, we save pre-impact state.
      if (eventsManager->getCurrentEventPtr()->getType() == "NonSmoothEvent")
      {
        dataPlot(k, 0) = s->getCurrentTime();
        dataPlot(k, 1) = (*ball->getQMemoryPtr()->getSiconosVector(1))(0);
        dataPlot(k, 2) = (*ball->getVelocityMemoryPtr()->getSiconosVector(1))(0);
        k++;
      }
      dataPlot(k, 0) = s->getCurrentTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*v)(0);
      dataPlot(k, 3) = (*lambda)(0);
      numberOfEvent++;
    }
    // --- Output files ---
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl << endl;
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    delete ball;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl << endl;
}

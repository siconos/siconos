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
  try
  {

    // --- Model loading from xml file ---
    Model * ball = new Model("./BallED.xml");
    cout << "\n *** BallED.xml file loaded ***" << endl;

    // --- Get and initialize the simulation ---
    Simulation* s;
    s = ball->getSimulationPtr();
    LagrangianLinearTIDS* lds = static_cast<LagrangianLinearTIDS*>(ball->getNonSmoothDynamicalSystemPtr()->getDynamicalSystemPtr(0));
    cout << "simulation initialization ..." << endl;
    s->initialize();

    // --- Get the time discretisation scheme ---
    //TimeDiscretisation* t = s->getTimeDiscretisationPtr();

    EventDriven * eventDriven = static_cast<EventDriven*>(s);

    cout << "End of simulation initialisation" << endl;

    int k = 0; // Current step
    int N = 100000;//t->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 4;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    // time

    dataPlot(k, 0) =  s->getCurrentTime();
    dataPlot(k, 1) = lds->getX()(0);
    dataPlot(k, 2) = lds->getX()(3);
    dataPlot(k, 3) = (ball->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);

    // --- Time loop ---
    cout << " ==== Start of Event Driven simulation - This may take a while ... ====" << endl;
    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    EventsManager * eventsManager = s->getEventsManagerPtr();
    unsigned int numberOfEvent = 0 ;
    while (eventDriven->hasNextEvent())
    {
      k++;
      eventDriven->advanceToEvent();

      eventDriven->processEvents();
      // If the treated event is non smooth, we save pre-impact state.
      if (eventsManager->getCurrentEventPtr()->getType() == "NonSmoothEvent")
      {
        dataPlot(k, 0) = s->getCurrentTime();
        dataPlot(k, 1) = (*lds->getQMemoryPtr()->getSiconosVector(1))(0);
        dataPlot(k, 2) = (*lds->getVelocityMemoryPtr()->getSiconosVector(1))(0);
        k++;
      }

      dataPlot(k, 0) = s->getCurrentTime();
      dataPlot(k, 1) = lds->getX()(0);
      dataPlot(k, 2) = lds->getX()(3);
      dataPlot(k, 3) = (ball->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0)->getLambda(1))(0);

      numberOfEvent++;
    }
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;
    // --- Output files ---
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl;
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
}

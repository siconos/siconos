/* Siconos-sample version 2.1.1, Copyright INRIA 2005-2007.
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

/*!\file BeadsColumnED.cpp
\brief \ref EMBeadsColumnED - C++ input file, Event-Driven version - F. Perignon.

A column of balls
Direct description of the model without XML input.
Simulation with an Event-Driven scheme.
Keywords: LagrangianLinearDS, LagrangianLinear relation, Event Driven, Lsodar, LCP.
*/

#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer t;
  t.restart();
  try
  {

    // User-defined main parameters
    unsigned int dsNumber = 10;      // the number of dynamical systems
    unsigned int nDof = 3;           // degrees of freedom for beads
    double increment_position = 1;   // initial position increment from one DS to the following
    double increment_velocity = 0;   // initial velocity increment from one DS to the following
    double t0 = 0;                   // initial computation time
    double T = 4;                   // final computation time
    double h = 0.005;                 // time step
    double position_init = 5;     // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double R = 0.1;                  // balls radius
    string solverName = "Lemke";     // solver algorithm used for non-smooth problem

    // ================= Creation of the model =======================

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    unsigned int i;
    // A set of DS that will handle all the "balls"
    DynamicalSystemsSet allDS;
    // mass matrix, set to identity
    SiconosMatrix *Mass = new SimpleMatrix(nDof, nDof);
    Mass->eye();
    (*Mass)(2, 2) = 3. / 5 * R * R; // m = 1

    // -- Initial positions and velocities --
    // q0[i] and v0[i] correspond to position and velocity of ball i.
    vector<SimpleVector *> q0;
    vector<SimpleVector *> v0;
    q0.resize(dsNumber, NULL);
    v0.resize(dsNumber, NULL);

    // External forces
    SiconosVector * gravity = new SimpleVector(nDof);
    double m = 1;   // beads mass
    double g = 9.8; // gravity
    (*gravity)(0) = - m * g;

    CheckInsertDS checkDS;
    for (i = 0; i < dsNumber; i++)
    {
      // Memory allocation for q0[i] and v0[i]
      q0[i] = new SimpleVector(nDof);
      v0[i] = new SimpleVector(nDof);
      // set values
      (*(q0[i]))(0) = position_init;
      (*(v0[i]))(0) = velocity_init;
      // Create and insert in allDS a new Lagrangian Linear Dynamical System ...
      checkDS = allDS.insert(new LagrangianLinearTIDS(i, *(q0[i]), *(v0[i]), *Mass));
      // Note that we now use a CheckInsertDS object: checkDS.first is
      // an iterator that points to the DS inserted above.
      //(static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("BeadsPlugin.so", "gravity");
      (static_cast<LagrangianDS*>(*(checkDS.first)))->setFExtPtr(gravity);
      position_init += increment_position;
      velocity_init += increment_velocity;
    }

    // ==> at this point, all the required dynamical systems are saved in allDS.

    // -------------------
    // --- Interactions---
    // -------------------

    InteractionsSet allInteractions;
    // The total number of Interactions
    int interactionNumber = dsNumber;

    // Interaction first bead and floor
    // A set for the systems handles by the "current" Interaction
    DynamicalSystemsSet dsConcerned;
    // Only the "bottom" bead is concerned by this first Interaction,
    // therefore DynamicalSystem number 0.
    dsConcerned.insert(allDS.getPtr(0));

    // -- Newton impact law --
    double e = 0.9;

    // Lagrangian Relation
    unsigned int interactionSize = 1; // y vector size
    SiconosMatrix *H = new SimpleMatrix(interactionSize, nDof);
    (*H)(0, 0) = 1.0;
    NonSmoothLaw * nslaw0 = new NewtonImpactNSL(e);
    SiconosVector *b = new SimpleVector(interactionSize);
    (*b)(0) = -R;
    Relation * relation0 = new LagrangianLinearR(*H, *b);
    unsigned int num = 0 ; // an id number for the Interaction
    Interaction * inter0 = new Interaction("bead-floor", dsConcerned, num, interactionSize, nslaw0, relation0);
    allInteractions.insert(inter0);

    // A list of names for the Interactions
    vector<string> id;
    id.resize(interactionNumber - 1);

    // Interactions ball-ball
    CheckInsertInteraction checkInter;
    // A vector that handles all the relations
    vector<Relation*> LLR(interactionNumber - 1);
    SiconosMatrix *H1 = new SimpleMatrix(1, 2 * nDof);
    if (dsNumber > 1)
    {
      (*H1)(0, 0) = -1.0;
      (*H1)(0, 3) = 1.0;
      for (i = 1; (int)i < interactionNumber; i++)
      {
        // The systems handled by the current Interaction ...
        dsConcerned.clear();
        dsConcerned.insert(allDS.getPtr(i - 1));
        dsConcerned.insert(allDS.getPtr(i));
        ostringstream ostr;
        ostr << i;
        id[i - 1] = ostr.str();
        // The relations
        // Since Ri=Rj and h=0, we do not need to set b.
        LLR[i - 1] = new LagrangianLinearR(*H1);
        checkInter = allInteractions.insert(new Interaction(id[i - 1], dsConcerned, i, interactionSize, nslaw0, LLR[i - 1]));
      }
    }

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    NonSmoothDynamicalSystem * nsds = new NonSmoothDynamicalSystem(allDS, allInteractions);

    // -------------
    // --- Model ---
    // -------------

    Model * multiBeads = new Model(t0, T);
    multiBeads->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    TimeDiscretisation * t = new TimeDiscretisation(h, multiBeads);

    EventDriven* s = new EventDriven(t);

    OneStepIntegrator * OSI = new Lsodar(allDS, s);

    // -- OneStepNsProblem --
    OneStepNSProblem * impact = new LCP(s, "impact", solverName, 10000, 0.0001);
    OneStepNSProblem * acceleration = new LCP(s, "acceleration", solverName, 10000, 0.0001);

    // =========================== End of model definition
    // ================================= Computation
    // --- Simulation initialization ---
    cout << "====> Simulation initialisation ..." << endl << endl;
    s->initialize();

    int k = 0; // Current step
    int N = 14692;// Number of Events

    // Prepare output and save value for the initial time
    unsigned int outputSize = dsNumber * 2 + 1;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // time
    dataPlot(k, 0) =   multiBeads->getT0();
    // Positions and velocities
    i = 0; // Remember that DS are sorted in a growing order according to their number.
    DSIterator it;
    for (it = allDS.begin(); it != allDS.end(); ++it)
    {
      dataPlot(k, (int)i * 2 + 1) = static_cast<LagrangianLinearTIDS*>(*it)->getQ()(0);
      dataPlot(k, (int)i * 2 + 2) = static_cast<LagrangianLinearTIDS*>(*it)->getVelocity()(0);
      i++;
      if ((*it)->getNumber() == 9)
        break;
    }

    // --- Time loop ---
    cout << "Start computation ... " << endl;
    EventsManager * eventsManager = s->getEventsManagerPtr();
    unsigned int numberOfEvent = 0 ;
    while (eventsManager->hasNextEvent())
    {
      k++;
      s->advanceToEvent();

      s->processEvents();

      if (eventsManager->getStartingEventPtr()->getType() == "NonSmoothEvent")
      {
        i = 0; // Remember that DS are sorted in a growing order according to their number.
        dataPlot(k, 0) = s->getStartingTime();
        for (it = allDS.begin(); it != allDS.end(); ++it)
        {
          dataPlot(k, (int)i * 2 + 1) = (*static_cast<LagrangianLinearTIDS*>(*it)->getQMemoryPtr()->getSiconosVector(1))(0);
          dataPlot(k, (int)i * 2 + 2) = (*static_cast<LagrangianLinearTIDS*>(*it)->getVelocityMemoryPtr()->getSiconosVector(1))(0);
          i++;
          if ((*it)->getNumber() == 9)
            break;
        }
        k++;
      }
      // Positions and velocities

      i = 0; // Remember that DS are sorted in a growing order according to their number.
      dataPlot(k, 0) = s->getStartingTime();
      for (it = allDS.begin(); it != allDS.end(); ++it)
      {
        dataPlot(k, (int)i * 2 + 1) = static_cast<LagrangianLinearTIDS*>(*it)->getQ()(0);
        dataPlot(k, (int)i * 2 + 2) = static_cast<LagrangianLinearTIDS*>(*it)->getVelocity()(0);
        i++;
        if ((*it)->getNumber() == 9)
          break;
      }
      numberOfEvent++;
    }
    cout << "===== End of Event Driven simulation. " << numberOfEvent << " events have been processed. ==== " << endl;
    cout << k - numberOfEvent << " impacts detected." << endl;
    // --- Output files ---
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    // --- Free memory ---
    delete impact;
    delete acceleration;
    delete OSI;
    delete t;
    delete s;
    delete multiBeads;
    delete nsds;

    for (i = 0; i < LLR.size(); i++)
      delete LLR[i];
    delete H1;
    delete inter0;
    delete relation0;
    delete nslaw0;
    delete b;
    delete H;
    delete gravity;
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in \'sample/MultiBeadsColumn\'" << endl;
  }
  cout << "Computation Time " << t.elapsed()  << endl;
}

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

/*!\file

\brief \ref EMBeadsColumn - C++ input file, Time-Stepping version - F. Perignon.

A column of balls
Direct description of the model without XML input.
Simulation with a Time-Stepping scheme.
Keywords: LagrangianLinearDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
*/

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  // Exception handling
  try
  {

    boost::timer time;
    time.restart();
    // User-defined main parameters
    unsigned int dsNumber = 10;      // the number of dynamical systems
    unsigned int nDof = 3;           // degrees of freedom for beads
    double increment_position = 1;   // initial position increment from one DS to the following
    double increment_velocity = 0;   // initial velocity increment from one DS to the following
    double t0 = 0;                   // initial computation time
    double T = 4.0;                   // final computation time
    double h = 0.005;                // time step
    double position_init = 5.5;     // initial position for lowest bead.
    double velocity_init = 0.0;      // initial velocity for lowest bead.
    double R = 0.1;                  // balls radius

    // ================= Creation of the model =======================

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." << endl << endl;
    unsigned int i;
    // A set of DS that will handle all the "balls"
    DynamicalSystemsSet allDS;
    // mass matrix, set to identity
    SP::SiconosMatrix Mass(new SimpleMatrix(nDof, nDof));
    Mass->eye();
    (*Mass)(2, 2) = 3. / 5 * R * R; // m = 1

    // -- Initial positions and velocities --
    // q0[i] and v0[i] correspond to position and velocity of ball i.
    vector<SP::SimpleVector> q0;
    vector<SP::SimpleVector> v0;
    q0.resize(dsNumber);
    v0.resize(dsNumber);

    // External forces
    SP::SimpleVector gravity(new SimpleVector(nDof));
    double m = 1;   // beads mass
    double g = 9.8; // gravity
    (*gravity)(0) = - m * g;

    CheckInsertDS checkDS;
    for (i = 0; i < dsNumber; i++)
    {
      // Memory allocation for q0[i] and v0[i]
      q0[i].reset(new SimpleVector(nDof));
      v0[i].reset(new SimpleVector(nDof));
      // set values
      (*(q0[i]))(0) = position_init;
      (*(v0[i]))(0) = velocity_init;
      // Create and insert in allDS a new Lagrangian Linear Dynamical System ...
      SP::LagrangianLinearTIDS ds(new LagrangianLinearTIDS(*(q0[i]), *(v0[i]), *Mass));
      checkDS = allDS.insert(ds);
      // Note that we now use a CheckInsertDS object: checkDS.first is
      // an iterator that points to the DS inserted above.
      //(static_cast<LagrangianDS*>(*(checkDS.first)))->setComputeFExtFunction("BeadsPlugin.so", "gravity");
      (boost::static_pointer_cast<LagrangianLinearTIDS>(*(checkDS.first)))->setFExtPtr(gravity);
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
    SP::SiconosMatrix H(new SimpleMatrix(interactionSize, nDof));
    (*H)(0, 0) = 1.0;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
    SP::SiconosVector b(new SimpleVector(interactionSize));
    (*b)(0) = -R;
    SP::Relation relation0(new LagrangianLinearTIR(*H, *b));
    unsigned int num = 0 ; // an id number for the Interaction
    SP::Interaction inter0(new Interaction("bead-floor", dsConcerned, num, interactionSize, nslaw0, relation0));
    allInteractions.insert(inter0);

    // A list of names for the Interactions
    vector<string> id;
    id.resize(interactionNumber - 1);

    // Interactions ball-ball
    CheckInsertInteraction checkInter;
    // A vector that handles all the relations
    vector<SP::Relation> LLR(interactionNumber - 1);
    SP::SiconosMatrix H1(new SimpleMatrix(1, 2 * nDof));
    (*b)(0) = -2 * R;
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
        LLR[i - 1].reset(new LagrangianLinearTIR(*H1, *b));
        SP::Interaction inter(new Interaction(id[i - 1], dsConcerned, i, interactionSize, nslaw0, LLR[i - 1]));
        checkInter = allInteractions.insert(inter);
      }
    }

    // --------------------------------
    // --- NonSmoothDynamicalSystem ---
    // --------------------------------
    SP::NonSmoothDynamicalSystem nsds(new NonSmoothDynamicalSystem(allDS, allInteractions));

    // -------------
    // --- Model ---
    // -------------

    SP::Model multiBeads(new Model(t0, T));
    multiBeads->setNonSmoothDynamicalSystemPtr(nsds); // set NonSmoothDynamicalSystem of this model

    // ----------------
    // --- Simulation ---
    // ----------------

    // -- Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    SP::TimeStepping s(new TimeStepping(t));

    // -- OneStepIntegrators --
    double theta = 0.5000001;
    SP::OneStepIntegrator OSI(new Moreau(allDS , theta));
    s->insertIntegrator(OSI);
    // That means that all systems in allDS have the same theta value.

    // -- OneStepNsProblem --
    IntParameters iparam(5);
    iparam[0] = 1000; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-5; // Tolerance
    string solverName = "Lemke" ;
    SP::NonSmoothSolver mySolver(new NonSmoothSolver(solverName, iparam, dparam));
    SP::LCP osnspb(new LCP(mySolver));
    s->insertNonSmoothProblem(osnspb);

    //    osnspb->getSolverPtr()->setSolverBlock(true);

    // =========================== End of model definition =================================
    // ================================= Computation =================================
    // --- Simulation initialization ---
    cout << "====> Model initialisation ..." << endl << endl;
    multiBeads->initialize(s);

    int k = 0;
    int N = (int)((T - t0) / h); // Number of time steps

    // Prepare output and save value for the initial time
    unsigned int outputSize = dsNumber * 2 + 1;
    SimpleMatrix dataPlot(N + 1, outputSize); // Output data matrix
    // time
    dataPlot(k, 0) = multiBeads->t0();
    // Positions and velocities
    i = 0; // Remember that DS are sorted in a growing order according to their number.
    DSIterator it;
    for (it = allDS.begin(); it != allDS.end(); ++it)
    {
      dataPlot(k, (int)i * 2 + 1) = (*boost::static_pointer_cast<LagrangianLinearTIDS>(*it)->q())(0);
      dataPlot(k, (int)i * 2 + 2) = (*boost::static_pointer_cast<LagrangianLinearTIDS>(*it)->velocity())(0);
      i++;
      if ((*it)->number() == 9)
        break;
    }
    double cpuTime = 0;
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    boost::timer tt;
    tt.restart();
    while (s->nextTime() <= multiBeads->finalT())
    {
      // solve ...
      try
      {
        s->computeOneStep();
        cpuTime += osnspb-> getCPUtime();
        //      osnspb-> display();
      }
      catch (SiconosException e)
      {
        cout << e.report() << endl;
        osnspb-> display();
        return 0;
        ///        ioMatrix io("M.dat", "ascii");
        //       io.write(osnspb->getM(),"noDim");       cout << k << endl;
        //       ioVector io2("q.dat", "ascii");
        //       io2.write((*osnspb->q()),"noDim");      cout << k << endl;
      }
      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      i = 0;
      for (it = allDS.begin(); it != allDS.end(); ++it)
      {
        dataPlot(k, (int)i * 2 + 1) = (*boost::static_pointer_cast<LagrangianLinearTIDS>(*it)->q())(0);
        dataPlot(k, (int)i * 2 + 2) = (*boost::static_pointer_cast<LagrangianLinearTIDS>(*it)->velocity())(0);
        i++;
        if ((*it)->number() == 9)
          break;
      }
      // transfer of state i+1 into state i and time incrementation
      s->nextStep();
      k++;
    }
    cout << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << tt.elapsed()  << "( " << cpuTime << " )" << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("result.dat", "ascii");
    io.write(dataPlot, "noDim");

    cout << "Total time:" << time.elapsed() << endl;
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

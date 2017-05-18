/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*!\file BouncingBallED.cpp
  \brief \ref EMBouncingBall - C++ input file, Event-Driven version - V. Acary, F. Perignon.

  A Ball bouncing on the ground.
  Direct description of the model.
  Simulation with an Event-Driven scheme.
*/

#include "SiconosKernel.hpp"
const double PI = 3.14159265;
const double g = 9.81; // Gravity

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {
    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDofBall = 1;            // degrees of freedom of ball 1
    double Height = 0.05;         // Distance between impactor balls and monodisperse balls
    // Balls in the tapered chain
    unsigned int NumberBalls  = 10;            // Number
    double R_base_taper = 0.005;         // Base radii of the tapered chain
    double q_taper = 0.0;               // Tapering factor
    // Material properties of balls
    double mass_density = 7780;           // mass density
    double Res_BallBall = 1.0;          // Restitution coefficient at contacts ball-ball
    double Res_BallWall = 1.0;          // Restitution coefficient at contacts ball-wall
    double YoungBall = 203.0e9;           // Young modulus of the balls
    double PoissonBall = 0.3;             // Poison coefficient of the balls
    double YoungWall = 203.0e9;           // Young modulus of the wall
    double PoissonWall = 0.3;             // Poison coefficient of the wall
    double PowCompLaw = 1.5;              // Power of the compliance law: 1.0 for linear contact and 3/2 for the Hertzian contact
    string TypeContactLaw = "BiStiffness"; // Type of compliance contact law
    // Parameters for the global simulation
    double t0 = 0;                         // initial computation time
    double T = 0.5;                        // final computation time
    double h = 0.001;                      // time step
    unsigned int Npointsave = 1000;        // Number of data points to be saved
    // For impact computation
    double DelPest = 1.0e-6;               // Step size estimated for multiple impacts computation
    unsigned int Nstep_save_impact = 100;   // Number of steps every which we save data during impact
    unsigned int Npoint_save_impact = 2000; // Number of points saved during impact
    unsigned int step_start = 0;
    unsigned int step_end = step_start + Npoint_save_impact * Nstep_save_impact;
    unsigned int Nstep_max_impact = 1000000; // Number maximal of steps allowed for impact computation
    string impact_data_name = "data_impact.dat";
    bool _IsSaveDataImpact = false;
    //---------------------------------------
    // ---- Configuration of chaines
    //--------------------------------------
    //************* Balls ******************
    double NumberContacts = NumberBalls ; // Number of contacts
    //(1) Radius of balls
    SP::SiconosVector RadiusBalls(new SiconosVector(NumberBalls));
    for (unsigned int k = 0; k < NumberBalls; ++k)
    {
      (*RadiusBalls)(k) = (pow((1.0 - q_taper), (int)(k + 1))) * R_base_taper;
    }
    // (2) Mass of balls
    SP::SiconosVector MassBalls(new SiconosVector(NumberBalls));
    for (unsigned int id = 0; id < NumberBalls; ++id)
    {
      (*MassBalls)(id) = (4.0 / 3.0) * PI * pow((*RadiusBalls)(id), 3) * mass_density;
    }
    // (3) Initial position of balls
    // For the impactor balls
    SP::SiconosVector InitPosBalls(new SiconosVector(NumberBalls));
    (*InitPosBalls)(0) = Height + (*RadiusBalls)(0);
    for (unsigned int j = 1; j < NumberBalls; ++j)
    {
      (*InitPosBalls)(j) = (*InitPosBalls)(j - 1) + (*RadiusBalls)(j - 1) + (*RadiusBalls)(j);
    }
    // (4) Initial velocity of balls
    SP::SiconosVector InitVelBalls(new SiconosVector(NumberBalls));
    for (unsigned int i = 0; i < NumberBalls; ++i)
    {
      (*InitVelBalls)(i) = 0.0;
    }
    //****************** Contacts ******************
    // (1) Restitution coefficient at contacts
    SP::SiconosVector ResCofContacts(new SiconosVector(NumberContacts));
    SP::SiconosVector ElasCofContacts(new SiconosVector(NumberContacts));
    for (unsigned int id = 0; id < NumberContacts; ++id)
    {
      if (id == 0) // contact ball-wall
      {
        (*ResCofContacts)(id) = Res_BallWall;
      }
      else  // contact ball-ball
      {
        (*ResCofContacts)(id) = Res_BallBall;
      }
      //
      (*ElasCofContacts)(id) = PowCompLaw;
    }
    // (2) Stiffness at contacts
    SP::SiconosVector StiffContacts(new SiconosVector(NumberContacts));
    double Rmoy, Emoy;
    for (unsigned int id = 0; id < NumberContacts; ++id)
    {
      // for ball-wall contact
      if (id == 0)
      {
        Emoy = (4.0 / 3.0) * ((YoungBall * YoungWall) / ((1.0 - pow(PoissonBall, 2)) * YoungWall + (1.0 - pow(PoissonWall, 2)) * YoungBall));
        Rmoy = (*RadiusBalls)(0);
      }
      // Ball-ball contact
      else
      {
        Emoy = (2.0 / 3.0) * (YoungBall / (1.0 - pow(PoissonBall, 2)));
        Rmoy = ((*RadiusBalls)(id - 1) * (*RadiusBalls)(id)) / ((*RadiusBalls)(id - 1) + (*RadiusBalls)(id));
      }
      (*StiffContacts)(id) = pow(Rmoy, 0.5) * Emoy;
    }
    // // Display and save the configuration of the chain simulated
    // cout << "Configuation of ball chains" << endl;
    // cout.precision(15);
    // cout << "Radius of balls: " << endl;
    // RadiusBalls->display();
    // cout << "Mass of balls: " << endl;
    // MassBalls->display();
    // cout << "Initial position of balls: " << endl;
    // InitPosBalls->display();
    // cout << "Initial velocity of balls: " << endl;
    // InitVelBalls->display();
    // cout << "Restitution coefficient at contacts:" << endl;
    // ResCofContacts->display();
    // cout<< "Stiffness at contacts: " << endl;
    // StiffContacts->display();
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;
    std::vector<SP::DynamicalSystem> VecOfallDS;
    SP::SiconosMatrix MassBall;
    SP::SiconosVector q0Ball;
    SP::SiconosVector v0Ball;
    SP::LagrangianLinearTIDS ball;
    SP::SiconosVector FextBall;
    double _Rball, _massBall, _Pos0Ball, _Vel0Ball ;
    // -------------
    // --- Model ---
    // -------------
    SP::Model BallChain(new Model(t0, T));
    // ----------------
    // --- Simulation ---
    // ----------------
    // -- (1) OneStepIntegrators --
    SP::OneStepIntegrator OSI(new LsodarOSI());

    for (unsigned int i = 0; i < NumberBalls; ++i)
    {
      _Rball = (*RadiusBalls)(i); // radius of the ball
      _massBall = (*MassBalls)(i); // mass of the ball
      _Pos0Ball = (*InitPosBalls)(i); // initial position of the ball
      _Vel0Ball = (*InitVelBalls)(i); // initial velocity of the ball
      // Declaration of the DS in Siconos
      MassBall =  SP::SiconosMatrix(new SimpleMatrix(nDofBall, nDofBall));
      (*MassBall)(0, 0) = _massBall;
      // -- Initial positions and velocities --
      q0Ball = SP::SiconosVector(new SiconosVector(nDofBall));
      v0Ball = SP::SiconosVector(new SiconosVector(nDofBall));
      (*q0Ball)(0) = _Pos0Ball;
      (*v0Ball)(0) = _Vel0Ball;
      // -- The dynamical system --
      ball = SP::LagrangianLinearTIDS(new LagrangianLinearTIDS(q0Ball, v0Ball, MassBall));
      // -- Set external forces (weight1) --
      FextBall = SP::SiconosVector(new SiconosVector(nDofBall));
      (*FextBall)(0) = -_massBall * g;
      ball->setFExtPtr(FextBall);
      //
      VecOfallDS.push_back(ball);
      BallChain->nonSmoothDynamicalSystem()->insertDynamicalSystem(ball);
    }
    // --------------------
    // --- Interactions ---
    // --------------------
    SP::SimpleMatrix H;
    SP::SiconosVector E;
    SP::NonSmoothLaw  nslaw;
    SP::Relation relation;
    SP::Interaction interaction;
    double ResCoef, Stiff, ElasPow;
    for (unsigned int j = 0; j < NumberContacts; ++j)
    {
      ResCoef = (*ResCofContacts)(j) ;
      Stiff = (*StiffContacts)(j);
      ElasPow = (*ElasCofContacts)(j);
      if (j == 0) // for contact wall-ball
      {
        H = SP::SimpleMatrix(new SimpleMatrix(1, nDofBall));
        (*H)(0, 0) = 1.0;
        E = SP::SiconosVector(new SiconosVector(1));
        (*E)(0) = -(*RadiusBalls)(j);
      }
      else // For ball-ball contact
      {
        H = SP::SimpleMatrix(new SimpleMatrix(1, (nDofBall + nDofBall)));
        (*H)(0, 0) = -1.0;
        (*H)(0, 1) = 1.0;
        E = SP::SiconosVector(new SiconosVector(1));
        (*E)(0) = -1.0 * ((*RadiusBalls)(j - 1) + (*RadiusBalls)(j));
      }
      //
      nslaw = SP::NonSmoothLaw(new MultipleImpactNSL(ResCoef, Stiff, ElasPow));
      relation = SP::Relation(new LagrangianLinearTIR(H, E));
      interaction = SP::Interaction(new Interaction(nslaw, relation));
      if (j == 0) // for contact wall-ball
        BallChain->nonSmoothDynamicalSystem()->link(interaction, VecOfallDS[j]);
      else // For ball-ball contact
        BallChain->nonSmoothDynamicalSystem()->link(interaction, VecOfallDS[j-1],VecOfallDS[j]);

    }
    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    // -- (3) Non smooth problem --
    SP::OneStepNSProblem impact(new OSNSMultipleImpact(TypeContactLaw, DelPest));
    //SP::OneStepNSProblem impact(new OSNSMultipleImpact(TypeContactLaw,NestImpact));
    SP::OSNSMultipleImpact multiple_impact = std11::dynamic_pointer_cast<OSNSMultipleImpact>(impact);
    multiple_impact->SetSaveData(_IsSaveDataImpact);
    multiple_impact->SetNameOutput(impact_data_name.c_str());
    multiple_impact->SetNstepSave(Nstep_save_impact);
    multiple_impact->SetNstepMax(Nstep_max_impact);
    multiple_impact->SetSizeDataSave(Npoint_save_impact);
    multiple_impact->SetStepMinMaxSave(step_start, step_end);
    SP::OneStepNSProblem acceleration(new LCP());
    // -- (4) Simulation setup with (1) (2) (3)
    SP::EventDriven s(new EventDriven(t));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    s->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);
    BallChain->setSimulation(s);
    
    // =========================== End of model definition ===========================
    //----------------------------------- Initialization-------------------------------
    s->setPrintStat(true);
    BallChain->initialize();
    SP::DynamicalSystemsGraph DSG0 = BallChain->nonSmoothDynamicalSystem()->topology()->dSG(0);
    SP::InteractionsGraph IndexSet0 = BallChain->nonSmoothDynamicalSystem()->topology()->indexSet(0);
    SP::InteractionsGraph IndexSet1 = BallChain->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    SP::InteractionsGraph IndexSet2 = BallChain->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    // // Display topology of the system
    // cout << "Number of vectices of IndexSet0: " << IndexSet0->size() << endl;
    // cout << "Number of vectices of IndexSet1: " << IndexSet1->size() << endl;
    // cout << "Number of vectices of IndexSet2: " << IndexSet2->size() << endl;
    // cout << "Number of vectices of DSG0: " << DSG0->size() << endl;
    //
    SP::EventsManager eventsManager = s->eventsManager();
    boost::progress_display show_progress(Npointsave);
    // ================================= Computation =================================
    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 2 * NumberBalls + 1;
    SimpleMatrix dataPlot(Npointsave, outputSize);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    bool nonSmooth = false;
    unsigned int NumberOfEvents = 0;
    unsigned int NumberOfNSEvents = 0;
    unsigned int k = 0;
    DynamicalSystemsGraph::VIterator ui, uiend;
    //====================================================================
    while ((k < Npointsave) & (s->hasNextEvent()))
    {
      dataPlot(k, 0) =  s->startingTime();
      // Save state of the balls
      unsigned int col_pos = 1;
      unsigned int col_vel = NumberBalls + 1;
      for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui)
      {
        SP::DynamicalSystem ds = DSG0->bundle(*ui);
        SP::LagrangianDS lag_ds = std11::dynamic_pointer_cast<LagrangianDS>(ds);
        SP::SiconosVector q = lag_ds->q();
        SP::SiconosVector v = lag_ds->velocity();
        dataPlot(k, col_pos) = (*q)(0);
        dataPlot(k, col_vel) = (*v)(0);
        col_pos++;
        col_vel++;
      }
      ++k;
      s->advanceToEvent(); // run simulation from one event to the next
      if (eventsManager->nextEvent()->getType() == 2)
      {
        nonSmooth = true;
      };
      //
      s->processEvents();  // process events
      if (nonSmooth)
      {
        // multiple_impact->display();
        dataPlot(k, 0) = s->startingTime();
        // Save state of the balls
        unsigned int col_pos = 1;
        unsigned int col_vel = NumberBalls + 1;
        for (boost::tie(ui, uiend) = DSG0->vertices(); ui != uiend; ++ui)
        {
          SP::DynamicalSystem ds = DSG0->bundle(*ui);
          SP::LagrangianDS lag_ds = std11::dynamic_pointer_cast<LagrangianDS>(ds);
          SP::SiconosVector q = lag_ds->qMemory()->getSiconosVector(1);
          SP::SiconosVector v = lag_ds->velocityMemory()->getSiconosVector(1);
          dataPlot(k, col_pos) = (*q)(0);
          dataPlot(k, col_vel) = (*v)(0);
          col_pos++;
          col_vel++;
        }
        nonSmooth = false;
        ++NumberOfNSEvents;
        ++NumberOfEvents;
        ++show_progress;
        ++k;
      }
      // --- Get values to be plotted ---
      ++NumberOfEvents;
      ++show_progress;
    }

    cout << "Computation Time " << time.elapsed()  << endl;
    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");
    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("resultED_LZBModel.ref", "ascii", dataPlotRef);
    double error = (dataPlot - dataPlotRef).normInf()/ dataPlotRef.normInf();
    std::cout << "Error = "<< error << std::endl;
    if (error > 1e-12)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      return 1;
    }
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught." << endl;
  }
  cout << "Computation Time: " << time.elapsed()  << endl;
}

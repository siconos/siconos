// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include "SiconosKernel.hpp"
#include <stdlib.h>
using namespace std;
#define PI 3.14159
#define GGearth  9.8100
//---------------------------------Decalre global variables ---------------------------------------------------
double LengthBlock = 0.2;        // Length of the rocking block
double HeightBlock = 0.1;        // Height of the rocking block
unsigned int Nfreedom = 3;       // Number of degrees of freedom
unsigned int Ncontact = 2;       // Number of contacts
double MassBlock = 1.0;          // Mass of the rocking block
double PosXiniPointA = 0.0;      // Initial coordinate X of the point A
double PosYiniPointA = 0.0;      // Initial coordinate Y of the point A
double AngleThetaIni = PI / 10.0; // Initial angle theta of the block
double VelXiniPointA = 0.0 ;     // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0 ;     // Initial relative velocity Vy of the point A
double RotVelBlockIni = -0.1;    // Initial angular velocity of the block
double e = 0.9;          // Restitution coefficient
double TimeInitial = 0.0;        // Initial time of the simulation
double TimeFinal =  2.0;     // Final time of the simulation
double StepSize = 0.005;         // Time step size
unsigned int NpointSave = 400;   //
unsigned int SizeOutput = 7;     //
double criterion = 0.05;
unsigned int maxIter = 20000;
//==========================================================================================================
//                                             Main function
//==========================================================================================================
int main(int argc, char* argv[])
{
  //---------------------------- calculate the computation time --------------------------------------------------
  boost::timer time;
  time.restart();
  try
  {
    //===========================================================================================================
    //                  I: Declare the dynamical systems
    //===========================================================================================================
    //1. Set the mass matrix
    SP::SiconosMatrix Mass(new SimpleMatrix(Nfreedom, Nfreedom));
    double InertiaBlock;
    InertiaBlock = (MassBlock / 12.0) * (pow(HeightBlock, 2) + pow(LengthBlock, 2)); // moment of inertia
    (*Mass)(0, 0) = MassBlock;
    (*Mass)(1, 1) = MassBlock;
    (*Mass)(2, 2) = InertiaBlock;
    //2. Set the initial position of the block in function of the initial position of the contact point A (left-hand contact)
    SP::SiconosVector PosIniBlock(new SiconosVector(Nfreedom));
    /*
    (*PosIniBlock)(0) = PosXiniPointA + 0.5*LengthBlock*cos(AngleThetaIni) - 0.5*HeightBlock*sin(AngleThetaIni);
    (*PosIniBlock)(1) = PosYiniPointA + 0.5*LengthBlock*sin(AngleThetaIni) + 0.5*HeightBlock*cos(AngleThetaIni);
    (*PosIniBlock)(2) = AngleThetaIni;
    */
    (*PosIniBlock)(0) = 0.5;
    (*PosIniBlock)(1) = 0.5;
    (*PosIniBlock)(2) = 0.0;

    //3. Set the initial velocity of the block in function of the initial relative velocity of the contact point A
    SP::SiconosVector VelIniBlock(new SiconosVector(Nfreedom));
    /*
    (*VelIniBlock)(0) = VelXiniPointA - (0.5*LengthBlock*sin(AngleThetaIni) + 0.5*HeightBlock*cos(AngleThetaIni))*RotVelBlockIni;
    (*VelIniBlock)(1) = VelYiniPointA + (0.5*LengthBlock*cos(AngleThetaIni) - 0.5*HeightBlock*sin(AngleThetaIni))*RotVelBlockIni;
    (*VelIniBlock)(2) = RotVelBlockIni;
    */
    (*VelIniBlock)(0) = 0.0;
    (*VelIniBlock)(1) = 0.0;
    (*VelIniBlock)(2) = 0.0;

    //4. Instantiate the object of "LagrangianTIDS"
    SP::LagrangianLinearTIDS RockingBlock(new LagrangianLinearTIDS(PosIniBlock, VelIniBlock, Mass));
    //5. Set the external force
    SP::SiconosVector ForceExtern(new SiconosVector(Nfreedom));
    (*ForceExtern)(1) = -MassBlock * GGearth;
    RockingBlock->setFExtPtr(ForceExtern);
    //
    //----------------------------- Display variables of the dynamical system---------------------------------------
    cout << "Initial position of the rocking block:" << endl;
    PosIniBlock->display();
    cout << "Initial velocity of the rocking block:" << endl;
    VelIniBlock->display();
    cout << "Mass matrix of the rocking block:" << endl;
    Mass->display();
    cout << "External force applied on the rocking block:"  << endl;
    ForceExtern->display();
    //==================================================================================================================
    //              II: Declare the relation et interaction between dynamical systems
    //==================================================================================================================
    //
    // Impact law
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    // Interaction at contact point 1
    //SP::Relation relation1(new LagrangianLinearTIR(H, E));
    SP::Relation relation1(new LagrangianScleronomousR("RockingBlockPlugin:h1", "RockingBlockPlugin:G1", "RockingBlockPlugin:G1dot"));
    SP::Interaction inter1(new Interaction(nslaw, relation1));
    // Interaction at contact point 2
    //SP::Relation relation2(new LagrangianLinearTIR(H, E));
    SP::Relation relation2(new LagrangianScleronomousR("RockingBlockPlugin:h2", "RockingBlockPlugin:G2", "RockingBlockPlugin:G2dot"));
    SP::Interaction inter2(new Interaction(nslaw, relation2));
    // Interactions for the whole dynamical system
    //================================================================================================================
    //            III. Create the "model" object
    //================================================================================================================
    SP::Model RoBlockModel(new Model(TimeInitial, TimeFinal));
    RoBlockModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(RockingBlock);
    RoBlockModel->nonSmoothDynamicalSystem()->link(inter1, RockingBlock);
    RoBlockModel->nonSmoothDynamicalSystem()->link(inter2, RockingBlock);
    
    //================================================================================================================
    //            IV. Create the simulation
    //================================================================================================================
    //1. Time discretization
    SP::TimeDiscretisation TimeDiscret(new TimeDiscretisation(TimeInitial, StepSize));
    //2. Integration solver for one step
    SP::OneStepIntegrator OSI(new LsodarOSI());
    //3. Nonsmooth problem
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem acceleration(new LCP());
    //4. Simulation with (1), (2), (3)
    SP::Simulation EDscheme(new EventDriven(TimeDiscret));
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);
    RoBlockModel->setSimulation(EDscheme); // initialize the model

    // bool check1 = EDscheme->hasOneStepNSProblem(impact);
    // bool check2 = EDscheme->hasOneStepNSProblem(acceleration);
    // cout << "Impact law included in the simulation: " << check1 << endl;
    // cout << "LCP at acceleration level included in the simulation: " << check2 << endl;
    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization ------------------------------------------------------
    cout << "====> Simulation initialisation ..." << endl << endl;
    RoBlockModel->initialize(); // initialize the model
    EDscheme->setPrintStat(true);
    SP::EventsManager eventsManager = EDscheme->eventsManager(); // ponters point to the "eventsManager" object
    SP::SiconosVector PosBlock = RockingBlock->q();              // pointer points to the position vector of the rocking block
    SP::SiconosVector VelBlock = RockingBlock->velocity();       // pointer points to the velocity of the rocking block
    //-------------------- Save the output during simulation ---------------------------------------------------------
    SimpleMatrix DataPlot(NpointSave, SizeOutput);
    //------------- At the initial time -----------------------------------------------------------------------------
    DataPlot(0, 0) = RoBlockModel->t0();
    DataPlot(0, 1) = (*PosBlock)(0); // Position X
    DataPlot(0, 2) = (*PosBlock)(1); // Position Y
    DataPlot(0, 3) = (*PosBlock)(2); // Angle theta
    DataPlot(0, 4) = (*VelBlock)(0); // Velocity Vx
    DataPlot(0, 5) = (*VelBlock)(1); // Velocity Vy
    DataPlot(0, 6) = (*VelBlock)(2); // Angular velocity
    //----------------------------------- Simulation starts ----------------------------------------------------------
    cout << "====> Start computation ... " << endl << endl;
    bool NSEvent = false;
    unsigned int NumberNSEvent = 0;
    unsigned int k = 1;
    boost::progress_display show_progress(NpointSave);
    while (EDscheme->hasNextEvent() && (k < NpointSave))
    {
      EDscheme->advanceToEvent(); // lead the simulation run from one event to the next
      //---------- detect the statue of the current event ------------------------------------
      if (eventsManager->nextEvent()->getType() == 2) // the current event is non-smooth
      {
        NSEvent = true;
      };
      EDscheme->processEvents();  // process the current event
      //------------------- get data at the beginning of non-smooth events ---------------------------
      if (NSEvent)
      {
        DataPlot(k, 0) = EDscheme->startingTime(); // instant at non-smooth event
        DataPlot(k, 1) = (*(RockingBlock->qMemory()->getSiconosVector(1)))(0);      // Position X
        DataPlot(k, 2) = (*(RockingBlock->qMemory()->getSiconosVector(1)))(1);      // Position Y
        DataPlot(k, 3) = (*(RockingBlock->qMemory()->getSiconosVector(1)))(2);      // Angle theta
        DataPlot(k, 4) = (*(RockingBlock->velocityMemory()->getSiconosVector(1)))(0); // Velocity Vx
        DataPlot(k, 5) = (*(RockingBlock->velocityMemory()->getSiconosVector(1)))(1); // Velocity Vy
        DataPlot(k, 6) = (*(RockingBlock->velocityMemory()->getSiconosVector(1)))(2); // Angular velocity
        //EDscheme->update(1);
        k++;
        ++NumberNSEvent;
        ++show_progress;
        NSEvent = false;                        // The next event is maybe smooth
      };
      //-------------------- get data at smooth events or at the end of non-smooth events ---------------
      DataPlot(k, 0) = EDscheme->startingTime();
      DataPlot(k, 1) = (*PosBlock)(0); //Position X
      DataPlot(k, 2) = (*PosBlock)(1); //Position Y
      DataPlot(k, 3) = (*PosBlock)(2); // Position theta
      DataPlot(k, 4) = (*VelBlock)(0); // Velocity Vx
      DataPlot(k, 5) = (*VelBlock)(1); // Velocity Vy
      DataPlot(k, 6) = (*VelBlock)(2); // Velocity Vtheta
      // go to the next time step
      k++;
      ++show_progress;
    };
    //----------------------- At the end of the simulation --------------------------
    cout << "End of the simulation" << endl;
    cout << "Number of events processed during simulation: " << (k + 1) << endl;
    cout << "Number of non-smooth events: " << NumberNSEvent << endl;
    cout << "====> Output file writing ..." << endl << endl;
    ioMatrix::write("result.dat", "ascii", DataPlot, "noDim");
  }
  //============================== Catch exceptions ===================================================================
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

// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include "SiconosKernel.hpp"
#include <stdlib.h>
using namespace std;
#define PI 3.141592653589793
#define GGearth  9.8100
//---------------------------------Decalre global variables ---------------------------------------------------
double LengthBlock = 1.0;        // Length of the rocking block
double HeightBlock = 0.5 * LengthBlock; // Height of the rocking block
unsigned int Nfreedom = 3;       // Number of degrees of freedom
unsigned int Ncontact = 2;       // Number of contacts
double MassBlock = 1.0;          // Mass of the rocking block
double PosXiniPointA = 0.0;      // Initial coordinate X of the point A
double PosYiniPointA = 0.0;      // Initial coordinate Y of the point A
double AngleThetaIni = PI / 3.0; // Initial angle theta of the block
double VelXiniPointA = 0.0 ;     // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0 ;     // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;    // Initial angular velocity of the block
// Parameters for impact model
double e1 = 0.8;     // Restitution coefficient at the contact A
double e2 = 0.8;                 // Restitution coefficient at the contact B
double K1 = 1.0e10;              // Stiffness at the contact A
double K2 = 1.0e10;              // Stiffness at the contact B
double eta = 1.0;                // Elasticity coefficient
double DelP = 5.0e-5;            // Step size used for impact resolution
string TypeContactLaw = "BiStiffness"; // Type of compliance contact law
unsigned int Nstep_save_impact = 100;   // Number of steps every which we save data during impact
unsigned int Npoint_save_impact = 1; // Number of points saved during impact
unsigned int step_start = 0;
unsigned int step_end = step_start + Npoint_save_impact * Nstep_save_impact;
unsigned int Nstep_max_impact = 1000000; // Number maximal of steps allowed for impact computation
string impact_data_name = "data_impact.dat";
bool _IsSaveDataImpact = false;
// Parameter for Event-Driven simulation
double TimeInitial = 0.0;        // Initial time of the simulation
double TimeFinal =  2.0;   // Final time of the simulation
double StepSize =   0.01;         // Time step size
unsigned int NpointSave = 252;   //
unsigned int SizeOutput = 13;     //
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
    (*PosIniBlock)(0) = PosXiniPointA + 0.5 * LengthBlock * cos(AngleThetaIni) - 0.5 * HeightBlock * sin(AngleThetaIni);
    (*PosIniBlock)(1) = PosYiniPointA + 0.5 * LengthBlock * sin(AngleThetaIni) + 0.5 * HeightBlock * cos(AngleThetaIni);
    (*PosIniBlock)(2) = AngleThetaIni;


    //3. Set the initial velocity of the block in function of the initial relative velocity of the contact point A
    SP::SiconosVector VelIniBlock(new SiconosVector(Nfreedom));
    (*VelIniBlock)(0) = VelXiniPointA - (0.5 * LengthBlock * sin(AngleThetaIni) + 0.5 * HeightBlock * cos(AngleThetaIni)) * RotVelBlockIni;
    (*VelIniBlock)(1) = VelYiniPointA + (0.5 * LengthBlock * cos(AngleThetaIni) - 0.5 * HeightBlock * sin(AngleThetaIni)) * RotVelBlockIni;
    (*VelIniBlock)(2) = RotVelBlockIni;


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
    /*
    SP::SiconosMatrix H(new SimpleMatrix(1,Nfreedom));
    (*H)(0,1) = 1.0;
    SP::SiconosVector E(new SiconosVector(1));
    (*E)(0) = -0.5*HeightBlock;
    */
    // Impact law
    SP::NonSmoothLaw nslaw1(new MultipleImpactNSL(e1, K1, eta));
    SP::NonSmoothLaw nslaw2(new MultipleImpactNSL(e2, K2, eta));
    // Interaction at contact point 1
    //SP::Relation relation1(new LagrangianLinearTIR(H, E));
    SP::Relation relation1(new LagrangianScleronomousR("RockingBlockPlugin:h1", "RockingBlockPlugin:G1", "RockingBlockPlugin:G1dot"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));
    // Interaction at contact point 2
    //SP::Relation relation2(new LagrangianLinearTIR(H, E));
    SP::Relation relation2(new LagrangianScleronomousR("RockingBlockPlugin:h2", "RockingBlockPlugin:G2", "RockingBlockPlugin:G2dot"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));
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
    SP::OneStepNSProblem impact(new OSNSMultipleImpact(TypeContactLaw, DelP));
    SP::OSNSMultipleImpact multiple_impact = std11::dynamic_pointer_cast<OSNSMultipleImpact>(impact);
    multiple_impact->SetSaveData(_IsSaveDataImpact);
    multiple_impact->SetNameOutput(impact_data_name.c_str());
    multiple_impact->SetNstepSave(Nstep_save_impact);
    multiple_impact->SetNstepMax(Nstep_max_impact);
    multiple_impact->SetSizeDataSave(Npoint_save_impact);
    multiple_impact->SetStepMinMaxSave(step_start, step_end);
    SP::OneStepNSProblem acceleration(new LCP());
    //4. Simulation with (1), (2), (3)
    SP::Simulation EDscheme(new EventDriven(TimeDiscret));
    EDscheme->insertIntegrator(OSI);
    EDscheme->insertNonSmoothProblem(impact, SICONOS_OSNSP_ED_IMPACT);
    EDscheme->insertNonSmoothProblem(acceleration, SICONOS_OSNSP_ED_SMOOTH_ACC);
    // bool check1 = EDscheme->hasOneStepNSProblem(impact);
    // bool check2 = EDscheme->hasOneStepNSProblem(acceleration);
    // cout << "Impact law included in the simulation: " << check1 << endl;
    // cout << "LCP at acceleration level included in the simulation: " << check2 << endl;
    RoBlockModel->setSimulation(EDscheme);
    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization ------------------------------------------------------
    cout << "====> Simulation initialisation ..." << endl << endl;
    EDscheme->setPrintStat(true);
    RoBlockModel->initialize(); // initialize the model


    //SP::LsodarOSI lsodar = std11::static_pointer_cast<LsodarOSI>(OSI);
    //lsodar->setMinMaxStepSizes(1.0e-3,1.0e-3);
    //lsodar->setTol(1,1.0e-3,1.0e-6);
    //lsodar->setMaxOrder(2, 2);


    SP::EventsManager eventsManager = EDscheme->eventsManager(); // ponters point to the "eventsManager" object
    SP::SiconosVector PosBlock = RockingBlock->q();              // pointer points to the position vector of the rocking block
    SP::SiconosVector VelBlock = RockingBlock->velocity();       // pointer points to the velocity of the rocking block
    SP::SiconosVector AcceBlock = RockingBlock->acceleration();       // pointer points to the velocity of the rocking block
    SP::SiconosVector GapCon1 = inter1->y(0);
    SP::SiconosVector GapCon2 = inter2->y(0);
    SP::SiconosVector VelCon1 = inter1->y(1);
    SP::SiconosVector VelCon2 = inter2->y(1);
    SP::SiconosVector LambdaCon1 = inter1->lambda(2);
    SP::SiconosVector LambdaCon2 = inter2->lambda(2);
    SP::InteractionsGraph indexSet0 = RoBlockModel->nonSmoothDynamicalSystem()->topology()->indexSet(0);
    SP::InteractionsGraph indexSet1 = RoBlockModel->nonSmoothDynamicalSystem()->topology()->indexSet(1);
    SP::InteractionsGraph indexSet2 = RoBlockModel->nonSmoothDynamicalSystem()->topology()->indexSet(2);
    cout << "Size of IndexSet0: " << indexSet0->size() << endl;
    cout << "Size of IndexSet1: " << indexSet1->size() << endl;
    cout << "Size of IndexSet2: " << indexSet2->size() << endl;

    InteractionsGraph::VIterator ui, uiend;
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
    DataPlot(0, 7) = (*GapCon1)(0);  // Gap at first contact
    DataPlot(0, 8) = (*GapCon2)(0);  // Gap at second contact
    DataPlot(0, 9) = (*VelCon1)(0);  // Relative velocity at first contact
    DataPlot(0, 10) = (*VelCon2)(0);  // Relative velocity at second contact
    DataPlot(0, 11) = (*LambdaCon1)(0); // Force at first contact
    DataPlot(0, 12) = (*LambdaCon2)(0); // Force at second contact
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
        DataPlot(k, 1) = (*PosBlock)(0); //Position X
        DataPlot(k, 2) = (*PosBlock)(1); //Position Y
        DataPlot(k, 3) = (*PosBlock)(2); // Position theta
        DataPlot(k, 4) = (*(RockingBlock->velocityMemory()->getSiconosVector(1)))(0); // Velocity Vx
        DataPlot(k, 5) = (*(RockingBlock->velocityMemory()->getSiconosVector(1)))(1); // Velocity Vy
        DataPlot(k, 6) = (*(RockingBlock->velocityMemory()->getSiconosVector(1)))(2); // Angular velocity
        DataPlot(k, 7) = (*GapCon1)(0);  // Gap at first contact
        DataPlot(k, 8) = (*GapCon2)(0);  // Gap at second contact
        DataPlot(k, 9) = (*VelCon1)(0);  // Relative velocity at first contact
        DataPlot(k, 10) = (*VelCon2)(0);  // Relative velocity at second contact
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
      DataPlot(k, 7) = (*GapCon1)(0);  // Gap at first contact
      DataPlot(k, 8) = (*GapCon2)(0);  // Gap at second contact
      DataPlot(k, 9) = (*VelCon1)(0);  // Relative velocity at first contact
      DataPlot(k, 10) = (*VelCon2)(0);  // Relative velocity at second contact
      DataPlot(k, 11) = (*LambdaCon1)(0); // Force at first contact
      DataPlot(k, 12) = (*LambdaCon2)(0); // Force at second contact
      // go to the next time step
      k++;
      ++show_progress;
      // // Display information
      // cout << "********At the end of integation step***************"<< (k - 1) << endl;
      // cout << "Information on Dynamical System" << endl;
      // cout << "Position: ";
      // PosBlock->display();
      // cout << "Velocity: ";
      // VelBlock->display();
      // cout << "Acceleration: ";
      // AcceBlock->display();
      // cout << "Information on contacts" << endl;
      // for(std11::tie(ui,uiend) = indexSet0->vertices(); ui!=uiend; ++ui)
      //   {
      //     SP::Interaction inter = indexSet0->bundle(*ui);
      //     cout << "Contact number: " << inter->number() << endl;
      //     cout << "Contact gap: ";
      //     inter->y(0)->display();
      //     cout << "Contact relative velocity: ";
      //     inter->y(1)->display();
      //     cout << "Contact Force: " << endl;
      //     inter->lambda(2)->display();
      //   }
    };


    //----------------------- At the end of the simulation --------------------------
    cout << "End of the simulation" << endl;
    cout << "Number of events processed during simulation: " << (k + 1) << endl;
    cout << "Number of non-smooth events: " << NumberNSEvent << endl;
    cout << "====> Output file writing ..." << endl << endl;
    ioMatrix::write("result.dat", "ascii", DataPlot, "noDim");

    std::cout << "Comparison with a reference file" << std::endl;
    SimpleMatrix dataPlotRef(DataPlot);
    dataPlotRef.zero();
    ioMatrix::read("resultED_LZBmodel_LsodarOSI.ref", "ascii", dataPlotRef);

    double error = (DataPlot - dataPlotRef).normInf()/dataPlotRef.normInf();
    std::cout << "Error = "<< error << std::endl;

    if (error > 1e-10)
    {
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      std::cout << (DataPlot - dataPlotRef).normInf() << std::endl;
      return 1;
    }
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

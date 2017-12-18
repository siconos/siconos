// This is the program to simulate the dynamic of a rocking block by using the Siconos platform
//==================================================================================================================
#include "SiconosKernel.hpp"
#include <stdlib.h>
using namespace std;
#define PI 3.14159
#define GGearth  9.8100
//---------------------------------Decalre global variables ---------------------------------------------------
double LengthBlock = 1.;        // Length of the rocking block
double HeightBlock = 1.5;        // Height of the rocking block
unsigned int Nfreedom = 3;       // Number of degrees of freedom
unsigned int Ncontact = 2;       // Number of contacts
double MassBlock = 1.0;          // Mass of the rocking block
double PosXiniPointA = 0.0;      // Initial coordinate X of the point A
double PosYiniPointA = 0.5;      // Initial coordinate Y of the point A
double AngleThetaIni = PI / 10.0; // Initial angle theta of the block
double VelXiniPointA = 0.0 ;     // Initial relative velocity Vx of the point A
double VelYiniPointA = 0.0 ;     // Initial relative velocity Vy of the point A
double RotVelBlockIni = 0.0;    // Initial angular velocity of the block
double e = 0.5;          // Restitution coefficient
double TimeInitial = 0.0;        // Initial time of the simulation
double TimeFinal =  2.0;     // Final time of the simulation
double StepSize = 0.01;         // Time step size
unsigned int NpointSave = 200;   //
unsigned int SizeOutput = 9;     //
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
    (*PosIniBlock)(0) = PosXiniPointA + 0.5 * LengthBlock * cos(AngleThetaIni) - 0.5 * HeightBlock * sin(AngleThetaIni);
    (*PosIniBlock)(1) = PosYiniPointA + 0.5 * LengthBlock * sin(AngleThetaIni) + 0.5 * HeightBlock * cos(AngleThetaIni);
    (*PosIniBlock)(2) = AngleThetaIni;

    (*PosIniBlock)(0) = 0.0;
    (*PosIniBlock)(1) = 1.;
    (*PosIniBlock)(2) = 0.2;

    //3. Set the initial velocity of the block in function of the initial relative velocity of the contact point A
    SP::SiconosVector VelIniBlock(new SiconosVector(Nfreedom));
    (*VelIniBlock)(0) = VelXiniPointA - (0.5 * LengthBlock * sin(AngleThetaIni) + 0.5 * HeightBlock * cos(AngleThetaIni)) * RotVelBlockIni;
    (*VelIniBlock)(1) = VelYiniPointA + (0.5 * LengthBlock * cos(AngleThetaIni) - 0.5 * HeightBlock * sin(AngleThetaIni)) * RotVelBlockIni;
    (*VelIniBlock)(2) = RotVelBlockIni;
    /*
    (*VelIniBlock)(0) = 0.0;
    (*VelIniBlock)(1) = 0.0;
    (*VelIniBlock)(2) = 0.0;
    */
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
    // Impact law
    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    // Interaction at contact point 1
    SP::Relation relation1(new LagrangianScleronomousR("RockingBlockPlugin:h1", "RockingBlockPlugin:G1"));
    SP::Interaction inter1(new Interaction(nslaw, relation1));
    // Interaction at contact point 2
    SP::Relation relation2(new LagrangianScleronomousR("RockingBlockPlugin:h2", "RockingBlockPlugin:G2"));
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
    SP::OneStepIntegrator OSI(new MoreauJeanCombinedProjectionOSI(0.50001));
    //3. Nonsmooth problem
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem impact_pos(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));
    //4. Simulation with (1), (2), (3)
    SP::TimeStepping TSscheme(new TimeSteppingCombinedProjection(TimeDiscret, OSI, impact, impact_pos));
    RoBlockModel->setSimulation(TSscheme); // initialize the model

    //==================================================================================================================
    //                    V. Process the simulation
    //==================================================================================================================
    // -------------------------------- Simulation initialization ------------------------------------------------------
    cout << "====> Simulation initialisation ..." << endl << endl;
    RoBlockModel->initialize(); // initialize the model
    SP::SiconosVector PosBlock = RockingBlock->q();        // pointer points to the position vector of the rocking block
    SP::SiconosVector VelBlock = RockingBlock->velocity(); // pointer points to the velocity of the rocking block
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

    SP::SiconosVector tmp(new SiconosVector(Nfreedom));
    prod(*Mass, *VelBlock, *tmp, true);
    double kineticEnergy = 0.5 * inner_prod(*VelBlock, *tmp);
    DataPlot(0, 7) = kineticEnergy;

    SP::SiconosVector PosRef(new SiconosVector(Nfreedom));
    (*PosRef)(0) = 0.0;
    (*PosRef)(1) = HeightBlock / 2.0;
    (*PosRef)(2) = 0.0;
    double potentialEnergy = -1.0 * inner_prod(*PosBlock - *PosRef, *ForceExtern);
    DataPlot(0, 8) = potentialEnergy;



    //----------------------------------- Simulation starts ----------------------------------------------------------
    cout << "====> Start computation ... " << endl << endl;
    unsigned int k = 1;
    boost::progress_display show_progress(NpointSave);
    while (k < NpointSave)
    {
      TSscheme->computeOneStep();
      DataPlot(k, 0) = TSscheme->nextTime();
      DataPlot(k, 1) = (*PosBlock)(0); //Position X
      DataPlot(k, 2) = (*PosBlock)(1); //Position Y
      DataPlot(k, 3) = (*PosBlock)(2); // Position theta
      DataPlot(k, 4) = (*VelBlock)(0); // Velocity Vx
      DataPlot(k, 5) = (*VelBlock)(1); // Velocity Vy
      DataPlot(k, 6) = (*VelBlock)(2); // Velocity Vtheta

      prod(*Mass, *VelBlock, *tmp, true);
      kineticEnergy = 0.5 * inner_prod(*VelBlock, *tmp);
      DataPlot(k, 7) = kineticEnergy;

      potentialEnergy = -1.0 * inner_prod(*PosBlock - *PosRef, *ForceExtern);
      DataPlot(k, 8) = potentialEnergy;
      // go to the next time step
      k++;
      ++show_progress;
      TSscheme->nextStep();
    };
    //----------------------- At the end of the simulation --------------------------
    cout << "End of the simulation" << endl;
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

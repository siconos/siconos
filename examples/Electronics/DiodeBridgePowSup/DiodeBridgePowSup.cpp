//-----------------------------------------------------------------------
//
//  DiodeBridgePowSup  : sample of an electrical circuit involving :
//  - a sinusoidal voltage source
//  - a non smooth system : a 4 diodes bridge used as a full wave rectifier
//        of the supplied voltage across the sinusoidal source, providing power
//        to the resistor load of the dynamical system
//  - a linear dynamical system LSDiodeBridgePowSup consisting of
//        a filtering capacitor in parallel with a load resistor
//
//  Expected behavior :
//      The non smooth system is a full wave rectifier :
//      each phase (positive and negative) of the supplied voltage allows current
//      to flow in a constant direction through the load.
//      The capacitor filter acts as a tank providing energy to the load resistor
//      when the voltage across the source weakens.
//      The load resistor consumes energy provided by the source.
//
//  State variable LSDiodeBridgePowSup:
//  - the voltage across the filtering capacitor
//
//  The interaction between the dynamical system and the source is defined by :
//  - complementarity laws between diodes current and voltage. Depending on
//        the diode position in the bridge, y stands for the reverse voltage across
//    the diode or for the diode current (see figure in the template file)
//  - a linear time invariant relation between the state variables and y and
//    lambda (derived from the Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 5e-3;           // Total simulation time
  double h_step = 1.0e-6;    // Time step
  double Rvalue = 1e3;       // load resistance
  double Cfilt  = 300.0e-9;  // filtering capacitor
  double VinitLS = 0.0;     // initial voltage Cfilt
  double DiodeThreshold = 0.21; // Guess what ???
  string Modeltitle = "DiodeBridgePowSup";
  double tinst;
  int k = 0;

  try
  {
    // --- Dynamical system creation ---
    // --- Linear system  (load and filter) specification ---
    SP::SiconosVector init_stateLS(new SiconosVector(1));
    (*init_stateLS)(0) = VinitLS;

    SP::SimpleMatrix LS_A(new SimpleMatrix(1, 1));
    (*LS_A)(0, 0) = -1.0 / (Rvalue * Cfilt);

    SP::FirstOrderLinearDS LSDiodeBridgePowSup(new FirstOrderLinearDS(init_stateLS, LS_A));

    // TODO: review this example with the new way to set the control.

    //    //  Source term "u" specification
    //    LSDiodeBridgePowSup->setUSize(1);
    //    LSDiodeBridgePowSup->setComputeBFunction("./SinPoPlugin","SinPo");

    // --- Interaction between linear system and non smooth system ---

    SP::SimpleMatrix Int_C(new SimpleMatrix(4, 1));
    (*Int_C)(0, 0) =  1.0;
    (*Int_C)(2, 0) =  1.0;

    SP::SimpleMatrix Int_D(new SimpleMatrix(4, 4));

    (*Int_D)(0, 1) = -1.0;
    (*Int_D)(1, 0) =  1.0;
    (*Int_D)(1, 2) =  1.0;
    (*Int_D)(1, 3) = -1.0;
    (*Int_D)(2, 1) = -1.0;
    (*Int_D)(3, 1) =  1.0;

    SP::SiconosVector Offset_y(new SiconosVector(4));
    (*Offset_y)(0) = 1.0;
    (*Offset_y)(2) = 1.0;
    (*Offset_y)(3) = 1.0;
    *Offset_y = -DiodeThreshold * (*Offset_y);

    SP::SiconosVector Offset_lambda(new SiconosVector(4));
    (*Offset_lambda)(1) = 1.0;
    *Offset_lambda = -DiodeThreshold * (*Offset_lambda);

    SP::SiconosVector Int_z(new SiconosVector(5));
    SP::SiconosVector tmp(new SiconosVector(4)) ;
    prod(*Int_D, *Offset_lambda, *tmp);
    *tmp -= *Offset_y;
    Int_z->setBlock(0, *tmp);

    LSDiodeBridgePowSup->setzPtr(Int_z);

    SP::SimpleMatrix Int_B(new SimpleMatrix(1, 4));
    (*Int_B)(0 , 0) = 1.0 / Cfilt;
    (*Int_B)(0 , 2) = 1.0 / Cfilt;

    SP::FirstOrderLinearR LTIRDiodeBridgePowSup(new FirstOrderLinearR(Int_C, Int_B));
    LTIRDiodeBridgePowSup->setDPtr(Int_D);
    LTIRDiodeBridgePowSup->setComputeEFunction("SinPoPlugin", "SinPo");

    SP::ComplementarityConditionNSL nslaw(new ComplementarityConditionNSL(4));

    SP::Interaction InterDiodeBridgePowSup(new Interaction(nslaw, LTIRDiodeBridgePowSup));

    // --- Model creation ---
    SP::Model DiodeBridgePowSup(new Model(t0, T, Modeltitle));
    // add the dynamical system in the non smooth dynamical system
    DiodeBridgePowSup->nonSmoothDynamicalSystem()->insertDynamicalSystem(LSDiodeBridgePowSup);
    // link the interaction and the dynamical system
    DiodeBridgePowSup->nonSmoothDynamicalSystem()->link(InterDiodeBridgePowSup, LSDiodeBridgePowSup);


    // ------------------
    // --- Simulation ---
    // ------------------


    // -- (1) OneStepIntegrators --
    double theta = 0.5;
    SP::EulerMoreauOSI aOSI(new EulerMoreauOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation aTiDisc(new TimeDiscretisation(t0, h_step));

    // -- (3) Non smooth problem
    SP::LCP aLCP(new LCP(SICONOS_LCP_NSQP));

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping aTS(new TimeStepping(aTiDisc, aOSI, aLCP));
    DiodeBridgePowSup->setSimulation(aTS);

    // Initialization
    cout << "====> Initialisation ..." << endl << endl;
    DiodeBridgePowSup->initialize();
    cout << " ---> End of initialization." << endl;



    k = 0;
    double h = aTS->timeStep();
    int N = ceil((T - t0) / h); // Number of time steps
    cout << "Number of time steps = " << N << endl;
    tinst = k * h_step;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int nbPlot = 9;
    SimpleMatrix dataPlot(N + 1, nbPlot);
    //    char buffer[30];
    double i_DF1, i_DR1, i_DF2, i_DR2;
    double v_DF1, v_DR1, v_DF2, v_DR2;

    // For the initial time step:

    i_DF1 = (InterDiodeBridgePowSup->getLambda(0))(2);
    i_DR1 = (InterDiodeBridgePowSup->getLambda(0))(0);
    i_DF2 = (*InterDiodeBridgePowSup->y(0))(1);
    i_DR2 = (InterDiodeBridgePowSup->getLambda(0))(3);

    v_DF1 = -(*InterDiodeBridgePowSup->y(0))(2) + DiodeThreshold;
    v_DR1 = -(*InterDiodeBridgePowSup->y(0))(0) + DiodeThreshold;
    v_DF2 = -(InterDiodeBridgePowSup->getLambda(0))(1) + DiodeThreshold;
    v_DR2 = -(*InterDiodeBridgePowSup->y(0))(3) + DiodeThreshold;

    // time
    dataPlot(k, 0) = DiodeBridgePowSup->t0();

    // source voltage
    dataPlot(k, 1) = (LSDiodeBridgePowSup->z())->getValue(4);

    // source current
    dataPlot(k, 2) = i_DF1 - i_DR2;

    // diode R1 current
    dataPlot(k, 3) = i_DR1;

    // diode R1 voltage
    dataPlot(k, 4) = v_DR1;

    // diode F2 voltage
    dataPlot(k, 5) = v_DF2;

    // diode F1 current
    dataPlot(k, 6) = i_DF1;

    // diode F2 current
    dataPlot(k, 7) = i_DF2;

    // diode R2 current
    dataPlot(k, 8) = i_DR2;


    // --- Compute elapsed time ---
    boost::timer t;
    t.restart();
    // --- Time loop  ---
    while (k < N - 1)
    {
      // get current time step
      k++;
      tinst = k * h_step;

      // solve ...
      aTS->computeOneStep();

      // --- Get values to be plotted ---
      i_DF1 = (InterDiodeBridgePowSup->getLambda(0))(2);
      i_DR1 = (InterDiodeBridgePowSup->getLambda(0))(0);
      i_DF2 = (*InterDiodeBridgePowSup->y(0))(1);
      i_DR2 = (InterDiodeBridgePowSup->getLambda(0))(3);

      v_DF1 = -(*InterDiodeBridgePowSup->y(0))(2) + DiodeThreshold;
      v_DR1 = -(*InterDiodeBridgePowSup->y(0))(0) + DiodeThreshold;
      v_DF2 = -(InterDiodeBridgePowSup->getLambda(0))(1) + DiodeThreshold;
      v_DR2 = -(*InterDiodeBridgePowSup->y(0))(3) + DiodeThreshold;

      // time
      dataPlot(k, 0) = aTS->nextTime();

      // source voltage
      dataPlot(k, 1) = (LSDiodeBridgePowSup->z())->getValue(4);

      // source current
      dataPlot(k, 2) = i_DF1 - i_DR2;

      // diode R1 current
      dataPlot(k, 3) = i_DR1;

      // diode R1 voltage
      dataPlot(k, 4) = v_DR1;

      // diode F2 voltage
      dataPlot(k, 5) = v_DF2;

      // diode F1 current
      dataPlot(k, 6) = i_DF1;

      // diode F2 current
      dataPlot(k, 7) = i_DF2;

      // diode R2 current
      dataPlot(k, 8) = i_DR2;

      aTS->nextStep();
    }


    // --- elapsed time computing ---
    cout << "time = " << t.elapsed() << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;


    // open the file
    char buffer[30];
    ofstream outFile("DiodeBridgePowSup.dat");          // checks that it's opened
    if (!outFile.is_open())
      SiconosMatrixException::selfThrow("function write error : Fail to open \"DiodeBridgePowSup.dat\"");
    for (int i = 0; i < N + 1; i++)
    {
      sprintf(buffer, "%1.10e ", dataPlot(i, 0)); // /!\ depends on machine precision
      outFile << buffer;
      for (int j = 1; j < (int)nbPlot; j++)
      {
        sprintf(buffer, "%1.6e ", dataPlot(i, j)); // /!\ depends on machine precision
        outFile << buffer;
      }
      outFile << endl;
    }
    outFile.close();

  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << "SiconosException at time step " << k << endl;
    cout << e.report() << endl;
  }
  catch (std::exception& e)
  {
    cout << "Exception: " << e.what() << endl;
    exit(-1);
  }
  catch (...)
  {
    cout << "Exception caught " << endl;
  }
}

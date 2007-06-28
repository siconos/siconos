/* Siconos-sample version 2.1.0, Copyright INRIA 2005-2006.
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
//-----------------------------------------------------------------------
//
//  DiodeBridgeCapFilter  : sample of an electrical circuit involving :
//  - a 1st linear dynamical system LSDiodeBridge1 consisting of
//        an LC oscillator (1 µF , 10 mH)
//  - a non smooth system : a 4 diodes bridge used as a full wave rectifier
//        of the supplied voltage across the LC oscillator, providing power
//    to the resistor load of the 2nd dynamical system
//      - a 2nd linear dynamical system LSDiodeBridge2 consisting of
//        a filtering capacitor in parallel with a load resistor
//
//  Expected behavior :
//  The initial state (Vc = 10 V , IL = 0) of the oscillator provides
//      an initial energy.
//  The oscillator period is 2 Pi sqrt(LC) ~ 0,628 ms.
//      The non smooth system is a full wave rectifier :
//  each phase (positive and negative) of the oscillation allows current
//      to flow in a constant direction through the load.
//      The capacitor filter acts as a tank providing energy to the load resistor
//      when the voltage across the oscillator weakens.
//      The load resistor consumes energy : the oscillation damps.
//
//  State variables LSDiodeBridge1:
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  State variable LSDiodeBridge2:
//  - the voltage across the filtering capacitor
//
//  The interaction between the two dynamical systems is defined by :
//  - complementarity laws between diodes current and voltage. Depending on
//        the diode position in the bridge, y stands for the reverse voltage across
//    the diode or for the diode current (see figure in the template file)
//  - a linear time invariant relation between the state variables and y and
//    lambda (derived from the Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include "SiconosKernel.h"

using namespace std;

int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 5e-3;           // Total simulation time
  double h_step = 1e-6;    // Time step
  double Lvalue = 1e-2;      // inductance
  double Cvalue = 1e-6;      // capacitance LC oscillator
  double Rvalue = 1e3;       // load resistance
  double Cfilt  = 300.0e-9;  // filtering capacitor
  double VinitLS1 = 10.0;    // initial voltage LC oscillator
  double VinitLS2 = 0.0;     // initial voltage Cfilt
  string Modeltitle = "DiodeBridgeCapFilter";

  try
  {

    // --- Linear system 1 (LC oscillator) specification ---
    SimpleVector init_stateLS1(2);
    init_stateLS1(0) = VinitLS1;

    SimpleMatrix LS1_A(2, 2);
    LS1_A(0 , 1) = -1.0 / Cvalue;
    LS1_A(1 , 0) = 1.0 / Lvalue;

    cout << " LS1 matrice A = " << endl;
    LS1_A.display();
    FirstOrderLinearDS* LS1DiodeBridgeCapFilter = new FirstOrderLinearDS(1, init_stateLS1, LS1_A);

    // --- Linear system 2 (load and filter) specification ---
    SimpleVector init_stateLS2(1);
    init_stateLS2(0) = VinitLS2;

    SimpleMatrix LS2_A(1, 1);
    LS2_A(0 , 0) = -1.0 / (Rvalue * Cfilt);

    cout << " LS2 matrice A = " << endl;
    LS2_A.display();
    FirstOrderLinearDS* LS2DiodeBridgeCapFilter = new FirstOrderLinearDS(2, init_stateLS2, LS2_A);

    // --- Interaction between linear systems and non smooth system ---

    DynamicalSystemsSet Inter_DS;
    Inter_DS.insert(LS1DiodeBridgeCapFilter);
    Inter_DS.insert(LS2DiodeBridgeCapFilter);

    SiconosMatrix* Int_C = new SimpleMatrix(4, 3);
    (*Int_C)(0 , 2) = 1.0;
    (*Int_C)(2 , 0) = -1.0;
    (*Int_C)(2 , 2) = 1.0;
    (*Int_C)(3 , 0) = 1.0;

    SiconosMatrix* Int_D = new SimpleMatrix(4, 4);
    (*Int_D)(0 , 1) = -1.0;
    (*Int_D)(1 , 0) = 1.0;
    (*Int_D)(1 , 2) = 1.0;
    (*Int_D)(1 , 3) = -1.0;
    (*Int_D)(2 , 1) = -1.0;
    (*Int_D)(3 , 1) = 1.0;

    SiconosMatrix* Int_B = new SimpleMatrix(3, 4);
    (*Int_B)(0 , 2) = -1.0 / Cvalue;
    (*Int_B)(0 , 3) = 1.0 / Cvalue;
    (*Int_B)(2 , 0) = 1.0 / Cfilt;
    (*Int_B)(2 , 2) = 1.0 / Cfilt;

    FirstOrderLinearTIR* LTIRDiodeBridgeCapFilter = new FirstOrderLinearTIR(*Int_C, *Int_B);
    LTIRDiodeBridgeCapFilter->setDPtr(Int_D);
    NonSmoothLaw * nslaw = new ComplementarityConditionNSL(4);

    Interaction* InterDiodeBridgeCapFilter = new Interaction("InterDiodeBridgeCapFilter", Inter_DS, 2, 4, nslaw, LTIRDiodeBridgeCapFilter);

    // --- Model creation ---
    Model DiodeBridgeCapFilter(t0, T, Modeltitle);

    // --- Dynamical system creation ---
    InteractionsSet allInteractions;
    allInteractions.insert(InterDiodeBridgeCapFilter);
    NonSmoothDynamicalSystem* NSDSDiodeBridgeCapFilter = new NonSmoothDynamicalSystem(Inter_DS, allInteractions, false);
    DiodeBridgeCapFilter.setNonSmoothDynamicalSystemPtr(NSDSDiodeBridgeCapFilter);

    // --- Simulation specification---

    TimeDiscretisation* TiDisc = new TimeDiscretisation(h_step, &DiodeBridgeCapFilter);

    TimeStepping* StratDiodeBridgeCapFilter = new TimeStepping(TiDisc);

    double theta = 1.0;
    //Moreau* OSI_LS1DiodeBridgeCapFilter = new Moreau(LS1DiodeBridgeCapFilter,theta,StratDiodeBridgeCapFilter);
    //Moreau* OSI_LS2DiodeBridgeCapFilter = new Moreau(LS2DiodeBridgeCapFilter,theta,StratDiodeBridgeCapFilter);
    Moreau* OSI_LS1DiodeBridgeCapFilter = new Moreau(Inter_DS, theta, StratDiodeBridgeCapFilter);

    LCP* LCP_DiodeBridgeCapFilter = new LCP(StratDiodeBridgeCapFilter, "LCP", "PGS", 10000, 1e-7);

    cout << " -----  Model description ------" << endl;
    DiodeBridgeCapFilter.display();


    StratDiodeBridgeCapFilter->initialize();


    int k = 0;
    int N = TiDisc->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 7);

    // For the initial time step:

    // time
    dataPlot(k, 0) = k * h_step;

    // inductor voltage
    dataPlot(k, 1) = (LS1DiodeBridgeCapFilter->getX())(0);

    // inductor current
    dataPlot(k, 2) = (LS1DiodeBridgeCapFilter->getX())(1);

    // diode R1 current
    dataPlot(k, 3) = (InterDiodeBridgeCapFilter->getLambda(0))(0);

    // diode R1 voltage
    dataPlot(k, 4) = -(InterDiodeBridgeCapFilter->getY(0))(0);

    // diode F2 voltage
    dataPlot(k, 5) = -(InterDiodeBridgeCapFilter->getLambda(0))(1);

    // diode F1 current
    dataPlot(k, 6) = (InterDiodeBridgeCapFilter->getLambda(0))(2);

    // --- Compute elapsed time ---
    double t1, t2, elapsed;
    struct timeval tp;
    int rtn;
    clock_t start, end;
    double elapsed2;
    start = clock();
    rtn = gettimeofday(&tp, NULL);
    t1 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;

    // --- Time loop  ---
    while (k < N)
    {
      // get current time step
      k++;

      // solve ...
      StratDiodeBridgeCapFilter->computeOneStep();

      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = k * h_step;

      // inductor voltage
      dataPlot(k, 1) = (LS1DiodeBridgeCapFilter->getX())(0);

      // inductor current
      dataPlot(k, 2) = (LS1DiodeBridgeCapFilter->getX())(1);

      // diode R1 current
      dataPlot(k, 3) = (InterDiodeBridgeCapFilter->getLambda(0))(0);

      // diode R1 voltage
      dataPlot(k, 4) = -(InterDiodeBridgeCapFilter->getY(0))(0);

      // diode F2 voltage
      dataPlot(k, 5) = -(InterDiodeBridgeCapFilter->getLambda(0))(1);

      // diode F1 current
      dataPlot(k, 6) = (InterDiodeBridgeCapFilter->getLambda(0))(2);

      StratDiodeBridgeCapFilter->nextStep();

    }


    // --- elapsed time computing ---
    end = clock();
    rtn = gettimeofday(&tp, NULL);
    t2 = (double)tp.tv_sec + (1.e-6) * tp.tv_usec;
    elapsed = t2 - t1;
    elapsed2 = (end - start) / (double)CLOCKS_PER_SEC;
    cout << "time = " << elapsed << " --- cpu time " << elapsed2 << endl;

    // Number of time iterations
    cout << "Number of iterations done: " << k << endl;

    // dataPlot (ascii) output
    ioMatrix io("DiodeBridgeCapFilter.dat", "ascii");
    io.write(dataPlot, "noDim");

    delete LCP_DiodeBridgeCapFilter;
    delete OSI_LS1DiodeBridgeCapFilter;
    delete TiDisc;
    delete StratDiodeBridgeCapFilter;
    delete LTIRDiodeBridgeCapFilter;
    delete InterDiodeBridgeCapFilter;
    delete Int_B ;
    delete Int_D ;
    delete Int_C;
    delete LS2DiodeBridgeCapFilter;
    delete LS1DiodeBridgeCapFilter;
    delete NSDSDiodeBridgeCapFilter;

  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << "SiconosException" << endl;
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


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
//-----------------------------------------------------------------------
//
//  CircuitRLCD  : sample of an electrical circuit involving :
//  - a linear dynamical system consisting of an LC oscillator (1 µF , 10 mH)
//  - a non smooth system (a 1000 Ohm resistor in series with a diode) in parallel
//    with the oscillator
//
//  Expected behavior :
//  The initial state of the oscillator provides an initial energy.
//  The period is 2 Pi sqrt(LC) ~ 0,628 ms.
//  A positive voltage across the capacitor allows current to flow
//  through the resistor-diode branch , resulting in an energy loss :
//  the oscillation damps.
//
//  State variables :
//  - the voltage across the capacitor (or inductor)
//  - the current through the inductor
//
//  Since there is only one dynamical system, the interaction is defined by :
//  - a complementarity law between diode current and voltage where y stands
//    for the reverse voltage across the diode and lambda stands for the
//    the diode current
//  - a linear time invariant relation between the state variables and
//    y and lambda (derived from Kirchhoff laws)
//
//-----------------------------------------------------------------------

#include "SiconosKernel.hpp"

using namespace std;

int main(int argc, char* argv[])
{
  double t0 = 0.0;
  double T = 5e-3;        // Total simulation time
  double h_step = 10.0e-7;// Time step
  double Lvalue = 1e-2;   // inductance
  double Cvalue = 1e-6;   // capacitance
  double Rvalue = 1e3;    // resistance
  double Vinit = 10.0;    // initial voltage
  string Modeltitle = "CircuitRLCD";

  try
  {

    // --- Dynamical system specification ---
    SimpleVector init_state(2);
    init_state(0) = Vinit;
    init_state(1) = 0.0;

    SimpleMatrix LS_A(2, 2);
    LS_A(0 , 1) = -1.0 / Cvalue;
    LS_A(1 , 0) = 1.0 / Lvalue;

    SP::FirstOrderLinearTIDS LSCircuitRLCD(new FirstOrderLinearTIDS(init_state, LS_A));

    // --- Interaction between linear system and non smooth system ---

    DynamicalSystemsSet Inter_DS;
    Inter_DS.insert(LSCircuitRLCD);

    SimpleMatrix Int_C(1, 2);
    Int_C(0 , 0) = -1.0;

    SimpleMatrix Int_D(1, 1);
    Int_D(0 , 0) = Rvalue;

    SimpleMatrix Int_B(2, 1);
    Int_B(0 , 0) = -1.0 / Cvalue;

    SP::FirstOrderLinearTIR LTIRCircuitRLCD(new FirstOrderLinearTIR(Int_C, Int_B));
    SP::NonSmoothLaw NSLaw(new ComplementarityConditionNSL(1));

    LTIRCircuitRLCD->setD(Int_D);

    SP::Interaction InterCircuitRLCD(new Interaction("InterCircuitRLCD", Inter_DS, 1, 1, NSLaw, LTIRCircuitRLCD));

    // --- Model creation ---
    SP::Model CircuitRLCD(new Model(t0, T, Modeltitle));

    // --- Non Smooth Dynamical system creation ---

    InteractionsSet allInteractions;
    allInteractions.insert(InterCircuitRLCD);

    SP::NonSmoothDynamicalSystem NSDSCircuitRLCD(new NonSmoothDynamicalSystem(Inter_DS, allInteractions, false));
    CircuitRLCD->setNonSmoothDynamicalSystemPtr(NSDSCircuitRLCD);

    // --- Simulation specification---

    SP::TimeDiscretisation TiDiscRLCD(new TimeDiscretisation(t0, h_step));

    SP::TimeStepping StratCircuitRLCD(new TimeStepping(TiDiscRLCD));

    double theta = 0.5000000000001;

    SP::Moreau OSI_RLCD(new Moreau(LSCircuitRLCD, theta));
    StratCircuitRLCD->recordIntegrator(OSI_RLCD);

    IntParameters iparam(5);
    iparam[0] = 101; // Max number of iteration
    DoubleParameters dparam(5);
    dparam[0] = 1e-8; // Tolerance
    string solverName = "Lemke" ;
    SP::NonSmoothSolver mySolver(new NonSmoothSolver(solverName, iparam, dparam));
    // -- OneStepNsProblem --

    SP::LCP LCP_RLCD(new LCP(mySolver));
    StratCircuitRLCD->recordNonSmoothProblem(LCP_RLCD);

    // Initialization
    CircuitRLCD->initialize(StratCircuitRLCD);
    cout << " -----> End of initialization." << endl;

    double h = StratCircuitRLCD->getTimeStep();
    int N = (int)((T - t0) / h); // Number of time steps
    int k = 0;

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 6);

    // For the initial time step:

    // time
    dataPlot(k, 0) = CircuitRLCD->getT0();

    // inductor voltage
    dataPlot(k, 1) = (LSCircuitRLCD->getX())(0);

    // inductor current
    dataPlot(k, 2) = (LSCircuitRLCD->getX())(1);

    // diode voltage
    dataPlot(k, 3) = - (InterCircuitRLCD->getY(0))(0);

    // diode current
    dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

    dataPlot(k, 5) = (LSCircuitRLCD->getR())(0);

    boost::timer t;
    t.restart();

    // --- Time loop  ---
    for (k = 1 ; k < N ; ++k)
    {
      // solve ...
      StratCircuitRLCD->computeOneStep();

      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = StratCircuitRLCD->getNextTime();

      // inductor voltage
      dataPlot(k, 1) = (LSCircuitRLCD->getX())(0);

      // inductor current
      dataPlot(k, 2) = (LSCircuitRLCD->getX())(1);

      // diode voltage
      dataPlot(k, 3) = - (InterCircuitRLCD->getY(0))(0);

      // diode current
      dataPlot(k, 4) = (InterCircuitRLCD->getLambda(0))(0);

      dataPlot(k, 5) = (LSCircuitRLCD->getR())(0);
      dataPlot(k, 5) = OSI_RLCD->computeResidu();
      // transfer of state i+1 into state i and time incrementation
      StratCircuitRLCD->nextStep();

    }
    // Number of time iterations
    cout << "Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << t.elapsed()  << endl;

    // dataPlot (ascii) output
    ioMatrix io("CircuitRLCD.dat", "ascii");
    io.write(dataPlot, "noDim");
  }

  // --- Exceptions handling ---
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught " << endl;
  }
}

/* Siconos-sample version 2.0.0, Copyright INRIA 2005-2006.
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
#include <SiconosKernel.h>

using namespace std;

int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 2e-3;    // Total simulation time
  double h_step = 5.0e-8;  // Time step
  double Lfvalue = 0.4e-3;   // inductance
  double Cfvalue = 2.2e-6;   // capacitance
  double Lrvalue = 150e-6;   // inductance
  double Crvalue = 68e-9;   // capacitance
  double Rvalue = 33.0; // resistance
  string Modeltitle = "PRC";
  string Author = "Ivan Merillas Santos";
  string Description = " ";
  string Date = "February 2006";
  string Bnamefunction = "Vgen";
  //double Vgen = sin(2*M_PI*55000*0)/fabs(sin(2*M_PI*55000*0));

  // TIMER
  boost::timer time;
  time.restart();

  try
  {

    // --- Dynamical system specification ---

    cout << "====> Model loading ..." << endl << endl;
    SimpleVector init_state(4);
    init_state(0) = 0.0;
    init_state(1) = 0.0;
    init_state(2) = 0.0;
    init_state(3) = 0.0;

    SimpleMatrix LS_A(4, 4);
    LS_A(0 , 1) = -1.0 / Lrvalue;
    LS_A(1 , 0) = 1.0 / Crvalue;
    LS_A(2 , 3) = -1.0 / Lfvalue;
    LS_A(3 , 2) = 1.0 / Cfvalue;
    LS_A(3 , 3) = -1.0 / (Rvalue * Cfvalue);

    FirstOrderLinearDS * LSPRC =  new FirstOrderLinearDS(1, init_state, LS_A);

    // Lrvalue is required in the plug-in, thus we set z[0] = 100.0/ Lrvalue.
    SiconosVector * z = new SimpleVector(1);

    // z[0] is used as a parameter in the plug-in.
    (*z)(0) = 1.0 / Lrvalue;
    LSPRC->setZPtr(z);
    LSPRC->setComputeBFunction("PRCPlugin.so", "computeU");

    // --- Interaction between linear system and non smooth system ---

    DynamicalSystemsSet dsConcerned;
    dsConcerned.insert(LSPRC);

    // -> Relation
    SiconosMatrix* Int_C = new SimpleMatrix(4, 4);
    (*Int_C)(0 , 1) = -1.0;
    (*Int_C)(1 , 1) = 1.0;
    (*Int_C)(2 , 2) = 1.0;
    (*Int_C)(3 , 2) = 1.0;

    SiconosMatrix* Int_D = new SimpleMatrix(4, 4);
    (*Int_D)(0 , 2) = 1.0;
    (*Int_D)(1 , 3) = 1.0;
    (*Int_D)(2 , 0) = -1.0;
    (*Int_D)(3 , 1) = -1.0;

    SiconosMatrix* Int_B = new SimpleMatrix(4, 4);
    (*Int_B)(1 , 0) = -1.0 / Crvalue;
    (*Int_B)(1 , 1) = 1.0 / Crvalue;
    (*Int_B)(2 , 2) = 1.0 / Lfvalue;
    (*Int_B)(2 , 3) = 1.0 / Lfvalue;
    FirstOrderLinearTIR* LTIRPRC = new FirstOrderLinearTIR(Int_C, Int_B);
    LTIRPRC->setDPtr(Int_D);

    // -> Non-smooth law
    NonSmoothLaw * nslaw = new ComplementarityConditionNSL(4);

    Interaction * InterPRC =  new Interaction("InterPRC", dsConcerned, 1, 4, nslaw, LTIRPRC);

    // --- Non Smooth Dynamical system  ---
    NonSmoothDynamicalSystem * NSDSPRC = new NonSmoothDynamicalSystem(LSPRC, InterPRC);

    // --- Model creation ---
    Model * PRC = new Model(t0, T, Modeltitle, Author, Description, Date);
    PRC->setNonSmoothDynamicalSystemPtr(NSDSPRC); // set NonSmoothDynamicalSystem of this model

    // --- Simulation specification---

    // -- Time discretisation --
    TimeDiscretisation* TiDiscRLCD = new TimeDiscretisation(h_step, PRC);
    TimeStepping* StratPRC = new TimeStepping(TiDiscRLCD);

    // -- OneStepIntegrators --
    double theta = 0.5000000000000001;
    OneStepIntegrator* OSI_RLCD = new Moreau(LSPRC, theta, StratPRC);

    // -- OneStepNsProblem --
    string solverName = "Lemke";
    OneStepNSProblem * LCP_RLCD = new LCP(StratPRC, "LCP", solverName, 10000, 0.000001, 1);
    // Note that not all the parameters are usefull, this depends on solverName. But it´s not a
    // problem to give all of them, the useless ones will be ignored.

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    cout << "====> Simulation initialisation ..." << endl << endl;
    StratPRC->initialize();

    int N = TiDiscRLCD->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N, 5);

    SiconosVector * x =  LSPRC->getXPtr();

    // For the initial time step:
    int k = 0;
    // time
    dataPlot(k, 0) = PRC->getCurrentT();

    // inductor voltage
    dataPlot(k, 1) = (*x)(0);

    // inductor current
    dataPlot(k, 2) = (*x)(1);

    // diode R1 current
    dataPlot(k, 3) = (*x)(2);

    // diode R1 voltage
    dataPlot(k, 4) = (*x)(3);

    // --- Time loop  ---
    cout << "====> Start computation ... " << endl << endl;
    for (k = 1 ; k < N ; ++k)
    {
      StratPRC->computeOneStep();
      // --- Get values to be plotted ---
      dataPlot(k, 0) = PRC->getCurrentT();
      dataPlot(k, 1) = (*x)(0);
      dataPlot(k, 2) = (*x)(1);
      dataPlot(k, 3) = (*x)(2);
      dataPlot(k, 4) = (*x)(3);
      // solve ...
      StratPRC->nextStep();
    }
    // Number of time iterations
    cout << "End of computation - Number of iterations done: " << k - 1 << endl;

    // dataPlot (ascii) output
    cout << "====> Output file writing ..." << endl;
    ioMatrix io("PRC.dat", "ascii");
    io.write(dataPlot, "noDim");

    delete LCP_RLCD;
    delete TiDiscRLCD;
    delete StratPRC;
    delete NSDSPRC;
    delete OSI_RLCD;
    delete nslaw;
    delete LTIRPRC;
    delete Int_B ;
    delete Int_D ;
    delete Int_C;
    delete InterPRC;
    delete LSPRC;
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
  cout << "Computation Time " << time.elapsed()  << endl;
}

/* Siconos-sample version 2.0.1, Copyright INRIA 2005-2006.
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

#include "SiconosKernel.h"
using namespace std;

int main(int argc, char* argv[])
{

  double t0 = 0.0;
  double T = 5e-7;        // Total simulation time
  double h_step = 5.0e-8;  // Time step
  double Lfvalue = 0.4e-3;   // inductance
  double Cfvalue = 2.2e-6;   // capacitance
  double Lrvalue = 150e-6;   // inductance
  double Crvalue = 68e-9;   // capacitance
  double Rvalue = 3;    // resistance
  string Modeltitle = "PRC";
  string Author = "Ivan Merillas Santos";
  string Description = " ";
  string Date = "February 2006";
  string Bnamefunction = "Vgen";
  int USize = 1;
  //double Vgen = sin(2*M_PI*55000*0)/fabs(sin(2*M_PI*55000*0));

  try
  {

    // --- Dynamical system specification ---

    SimpleVector init_state(4);
    init_state(0) = 10.0;
    init_state(1) = 10.0;
    init_state(2) = 10.0;
    init_state(3) = 10.0;

    SimpleMatrix LS_A(4, 4);
    LS_A(0 , 0) = 0.0;
    LS_A(0 , 1) = -1.0 / Lrvalue;
    LS_A(0 , 2) = 0.0;
    LS_A(0 , 3) = 0.0;
    LS_A(1 , 0) = 1.0 / Crvalue;
    LS_A(1 , 1) = 0.0;
    LS_A(1 , 2) = 0.0;
    LS_A(1 , 3) = 0.0;
    LS_A(2 , 0) = 0.0;
    LS_A(2 , 1) = 0.0;
    LS_A(2 , 2) = 0.0;
    LS_A(2 , 3) = -1.0 / Lfvalue;
    LS_A(3 , 0) = 0.0;
    LS_A(3 , 1) = 0.0;
    LS_A(3 , 2) = 1.0 / Cfvalue;
    LS_A(3 , 3) = -1.0 / (Rvalue * Cfvalue);

    SiconosMatrix* LS_T = new SimpleMatrix(4, 1);
    LS_T->zero();
    (*LS_T)(0 , 0) = 100.0 / Lrvalue;

    cout << " matrice A = " << endl;
    LS_A.display();
    FirstOrderLinearDS * LSPRC =  new FirstOrderLinearDS(1, init_state, LS_A);

    LSPRC->setUSize(USize);
    LSPRC->setComputeUFunction("PRCPlugin.so", "computeU");
    LSPRC->setTPtr(LS_T);

    //    LSPRC->computeU(0,);

    // --- Interaction between linear system and non smooth system ---

    DynamicalSystemsSet dsConcerned;
    dsConcerned.insert(LSPRC);

    // -> Relation
    SiconosMatrix* Int_C = new SimpleMatrix(4, 4);
    Int_C->zero();
    (*Int_C)(0 , 1) = -1.0;
    (*Int_C)(1 , 1) = 1.0;
    (*Int_C)(2 , 2) = 1.0;
    (*Int_C)(3 , 2) = 1.0;

    SiconosMatrix* Int_D = new SimpleMatrix(4, 4);
    Int_D->zero();
    (*Int_D)(0 , 2) = 1.0;
    (*Int_D)(1 , 3) = 1.0;
    (*Int_D)(2 , 0) = -1.0;
    (*Int_D)(3 , 1) = -1.0;

    SiconosMatrix* Int_B = new SimpleMatrix(4, 4);
    Int_B->zero();
    (*Int_B)(1 , 0) = -1.0 / Crvalue;
    (*Int_B)(1 , 1) = 1.0 / Crvalue;
    (*Int_B)(2 , 2) = 1.0 / Lfvalue;
    (*Int_B)(2 , 3) = 1.0 / Lfvalue;
    LinearTIR* LTIRPRC = new LinearTIR(*Int_C, *Int_B);
    LTIRPRC->setDPtr(Int_D);

    // -> Non-smooth law
    NonSmoothLaw * nslaw = new ComplementarityConditionNSL();

    Interaction * InterPRC =  new Interaction("InterPRC", dsConcerned, 1, 4, nslaw, LTIRPRC);

    // --- Non Smooth Dynamical system  ---
    bool isBVP = 0;
    NonSmoothDynamicalSystem * NSDSPRC = new NonSmoothDynamicalSystem(LSPRC, InterPRC, isBVP);

    // --- Model creation ---
    Model PRC(t0, T, Modeltitle, Author, Description, Date);
    PRC.setNonSmoothDynamicalSystemPtr(NSDSPRC); // set NonSmoothDynamicalSystem of this model

    // --- Simulation specification---

    // -- Time discretisation --
    TimeDiscretisation* TiDiscRLCD = new TimeDiscretisation(h_step, &PRC);
    TimeStepping* StratPRC = new TimeStepping(TiDiscRLCD);

    // -- OneStepIntegrators --
    double theta = 0.5;
    OneStepIntegrator* OSI_RLCD = new Moreau(LSPRC, theta, StratPRC);

    // -- OneStepNsProblem --
    string solverName = "NSQP";
    OneStepNSProblem * LCP_RLCD = new LCP(StratPRC, "LCP", solverName, 101, 0.0001, "max", 0.6);
    // Note that not all the parameters are usefull, this depends on solverName. But it´s not a
    // problem to give all of them, the useless ones will be ignored.

    cout << "=== End of model loading === " << endl;

    cout << " -----  Model description ------" << endl;
    PRC.display();

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    StratPRC->initialize();

    int k = 0;
    Int_C->zero();
    int N = TiDiscRLCD->getNSteps(); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    SimpleMatrix dataPlot(N + 1, 5);

    // For the initial time step:

    // time
    dataPlot(k, 0) = k * h_step;

    // inductor voltage
    dataPlot(k, 1) = (LSPRC->getX())(0);

    // inductor current
    dataPlot(k, 2) = (LSPRC->getX())(1);

    // diode R1 current
    dataPlot(k, 3) = (LSPRC->getX())(2);

    // diode R1 voltage
    dataPlot(k, 4) = (LSPRC->getX())(3);

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
      StratPRC->computeOneStep();

      // --- Get values to be plotted ---
      // time
      dataPlot(k, 0) = k * h_step;

      // inductor voltage
      dataPlot(k, 1) = (LSPRC->getX())(0);

      // inductor current
      dataPlot(k, 2) = (LSPRC->getX())(1);

      // diode R1 current
      dataPlot(k, 3) = (LSPRC->getX())(2);

      // diode R1 voltage
      dataPlot(k, 4) = (LSPRC->getX())(3);

      //cout << (LSPRC->getU())(0) << endl;

      StratPRC->nextStep();

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
    delete LS_T;

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

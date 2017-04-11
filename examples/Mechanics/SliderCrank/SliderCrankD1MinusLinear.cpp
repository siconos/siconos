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


/*!\file
  C++ input file, D1MinusLinearOSI-Time-Stepping version
  T. Schindler, V. Acary

  Slider-crank simulation with a D1MinusLinearOSI-Time-Stepping scheme

  see Flores/Leine/Glocker : Modeling and analysis of planar rigid multibody systems with
  translational clearance joints based on the non-smooth dynamics approach
  */
#include "SiconosKernel.hpp"
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;

int main(int argc, char* argv[])
{
  try
  {
    // ================= Creation of the model =======================

    // parameters according to Table 1
    unsigned int nDof = 3; // degrees of freedom for robot arm
    double t0 = 0.0;         // initial computation time
    double T = 0.2;       // final computation time
    //T=0.00375;
    double h = 1e-5;       // time step : do not decrease, because of strong penetrations

    // geometrical characteristics
    double l1 = 0.1530;
    double l2 = 0.3060;
    double a = 0.05;
    double b = 0.025;
    double c = 0.001;

    // contact parameters
    double e1 = 0.4;
    double e2 = 0.4;
    double e3 = 0.4;
    double e4 = 0.4;
    e1 = 0.1;
    e2 = 0.1;
    e3 = 0.1;
    e4 = 0.1;
    //double mu1 = 0.01;
    //double mu2 = 0.01;
    //double mu3 = 0.01;
    //double mu4 = 0.01;

    // initial conditions
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();
    (*v0)(0) = 150.;
    (*v0)(1) = -75.;

    // t0 = 7e-5;
    // (*q0)(0)=  1.129178e-02;
    // (*q0)(1)= -5.777764e-03;
    // (*q0)(2)=  0.000000e+00;

    // (*v0)(0) = 1.971606e+02 ;
    // (*v0)(1) = -1.064301e+02;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------
    cout << "====> Model loading ..." << endl << endl;

    SP::LagrangianDS slider(new LagrangianDS(q0, v0, "SliderCrankPlugin:mass"));
    slider->setComputeFGyrFunction("SliderCrankPlugin", "FGyr");
    slider->setComputeJacobianFGyrqFunction("SliderCrankPlugin", "jacobianFGyrq");
    slider->setComputeJacobianFGyrqDotFunction("SliderCrankPlugin", "jacobianFGyrqDot");
    slider->setComputeFIntFunction("SliderCrankPlugin", "FInt");
    slider->setComputeJacobianFIntqFunction("SliderCrankPlugin", "jacobianFIntq");
    slider->setComputeJacobianFIntqDotFunction("SliderCrankPlugin", "jacobianFIntqDot");

    // -------------------
    // --- Interactions---
    // -------------------
    // -- corner 1 --
    SP::NonSmoothLaw nslaw1(new NewtonImpactNSL(e1));
    SP::Relation relation1(new LagrangianScleronomousR("SliderCrankPlugin:g1", "SliderCrankPlugin:W1", "SliderCrankPlugin:W1dot"));
    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    // -- corner 2 --
    SP::NonSmoothLaw nslaw2(new NewtonImpactNSL(e2));
    SP::Relation relation2(new LagrangianScleronomousR("SliderCrankPlugin:g2", "SliderCrankPlugin:W2", "SliderCrankPlugin:W2dot"));
    SP::Interaction inter2(new Interaction(nslaw2, relation2));

    // -- corner 3 --
    SP::NonSmoothLaw nslaw3(new NewtonImpactNSL(e3));
    SP::Relation relation3(new LagrangianScleronomousR("SliderCrankPlugin:g3", "SliderCrankPlugin:W3", "SliderCrankPlugin:W3dot"));
    SP::Interaction inter3(new Interaction(nslaw3, relation3));

    // -- corner 4 --
    SP::NonSmoothLaw nslaw4(new NewtonImpactNSL(e4));
    SP::Relation relation4(new LagrangianScleronomousR("SliderCrankPlugin:g4", "SliderCrankPlugin:W4", "SliderCrankPlugin:W4dot"));
    SP::Interaction inter4(new Interaction(nslaw4, relation4));

    // -------------
    // --- Model ---
    // -------------
    SP::Model sliderWithClearance(new Model(t0, T));
    sliderWithClearance->nonSmoothDynamicalSystem()->insertDynamicalSystem(slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter1, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter2, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter3, slider);
    sliderWithClearance->nonSmoothDynamicalSystem()->link(inter4, slider);

    // ----------------
    // --- Simulation ---
    // ----------------
    SP::D1MinusLinearOSI OSI(new D1MinusLinearOSI());
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem force(new LCP());

    SP::TimeSteppingD1Minus s(new TimeSteppingD1Minus(t, 2));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, SICONOS_OSNSP_TS_VELOCITY + 1);
    sliderWithClearance->setSimulation(s);


    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---
    cout << "====> Initialisation ..." << endl << endl;

    sliderWithClearance->initialize();
    int N = ceil((T - t0) / h) + 1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 35;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q = slider->q();
    SP::SiconosVector v = slider->velocity();

    SP::SiconosVector lambda1old = (inter1->lambdaMemory(1))->getSiconosVector(0);
    (*lambda1old)(0);

    int k =0;
    dataPlot(k, 0) = sliderWithClearance->t0();
    dataPlot(k, 1) = (*q)(0) / (2.*M_PI); // crank revolution
    dataPlot(k, 2) = (*q)(1);
    dataPlot(k, 3) = (*q)(2);
    dataPlot(k, 4) = (*v)(0);
    dataPlot(k, 5) = (*v)(1);
    dataPlot(k, 6) = (*v)(2);
    // std::cout << "(*q)(0)= " << (*q)(0)<< std::endl;
    // std::cout << "(*q)(1)= " << (*q)(1)<< std::endl;


    dataPlot(k, 7) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 1 (normalized)
    dataPlot(k, 8) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 2 (normalized)
    dataPlot(k, 9) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (c); // y corner 3 (normalized)
    dataPlot(k, 10) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (c); // y corner 4 (normalized)


    dataPlot(k, 11) = (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1; // x slider (normalized)
    dataPlot(k, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c; // y slider (normalized)

    dataPlot(k, 13) = (*inter1->y(0))(0) ; // g1
    dataPlot(k, 14) = (*inter2->y(0))(0) ; // g2
    dataPlot(k, 15) = (*inter3->y(0))(0) ; // g3
    dataPlot(k, 16) = (*inter4->y(0))(0) ; // g4
    dataPlot(k, 17) = (*inter1->y(1))(0) ; // dot g1
    dataPlot(k, 18) = (*inter2->y(1))(0) ; // dot g2
    dataPlot(k, 19) = (*inter3->y(1))(0) ; // dot g3
    dataPlot(k, 20) = (*inter4->y(1))(0) ; // dot g4
    dataPlot(k, 21) = (*inter1->lambda(1))(0) ; // lambda1
    dataPlot(k, 22) = (*inter2->lambda(1))(0) ; // lambda2
    dataPlot(k, 23) = (*inter3->lambda(1))(0) ; // lambda3
    dataPlot(k, 24) = (*inter4->lambda(1))(0) ; // lambda4
    dataPlot(k, 25) = 0;
    dataPlot(k, 26) = 0;
    dataPlot(k, 27) = (*inter1->lambda(2))(0) ; // lambda1_{k+1}^-
    dataPlot(k, 28) = (*inter2->lambda(2))(0) ; // lambda1_{k+1}^-
    dataPlot(k, 29) = (*inter3->lambda(2))(0) ; // lambda1_{k+1}^-
    dataPlot(k, 30) = (*inter4->lambda(2))(0) ; // lambda1_{k+1}^-



    dataPlot(k, 31) = ( *((inter1->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda1_k^+
    dataPlot(k, 32) = ( *((inter2->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda2_k^+
    dataPlot(k, 33) = ( *((inter3->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda3_k^+
    dataPlot(k, 34) = ( *((inter4->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda4_k^+





    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;

    // ==== Simulation loop - Writing without explicit event handling =====
    k++;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();


//    while ((s->hasNextEvent()) && (k <= 500))
 while ((s->hasNextEvent()))
    {

      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;
      // std::cout <<"Iteration k = " << k <<std::endl;
      // std::cout <<"s->nextTime() = " <<s->nextTime()  <<std::endl;
      // std::cout <<"=====================================================" <<std::endl;

      //std::cout << k << std::endl;
      s->advanceToEvent();

      // --- Get values to be plotted ---
      dataPlot(k, 0) = s->nextTime();
      dataPlot(k, 1) = (*q)(0) / (2.*M_PI); // crank revolution
      dataPlot(k, 2) = (*q)(1);
      dataPlot(k, 3) = (*q)(2);
      dataPlot(k, 4) = (*v)(0);
      dataPlot(k, 5) = (*v)(1);
      dataPlot(k, 6) = (*v)(2);
      dataPlot(k, 7) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 1 (normalized)
      dataPlot(k, 8) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) + b * cos((*q)(2)) - b) / c; // y corner 2 (normalized)
      dataPlot(k, 9) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) - a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (c); // y corner 3 (normalized)
      dataPlot(k, 10) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1)) + a * sin((*q)(2)) - b * cos((*q)(2)) + b) / (c); // y corner 4 (normalized)
      dataPlot(k, 11) = (l1 * cos((*q)(0)) + l2 * cos((*q)(1)) - l2) / l1; // x slider (normalized)
      dataPlot(k, 12) = (l1 * sin((*q)(0)) + l2 * sin((*q)(1))) / c; // y slider (normalized)
      dataPlot(k, 13) = (*inter1->y(0))(0) ; // g1
      dataPlot(k, 14) = (*inter2->y(0))(0) ; // g2
      dataPlot(k, 15) = (*inter3->y(0))(0) ; // g3
      dataPlot(k, 16) = (*inter4->y(0))(0) ; // g4
      dataPlot(k, 17) = (*inter1->y(1))(0) ; // dot g1
      dataPlot(k, 18) = (*inter2->y(1))(0) ; // dot g2
      dataPlot(k, 19) = (*inter3->y(1))(0) ; // dot g3
      dataPlot(k, 20) = (*inter4->y(1))(0) ; // dot g4
      dataPlot(k, 21) = (*inter1->lambda(1))(0) ; // lambda1
      dataPlot(k, 22) = (*inter2->lambda(1))(0) ; // lambda1
      dataPlot(k, 23) = (*inter3->lambda(1))(0) ; // lambda3
      dataPlot(k, 24) = (*inter4->lambda(1))(0) ; // lambda4
      dataPlot(k, 25) = 0;
      dataPlot(k, 26) = 0;
      dataPlot(k, 27) = (*inter1->lambda(2))(0) ; // lambda1_{k+1}^-
      dataPlot(k, 28) = (*inter2->lambda(2))(0) ; // lambda1_{k+1}^-
      dataPlot(k, 29) = (*inter3->lambda(2))(0) ; // lambda1_{k+1}^-
      dataPlot(k, 30) = (*inter4->lambda(2))(0) ; // lambda1_{k+1}^-



      dataPlot(k, 31) = ( *((inter1->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda1_k^+
      dataPlot(k, 32) = ( *((inter2->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda2_k^+
      dataPlot(k, 33) = ( *((inter3->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda3_k^+
      dataPlot(k, 34) = ( *((inter4->lambdaMemory(2))->getSiconosVector(0) )) (0) ; // lambda4_k^+

      // std::cout << "dataPlot(k, 27)" << dataPlot(k, 27)  << std::endl;
      // std::cout << "dataPlot(k, 31)" << dataPlot(k, 31)  << std::endl;



      // std::cout <<" q->display()" <<  std::endl;
      // q->display();
      // std::cout <<" v->display()" <<  std::endl;
      // v->display();


      s->processEvents();
      ++show_progress;
      k++;
    }

    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("result.dat", "ascii", dataPlot, "noDim");


    // Comparison with a reference file
    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("SliderCrankD1MinusLinearOSI.ref", "ascii", dataPlotRef);

    SP::SiconosVector err(new SiconosVector(dataPlot.size(1)));
    (dataPlot - dataPlotRef).normInfByColumn(err);
    err->display();
    double error = 0.0;
    for (unsigned int i = 0; i < err->size(); ++i)
    {
      if (error < (*err)(i)) 
        error = (*err)(i);
    }

    std::cout << "Error = "<< error << std::endl;

    if (error > 1e-11)
    {
    //  (dataPlot - dataPlotRef).display();

      std::cout << "Warning. The result is rather different from the reference file." << std::endl;
      std::cout << "Error = "<< error << std::endl;
      return 1;
    }


  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in SliderCrankD1MinusLinear.cpp" << endl;
  }
}

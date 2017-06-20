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

/*!\file ImpactingBarD1MinusLinear.cpp
  V. Acary

  A Bar bouncing on the ground
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"
#define TS_VELOCITY_LEVEL
using namespace std;

int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

    // User-defined main parameters
#include "UserDefinedParameter.hpp"
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    cout << "====> Model loading ..." <<endl<<endl;

    double l = L/nDof; // length of an element

    cout << "bar length" << l <<  endl;

    SP::SiconosMatrix SparseMass(new SimpleMatrix(nDof,nDof,Siconos::SPARSE,nDof));
    SP::SiconosMatrix SparseStiffness(new SimpleMatrix(nDof,nDof,Siconos::SPARSE,3*nDof));


    SparseMass->setValue(0,0,1.0/3.0);
    SparseMass->setValue(0,1,1.0/6.0);
    SparseStiffness->setValue(0, 0, 1.0);
    SparseStiffness->setValue(0, 1, -1.0);

    for(unsigned int i = 1; i < nDof-1; i++)
    {
      SparseMass->setValue(i,i,2.0/3.0);
      SparseMass->setValue(i,i-1,1.0/6.0);
      SparseMass->setValue(i,i+1,1.0/6.0);


      SparseStiffness->setValue(i,i,2.0);
      SparseStiffness->setValue(i,i-1,-1.0);
      SparseStiffness->setValue(i,i+1,-1.0);
    }

    SparseMass->setValue(nDof-1,nDof-1,1.0/3.0);
    SparseMass->setValue(nDof-1,nDof-2,1.0/6.0);


    SparseStiffness->setValue(nDof-1, nDof-2,-1.0);
    SparseStiffness->setValue(nDof-1,nDof-1,1.0);


    SP::SiconosMatrix SparseDamping(new SimpleMatrix(*SparseStiffness));



    *SparseMass  *= rho*S*l;
    *SparseStiffness  *= E*S/l;
    double xsi = 1000.0;
    std::cout <<  xsi*S/l << std::endl;
    *SparseDamping *= xsi*S/l;


//      SparseMass->display();
//      SparseStiffness->display();


    // -- Initial positions and velocities --
    SP::SiconosVector q0(new SiconosVector(nDof,position_init));
    SP::SiconosVector v0(new SiconosVector(nDof,velocity_init));



    // -- The dynamical system --
    SP::SiconosMatrix SparseMassforDS(new SimpleMatrix(*SparseMass));
    SP::LagrangianLinearTIDS bar(new LagrangianLinearTIDS(q0,v0,SparseMassforDS));

    // -- Set stiffness matrix (weight) --
    bar->setKPtr(SparseStiffness);
    bar->setCPtr(SparseDamping);



    // -- Set external forces (weight) --
    //SP::SiconosVector weight(new SiconosVector(nDof,-g*rho*S/l));
    SP::SiconosVector weight(new SiconosVector(nDof,0.0));
    bar->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction bar-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(1,nDof));
    (*H)(0,0) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::Model impactingBar(new Model(t0, T));

    // add the dynamical system in the non smooth dynamical system
    impactingBar->nonSmoothDynamicalSystem()->insertDynamicalSystem(bar);

    // link the interaction and the dynamical system
    impactingBar->nonSmoothDynamicalSystem()->link(inter,bar);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --

#ifdef TS_VELOCITY_LEVEL
    SP::D1MinusLinearOSI OSI(new D1MinusLinearOSI(D1MinusLinearOSI::halfexplicit_velocity_level));
#else
    SP::D1MinusLinearOSI OSI(new D1MinusLinearOSI());
#endif

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0,h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::OneStepNSProblem impact(new LCP());
    SP::OneStepNSProblem force(new LCP());

    SP::TimeSteppingD1Minus s(new TimeSteppingD1Minus(t, 2));
    s->insertIntegrator(OSI);
    s->insertNonSmoothProblem(impact, SICONOS_OSNSP_TS_VELOCITY);
    s->insertNonSmoothProblem(force, SICONOS_OSNSP_TS_VELOCITY + 1);

    impactingBar->setSimulation(s);

    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout <<"====> Initialisation ..." <<endl<<endl;
    impactingBar->initialize();
    int N = floor((T-t0)/h) +1; // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 13;
    SimpleMatrix dataPlot(N,outputSize);

    SP::SiconosVector q = bar->q();
    SP::SiconosVector v = bar->velocity();
    SP::SiconosVector p = bar->p(1);
    SP::SiconosVector Lambda = inter->lambda(1);
    SP::SiconosVector lambdaminus = inter->lambda(2);
    SP::SiconosVector lambdaplus = ((inter->lambdaMemory(2))->getSiconosVector(0) );
    SP::SiconosVector y = inter->y(0);
    int k = 0;
    dataPlot(k,0) = impactingBar->t0();
    dataPlot(k,1) = (*q)(0);
    dataPlot(k,2) = (*v)(0);
    dataPlot(k,3) = (*p)(0);
    dataPlot(k,4) = (*Lambda)(0);
    dataPlot(k,11) = (*lambdaminus)(0); // lambda1_{k+1}^-
    dataPlot(k,12) = (*lambdaplus)(0);;

    dataPlot(k,7) = (*q)(nDof-1);
    dataPlot(k,8) = (*v)(nDof-1);
    dataPlot(k,9) = (*q)((nDof)/2);
    dataPlot(k,10) = (*v)((nDof)/2);


    SP::SiconosVector tmp(new SiconosVector(nDof));

    prod(*SparseStiffness, *q, *tmp, true);
    double potentialEnergy = inner_prod(*q,   *tmp);
    prod(*SparseMass, *v, *tmp, true);
    double kineticEnergy = inner_prod(*v,*tmp);

    dataPlot(k, 5) = potentialEnergy;
    dataPlot(k, 6) = kineticEnergy;

//    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//     std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



    // --- Time loop ---
    cout << "====> Start computation ... " <<endl<<endl;
    // ==== Simulation loop - Writing without explicit event handling =====

    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

//    while (s->nextTime() < T)
//    while(k < N)
    // while ((s->hasNextEvent()) && (k <= 505))

      while ((s->hasNextEvent()))
    {
      s->advanceToEvent();
//      std::cout << "k = "  << k << std::endl;
//       std::cout << "position"  << std::endl;
//       q->display();
//       std::cout << "velocity"  << std::endl;
//       v->display();

// --- Get values to be plotted ---
      dataPlot(k,0) =s->nextTime();
      dataPlot(k,1) = (*q)(0);
      dataPlot(k,2) = (*v)(0);
      dataPlot(k,3) = (*p)(0);
      dataPlot(k,4) = (*Lambda)(0);

      dataPlot(k,11) = (*lambdaminus)(0); // lambda1_{k+1}^-
      dataPlot(k,12) = (*lambdaplus)(0);;

      dataPlot(k,7) = (*q)(nDof-1);
      dataPlot(k,8) = (*v)(nDof-1);
      dataPlot(k,9) = (*q)((nDof)/2);
      dataPlot(k,10) = (*v)((nDof)/2);

      prod(*SparseStiffness, *q, *tmp, true);
      potentialEnergy = inner_prod(*q,   *tmp);
      prod(*SparseMass, *v, *tmp, true);
      kineticEnergy = inner_prod(*v,*tmp);

      dataPlot(k, 5) = potentialEnergy;
      dataPlot(k, 6) = kineticEnergy;

//      std::cout << "q" << std::endl;
//       q->display();


//       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



      s->processEvents();
      ++show_progress;
      k++;
    }
    cout<<endl << "End of computation - Number of iterations done: "<<k-1<<endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout<<"====> Output file writing ..."<<endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("ImpactingBarD1MinusLinear.dat", "ascii", dataPlot,"noDim");
    // // Comparison with a reference file
    // SimpleMatrix dataPlotRef(dataPlot);
    // dataPlotRef.zero();
    // ioMatrix::read("ImpactingBar.ref", "ascii", dataPlotRef);

    // if ((dataPlot - dataPlotRef).normInf() > 1e-12)
    // {

    //   std::cout << "Warning. The result is rather different from the reference file." << std::endl;
    //   std::cout << "Error = "<< (dataPlot - dataPlotRef).normInf()<<std::endl;
    //   return 1;
    // }


  }

  catch(SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch(...)
  {
    cout << "Exception caught in ImpactingBarTS.cpp" << endl;
  }



}

/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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

#include "SiconosKernel.hpp"
using namespace std;

#include "PunchLagrangianLinearTIDS.hpp"



int main(int argc, char* argv[])
{
  try
  {

    // ================= Creation of the model =======================

#include "UserDefinedParameter.hpp"


    
    // -- The dynamical system --
    const char *mesh_file = "./data/beam-hex.mesh";
    LinearElacticMaterial *  mat = new LinearElacticMaterial(E, nu, rho);
    SP::PunchLagrangianLinearTIDS  punch (new PunchLagrangianLinearTIDS(mesh_file, *mat));

    unsigned int nDof = punch->dimension();

    // -- Set external forces (weight) --
    //SP::SiconosVector weight(new SiconosVector(nDof,-g*rho*S/l));
    SP::SiconosVector weight(new SiconosVector(nDof,0.0));
    punch->setFExtPtr(weight);

    // --------------------
    // --- Interactions ---
    // --------------------

    // -- nslaw --
    double e = 0.0;

    // Interaction punch-floor
    //
    SP::SimpleMatrix H(new SimpleMatrix(1,nDof));
    (*H)(0,0) = 1.0;

    SP::NonSmoothLaw nslaw(new NewtonImpactNSL(e));
    SP::Relation relation(new LagrangianLinearTIR(H));

    SP::Interaction inter(new Interaction(nslaw, relation));

    // -------------
    // --- Model ---
    // -------------
    SP::NonSmoothDynamicalSystem impactingPunch(new NonSmoothDynamicalSystem(t0, T));

    // add the dynamical system in the non smooth dynamical system
    impactingPunch->insertDynamicalSystem(punch);

    // link the interaction and the dynamical system
    impactingPunch->link(inter,punch);


    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --

    SP::MoreauJeanOSI OSI(new MoreauJeanOSI(theta,0.5));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0,h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new LCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(impactingPunch, t,OSI,osnspb));


    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---


    int N = floor((T-t0)/h); // Number of time steps

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 11;
    SimpleMatrix dataPlot(N,outputSize);

    SP::SiconosVector q = punch->q();
    SP::SiconosVector v = punch->velocity();
    SP::SiconosVector p = punch->p(1);
    SP::SiconosVector lambda = inter->lambda(1);

    dataPlot(0, 0) = impactingPunch->t0();
    dataPlot(0,1) = (*q)(0);
    dataPlot(0,2) = (*v)(0);
    dataPlot(0, 3) = (*p)(0);
    dataPlot(0, 4) = (*lambda)(0);
    dataPlot(0,7) = (*q)(nDof-1);
    dataPlot(0,8) = (*v)(nDof-1);
    dataPlot(0,9) = (*q)((nDof)/2);
    dataPlot(0,10) = (*v)((nDof)/2);


    SP::SiconosVector tmp(new SiconosVector(nDof));

    // prod(*SparseStiffness, *q, *tmp, true);
    // double potentialEnergy = inner_prod(*q,   *tmp);
    // prod(*SparseMass, *v, *tmp, true);
    // double kineticEnergy = inner_prod(*v,*tmp);

    // dataPlot(0, 5) = potentialEnergy;
    // dataPlot(0, 6) = kineticEnergy;

//    std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//     std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



    // --- Time loop ---
    cout << "====> Start computation ... " <<endl<<endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 1;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();

//    while (s->nextTime() < T)
    while(k < N)
    {
      s->computeOneStep();

//       std::cout << "position"  << std::endl;
//       q->display();
//       std::cout << "velocity"  << std::endl;
//       v->display();

// --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k,1) = (*q)(0);
      dataPlot(k,2) = (*v)(0);
      dataPlot(k, 3) = (*p)(0)/h;
      dataPlot(k, 4) = (*lambda)(0);
      dataPlot(k,7) = (*q)(nDof-1);
      dataPlot(k,8) = (*v)(nDof-1);
      dataPlot(k,9) = (*q)((nDof)/2);
      dataPlot(k,10) = (*v)((nDof)/2);



      // prod(*SparseStiffness, *q, *tmp, true);
      // potentialEnergy = inner_prod(*q,   *tmp);
      // prod(*SparseMass, *v, *tmp, true);
      // kineticEnergy = inner_prod(*v,*tmp);

      // dataPlot(k, 5) = potentialEnergy;
      // dataPlot(k, 6) = kineticEnergy;

//      std::cout << "q" << std::endl;
//       q->display();


//       std::cout <<"potentialEnergy ="<<potentialEnergy << std::endl;
//       std::cout <<"kineticEnergy ="<<kineticEnergy << std::endl;



      s->nextStep();
      ++show_progress;
      k++;
    }
    cout<<endl << "End of computation - Number of iterations done: "<<k-1<<endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout<<"====> Output file writing ..."<<endl;
    ioMatrix::write("ImpactingPunch.dat", "ascii", dataPlot,"noDim");
//     cout << " Comparison with a reference file" << endl;
//     SimpleMatrix dataPlotRef(dataPlot);
//     dataPlotRef.zero();
//     ioMatrix::read("ImpactingPunch.ref", "ascii", dataPlotRef);

//     double error = (dataPlot - dataPlotRef).normInf() ;
// cout << "Error = " << error << endl;
//     if (error > 1e-11)
//     {
//       std::cout << "Warning. The result is rather different from the reference file." << std::endl;
//       std::cout << "Error = "<< (dataPlot - dataPlotRef).normInf()<<std::endl;
//       return 1;
//     }


  }

  catch(SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch(...)
  {
    cout << "Exception caught in ImpactingPunchTS.cpp" << endl;
  }



}

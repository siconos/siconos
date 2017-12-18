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

/*!\file NE....cpp
  \brief \ref EMNE_MULTIBODY - C++ input file, Time-Stepping version - O.B.

  A multibody example.
  Direct description of the model.
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"
#include "KneeJointR.hpp"
#include "PrismaticJointR.hpp"
#include <boost/math/quaternion.hpp>
using namespace std;

/* Given a position of a point in the Inertial Frame and the configuration vector q of a solid
 * returns a position in the spatial frame.
 */
void fromInertialToSpatialFrame(double *positionInInertialFrame, double *positionInSpatialFrame, SP::SiconosVector  q  )
{
double q0 = q->getValue(3);
double q1 = q->getValue(4);
double q2 = q->getValue(5);
double q3 = q->getValue(6);

::boost::math::quaternion<double>    quatQ(q0, q1, q2, q3);
::boost::math::quaternion<double>    quatcQ(q0, -q1, -q2, -q3);
::boost::math::quaternion<double>    quatpos(0, positionInInertialFrame[0], positionInInertialFrame[1], positionInInertialFrame[2]);
::boost::math::quaternion<double>    quatBuff;

//perform the rotation
quatBuff = quatQ * quatpos * quatcQ;

positionInSpatialFrame[0] = quatBuff.R_component_2()+q->getValue(0);
positionInSpatialFrame[1] = quatBuff.R_component_3()+q->getValue(1);
positionInSpatialFrame[2] = quatBuff.R_component_4()+q->getValue(2);

}
void tipTrajectories(SP::SiconosVector  q, double * traj, double length)
{
  double positionInInertialFrame[3];
  double positionInSpatialFrame[3];
  // Output the position of the tip of beam1
  positionInInertialFrame[0]=length/2;
  positionInInertialFrame[1]=0.0;
  positionInInertialFrame[2]=0.0;

  fromInertialToSpatialFrame(positionInInertialFrame, positionInSpatialFrame, q  );
  traj[0] = positionInSpatialFrame[0];
  traj[1] = positionInSpatialFrame[1];
  traj[2] = positionInSpatialFrame[2];


  // std::cout <<  "positionInSpatialFrame[0]" <<  positionInSpatialFrame[0]<<std::endl;
  // std::cout <<  "positionInSpatialFrame[1]" <<  positionInSpatialFrame[1]<<std::endl;
  // std::cout <<  "positionInSpatialFrame[2]" <<  positionInSpatialFrame[2]<<std::endl;

  positionInInertialFrame[0]=-length/2;
  fromInertialToSpatialFrame(positionInInertialFrame, positionInSpatialFrame, q  );
  traj[3]= positionInSpatialFrame[0];
  traj[4] = positionInSpatialFrame[1];
  traj[5] = positionInSpatialFrame[2];
}





int main(int argc, char* argv[])
{
  try
  {


    // ================= Creation of the model =======================

    // User-defined main parameters
    unsigned int nDof = 3;
    unsigned int qDim = 7;
    unsigned int nDim = 6;
    double t0 = 0;                   // initial computation time
    double T = 10.0;                  // final computation time
    double h = 0.01;                // time step
    int N = 1000;
    double L1 = 1.0;
    double L2 = 1.0;
    double L3 = 1.0;
    double theta = 1.0;              // theta for MoreauJeanOSI integrator
    double g = 9.81; // Gravity
    double m = 1.;

    // -------------------------
    // --- Dynamical systems ---
    // -------------------------

    FILE * pFile;
    pFile = fopen("data.h", "w");
    if (pFile == NULL)
    {
      printf("fopen exampleopen filed!\n");
      fclose(pFile);
    }


    cout << "====> Model loading ..." << endl << endl;
    // -- Initial positions and velocities --

    //First DS
    SP::SiconosVector q10(new SiconosVector(qDim));
    SP::SiconosVector v10(new SiconosVector(nDim));
    SP::SimpleMatrix I1(new SimpleMatrix(3, 3));
    v10->zero();
    I1->eye();
    I1->setValue(0, 0, 0.1);
    // Initial position of the center of gravity CG1
    (*q10)(0) = 1.0;0.5 * L1 / sqrt(2.0);
    (*q10)(1) = 0;
    (*q10)(2) = -1.0;-0.5 * L1 / sqrt(2.0);
    // Initial orientation (a quaternion that gives the rotation w.r.t the spatial frame)
    // angle of the rotation Pi/4
    double angle = M_PI / 4;
    angle =0.0;
    SiconosVector V1(3);
    V1.zero();
    // vector of the rotation (Y-axis)
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    // construction of the quaternion
    q10->setValue(3, cos(angle / 2));
    q10->setValue(4, V1.getValue(0)*sin(angle / 2));
    q10->setValue(5, V1.getValue(1)*sin(angle / 2));
    q10->setValue(6, V1.getValue(2)*sin(angle / 2));

    // -- The dynamical system --
    SP::NewtonEulerDS beam1(new NewtonEulerDS(q10, v10, m, I1));
    // -- Set external forces (weight) --
    SP::SiconosVector weight(new SiconosVector(nDof));
    (*weight)(2) = -m * g;
    beam1->setFExtPtr(weight);
    //second DS
    SP::SiconosVector q02(new SiconosVector(qDim));
    SP::SiconosVector v02(new SiconosVector(nDim));
    SP::SimpleMatrix I2(new SimpleMatrix(3, 3));
    v02->zero();
    I2->eye();
    I2->setValue(0, 0, 0.1);
    (*q02)(0) = L1 / sqrt(2.0) - 0.5 * L2 / sqrt(2.0);
    (*q02)(1) = 0;
    (*q02)(2) = -L1 / sqrt(2.0) - 0.5 * L2 / sqrt(2.0);

    // angle = -M_PI / 4;
    // V1.zero();
    // V1.setValue(0, 0);
    // V1.setValue(1, 1);
    // V1.setValue(2, 0);
    // q02->setValue(3, cos(angle / 2));
    // q02->setValue(4, V1.getValue(0)*sin(angle / 2));
    // q02->setValue(5, V1.getValue(1)*sin(angle / 2));
    // q02->setValue(6, V1.getValue(2)*sin(angle / 2));

    // SP::NewtonEulerDS beam2(new NewtonEulerDS(q02, v02, m, I2));
    // // -- Set external forces (weight) --
    // SP::SiconosVector weight2(new SiconosVector(nDof));
    // (*weight2)(2) = -m * g;
    // beam2->setFExtPtr(weight2);


    // SP::SiconosVector q03(new SiconosVector(qDim));
    // SP::SiconosVector v03(new SiconosVector(nDim));
    // SP::SimpleMatrix I3(new SimpleMatrix(3, 3));
    // v03->zero();
    // I3->eye();
    // I3->setValue(0, 0, 0.1);
    // q03->zero();
    // (*q03)(2) = -L1 * sqrt(2.0) - L1 / 2;

    // angle = M_PI / 2;
    // V1.zero();
    // V1.setValue(0, 0);
    // V1.setValue(1, 1);
    // V1.setValue(2, 0);
    // q03->setValue(3, cos(angle / 2));
    // q03->setValue(4, V1.getValue(0)*sin(angle / 2));
    // q03->setValue(5, V1.getValue(1)*sin(angle / 2));
    // q03->setValue(6, V1.getValue(2)*sin(angle / 2));

    // SP::NewtonEulerDS beam3(new NewtonEulerDS(q03, v03, m, I3));
    // // -- Set external forces (weight) --
    // SP::SiconosVector weight3(new SiconosVector(nDof));
    // (*weight3)(2) = -m * g;
    // beam3->setFExtPtr(weight3);





    beam1->setComputeFIntFunction("SimplePlugin", "fInt_beam1");
    beam1->setComputeJacobianFIntqFunction("SimplePlugin", "jacobianFIntq_beam1");
    beam1->setComputeJacobianFIntvFunction("SimplePlugin", "jacobianFIntv_beam1");


    // --------------------
    // --- Interactions ---
    // --------------------

    // Interaction with the floor
    double e = 0.9;
    SP::SimpleMatrix H(new SimpleMatrix(1, qDim));
    SP::SiconosVector eR(new SiconosVector(1));
    eR->setValue(0, 2.3);
    H->zero();
    (*H)(0, 2) = 1.0;
    SP::NonSmoothLaw nslaw0(new NewtonImpactNSL(e));
    SP::NewtonEulerR relation0(new NewtonEulerR());
    relation0->setJachq(H);
    relation0->setE(eR);
    cout << "main jacQH" << endl;
    relation0->jachq()->display();


    // Interactions



    //SP::NonSmoothLaw nslaw3(new EqualityConditionNSLKneeJointR::numberOfConstraints()());
    SP::SiconosVector P(new SiconosVector(3));
    P->zero();
    // Building the first knee joint for beam1
    // input  - the concerned DS : beam1
    //        - a point in the spatial frame (absolute frame) where the knee is defined P
    SP::KneeJointR relation1(new KneeJointR(P, true, beam1));

    // // relation4->setJachq(H4);
    SP::NonSmoothLaw nslaw1(new EqualityConditionNSL(relation1->numberOfConstraints()));

    SP::Interaction inter1(new Interaction(nslaw1, relation1));

    // -------------
    // --- Model ---
    // -------------
    SP::Model myModel(new Model(t0, T));
    // add the dynamical system in the non smooth dynamical system
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam1);
    // link the interaction and the dynamical system
    myModel->nonSmoothDynamicalSystem()->link(inter1, beam1);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::MoreauJeanOSI OSI1(new MoreauJeanOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new MLCP());
    osnspb->numericsSolverOptions()->dparam[0]=1e-10;

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI1, osnspb));

//    s->insertIntegrator(OSI2);
    //s->insertIntegrator(OSI3);
    //    s->setComputeResiduY(true);
    //  s->setUseRelativeConvergenceCriteron(false);
    myModel->setSimulation(s);



    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    myModel->initialize();


    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15 + 7;
    SimpleMatrix dataPlot(N, outputSize);
    SimpleMatrix beam1Plot(2,3*N);
    SimpleMatrix beam2Plot(2,3*N);
    SimpleMatrix beam3Plot(2,3*N);

    SP::SiconosVector q1 = beam1->q();
    SP::SiconosVector q2 = beam1->q();
    SP::SiconosVector q3 = beam1->q();
    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    SP::SiconosVector yAux(new SiconosVector(3));
    yAux->setValue(0, 1);
    SP::SimpleMatrix Jaux(new SimpleMatrix(3, 3));
    Index dimIndex(2);
    Index startIndex(4);
    fprintf(pFile, "double T[%d*%d]={", N + 1, outputSize);
    double beamTipTrajectories[6];

    for (k = 0; k < N; k++)
    {
      // solve ...
      s->newtonSolve(1e-10, 4);



      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();

      dataPlot(k, 1) = (*q1)(0);
      dataPlot(k, 2) = (*q1)(1);
      dataPlot(k, 3) = (*q1)(2);
      dataPlot(k, 4) = (*q1)(3);
      dataPlot(k, 5) = (*q1)(4);
      dataPlot(k, 6) = (*q1)(5);
      dataPlot(k, 7) = (*q1)(6);

      dataPlot(k, 8) = (*q2)(0);
      dataPlot(k, 9) = (*q2)(1);
      dataPlot(k, 10) = (*q2)(2);
      dataPlot(k, 11) = (*q2)(3);
      dataPlot(k, 12) = (*q2)(4);
      dataPlot(k, 13) = (*q2)(5);
      dataPlot(k, 14) = (*q2)(6);

      dataPlot(k, 15) = (*q3)(0);
      dataPlot(k, 16) = (*q3)(1);
      dataPlot(k, 17) = (*q3)(2);
      dataPlot(k, 18) = (*q3)(3);
      dataPlot(k, 19) = (*q3)(4);
      dataPlot(k, 20) = (*q3)(5);
      dataPlot(k, 21) = (*q3)(6);

      tipTrajectories(q1,beamTipTrajectories,L1);
      beam1Plot(0,3*k) = beamTipTrajectories[0];
      beam1Plot(0,3*k+1) = beamTipTrajectories[1];
      beam1Plot(0,3*k+2) = beamTipTrajectories[2];
      beam1Plot(1,3*k) = beamTipTrajectories[3];
      beam1Plot(1,3*k+1) = beamTipTrajectories[4];
      beam1Plot(1,3*k+2) = beamTipTrajectories[5];


      //printf("reaction1:%lf \n", interFloor->lambda(1)->getValue(0));

      for (unsigned int jj = 0; jj < outputSize; jj++)
      {
        if ((k || jj))
          fprintf(pFile, ",");
        fprintf(pFile, "%f", dataPlot(k, jj));
      }
      fprintf(pFile, "\n");
      s->nextStep();
      ++show_progress;
    }
    fprintf(pFile, "};");
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    cout << "====> Output file writing ..." << endl;
    ioMatrix::write("NE_1DS_1Knee_MLCP.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write("NE_1DS_1Knee_MLCP_beam1.dat", "ascii", beam1Plot, "noDim");
    
    std::cout << "====> Comparison with reference file ..." << std::endl;

    SimpleMatrix dataPlotRef(dataPlot);
    dataPlotRef.zero();
    ioMatrix::read("NE_1DS_1Knee_MLCP.ref", "ascii", dataPlotRef);
    std::cout << "Error w.r.t. reference file : " << (dataPlot - dataPlotRef).normInf() << std::endl;
    if ((dataPlot - dataPlotRef).normInf() > 1e-7)
    {
      (dataPlot - dataPlotRef).display();
      std::cout << "Warning. The results is rather different from the reference file." << std::endl;
      return 1;
    }

    fclose(pFile);
  }

  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in NE_...cpp" << endl;
  }

}

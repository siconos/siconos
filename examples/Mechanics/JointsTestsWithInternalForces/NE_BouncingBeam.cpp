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
    SP::SiconosVector q03(new SiconosVector(qDim));
    SP::SiconosVector v03(new SiconosVector(nDim));
    SP::SimpleMatrix I3(new SimpleMatrix(3, 3));
    v03->zero();
    I3->eye();
    I3->setValue(0, 0, 0.1);
    q03->zero();
    (*q03)(2) = -L1 * sqrt(2.0) - L1 / 2;

    double angle = M_PI / 2;
    SiconosVector V1(3);
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q03->setValue(3, cos(angle / 2));
    q03->setValue(4, V1.getValue(0)*sin(angle / 2));
    q03->setValue(5, V1.getValue(1)*sin(angle / 2));
    q03->setValue(6, V1.getValue(2)*sin(angle / 2));

    SP::NewtonEulerDS bouncingbeam(new NewtonEulerDS(q03, v03, m, I3));
    // -- Set external forces (weight) --
    SP::SiconosVector weight3(new SiconosVector(nDof));
    (*weight3)(2) = -m * g;
    bouncingbeam->setFExtPtr(weight3);
    bouncingbeam->setComputeFIntFunction("SimplePlugin", "fInt_beam1");
    bouncingbeam->setComputeJacobianFIntqFunction("SimplePlugin", "jacobianFIntq_beam1");
    bouncingbeam->setComputeJacobianFIntvFunction("SimplePlugin", "jacobianFIntv_beam1");


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



     SP::SiconosVector axe1(new SiconosVector(3));
    axe1->zero();
    axe1->setValue(2, 1);

    SP::PrismaticJointR relation4(new PrismaticJointR(axe1, false, bouncingbeam));
    SP::NonSmoothLaw nslaw4(new EqualityConditionNSL(relation4->numberOfConstraints()));
    SP::Interaction inter4(new Interaction(nslaw4, relation4));
    SP::Interaction interFloor(new Interaction(nslaw0, relation0));

    // -------------
    // --- Model ---
    // -------------
    SP::Model myModel(new Model(t0, T));
    // add the dynamical system in the non smooth dynamical system
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(bouncingbeam);
    // link the interaction and the dynamical system

    myModel->nonSmoothDynamicalSystem()->link(inter4, bouncingbeam);
    myModel->nonSmoothDynamicalSystem()->link(interFloor, bouncingbeam);
    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --

    SP::MoreauJeanOSI OSI3(new MoreauJeanOSI(theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new MLCP());

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI3, osnspb));
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
    SimpleMatrix bouncingbeamPlot(2,3*N);

    SP::SiconosVector q3 = bouncingbeam->q();
    SP::SiconosVector y= interFloor->y(0);
    SP::SiconosVector ydot= interFloor->y(1);

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
      s->newtonSolve(1e-4, 50);



      // --- Get values to be plotted ---
      dataPlot(k, 0) =  s->nextTime();

      dataPlot(k, 1) = (*q3)(0);
      dataPlot(k, 2) = (*q3)(1);
      dataPlot(k, 3) = (*q3)(2);
      dataPlot(k, 4) = (*q3)(3);
      dataPlot(k, 5) = (*q3)(4);
      dataPlot(k, 6) = (*q3)(5);
      dataPlot(k, 7) = (*q3)(6);

      dataPlot(k, 8) = y->norm2();
      dataPlot(k, 9) = ydot->norm2();



      tipTrajectories(q3,beamTipTrajectories,L3);
      bouncingbeamPlot(0,3*k) = beamTipTrajectories[0];
      bouncingbeamPlot(0,3*k+1) = beamTipTrajectories[1];
      bouncingbeamPlot(0,3*k+2) = beamTipTrajectories[2];
      bouncingbeamPlot(1,3*k) = beamTipTrajectories[3];
      bouncingbeamPlot(1,3*k+1) = beamTipTrajectories[4];
      bouncingbeamPlot(1,3*k+2) = beamTipTrajectories[5];

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
    ioMatrix::write("NE_BouncingBeam.dat", "ascii", dataPlot, "noDim");
    ioMatrix::write("NE_BouncingBeam_beam.dat", "ascii", bouncingbeamPlot, "noDim");

    // SimpleMatrix dataPlotRef(dataPlot);
    // dataPlotRef.zero();
    // ioMatrix::read("NE_BoundingBeam.ref", "ascii", dataPlotRef);
    // if ((dataPlot - dataPlotRef).normInf() > 1e-7)
    // {
    //   (dataPlot - dataPlotRef).display();
    //   std::cout << "Warning. The results is rather different from the reference file." << std::endl;
    //   return 1;
    // }

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

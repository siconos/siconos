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

/*!\file NE....cpp
  \brief \ref EMNE_MULTIBIDY - C++ input file, Time-Stepping version - O.B.

  A multiby example.
  Direct description of the model without XML input.
  Simulation with a Time-Stepping scheme.
*/

#include "SiconosKernel.hpp"
#include "../../../Kernel/src/multiBodyTools/joints/KneeJoint.hpp"
#include "../../../Kernel/src/multiBodyTools/joints/PrismaticJoint.hpp"
using namespace std;

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
    double L1 = 1.0;
    double L2 = 1.0;
    double L3 = 1.0;
    double theta = 1.0;              // theta for Moreau integrator
    double g = 9.81; // Gravity
    double m = 1.;
    double wx = 0.0;
    double wz = 0.0;
    double wy = 0.0;
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
    DynamicalSystemsSet allDS1; // the list of DS
    DynamicalSystemsSet allDS2; // the list of DS
    DynamicalSystemsSet allDS3; // the list of DS
    DynamicalSystemsSet allDS4; // the list of DS
    DynamicalSystemsSet allDS_With_Floor;

    // -- Initial positions and velocities --

    //First DS
    SP::SimpleVector q10(new SimpleVector(qDim));
    SP::SimpleVector v10(new SimpleVector(nDim));
    SP::SimpleMatrix I1(new SimpleMatrix(3, 3));
    v10->zero();
    I1->eye();
    I1->setValue(0, 0, 0.1);
    (*q10)(0) = 0.5 * L1 / sqrt(2);
    (*q10)(1) = 0;
    (*q10)(2) = -0.5 * L1 / sqrt(2);
    double angle = M_PI / 4;
    SimpleVector V1(3);
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q10->setValue(3, cos(angle / 2));
    q10->setValue(4, V1.getValue(0)*sin(angle / 2));
    q10->setValue(5, V1.getValue(1)*sin(angle / 2));
    q10->setValue(6, V1.getValue(2)*sin(angle / 2));
    // -- The dynamical system --
    SP::NewtonEulerDS beam1(new NewtonEulerDS(q10, v10, m, I1));
    allDS1.insert(beam1);
    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof));
    (*weight)(2) = -m * g;
    beam1->setFExtPtr(weight);


    //second DS
    SP::SimpleVector q02(new SimpleVector(qDim));
    SP::SimpleVector v02(new SimpleVector(nDim));
    SP::SimpleMatrix I2(new SimpleMatrix(3, 3));
    v02->zero();
    I2->eye();
    I2->setValue(0, 0, 0.1);
    (*q02)(0) = L1 / sqrt(2) - 0.5 * L2 / sqrt(2);
    (*q02)(1) = 0;
    (*q02)(2) = -L1 / sqrt(2) - 0.5 * L2 / sqrt(2);

    angle = -M_PI / 4;
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q02->setValue(3, cos(angle / 2));
    q02->setValue(4, V1.getValue(0)*sin(angle / 2));
    q02->setValue(5, V1.getValue(1)*sin(angle / 2));
    q02->setValue(6, V1.getValue(2)*sin(angle / 2));

    SP::NewtonEulerDS beam2(new NewtonEulerDS(q02, v02, m, I2));
    allDS2.insert(beam1);
    allDS2.insert(beam2);
    // -- Set external forces (weight) --
    SP::SimpleVector weight2(new SimpleVector(nDof));
    (*weight2)(2) = -m * g;
    beam2->setFExtPtr(weight2);


    SP::SimpleVector q03(new SimpleVector(qDim));
    SP::SimpleVector v03(new SimpleVector(nDim));
    SP::SimpleMatrix I3(new SimpleMatrix(3, 3));
    v03->zero();
    I3->eye();
    I3->setValue(0, 0, 0.1);
    q03->zero();
    (*q03)(2) = -L1 * sqrt(2) - L1 / 2;

    angle = M_PI / 2;
    V1.zero();
    V1.setValue(0, 0);
    V1.setValue(1, 1);
    V1.setValue(2, 0);
    q03->setValue(3, cos(angle / 2));
    q03->setValue(4, V1.getValue(0)*sin(angle / 2));
    q03->setValue(5, V1.getValue(1)*sin(angle / 2));
    q03->setValue(6, V1.getValue(2)*sin(angle / 2));

    SP::NewtonEulerDS beam3(new NewtonEulerDS(q03, v03, m, I3));
    allDS3.insert(beam2);
    allDS3.insert(beam3);
    // -- Set external forces (weight) --
    SP::SimpleVector weight3(new SimpleVector(nDof));
    (*weight3)(2) = -m * g;
    beam3->setFExtPtr(weight3);
    allDS_With_Floor.insert(beam3);




    // --------------------
    // --- Interactions ---
    // --------------------

    InteractionsSet allInteractions;


    // Interaction with the floor
    double e = 0.9;
    SP::SiconosMatrix H(new SimpleMatrix(1, qDim));
    SP::SimpleVector eR(new SimpleVector(1));
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
    //
    SP::SiconosMatrix H1(new SimpleMatrix(KneeJointR::_sNbEqualities, qDim));
    H1->zero();
    SP::SiconosMatrix H2(new SimpleMatrix(KneeJointR::_sNbEqualities, 2 * qDim));
    SP::SiconosMatrix H3(new SimpleMatrix(KneeJointR::_sNbEqualities, 2 * qDim));
    H2->zero();
    H3->zero();
    SP::NonSmoothLaw nslaw1(new EqualityConditionNSL(KneeJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw2(new EqualityConditionNSL(KneeJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw3(new EqualityConditionNSL(KneeJointR::_sNbEqualities));

    //SP::NonSmoothLaw nslaw3(new EqualityConditionNSLKneeJointR::_sNbEqualities());
    SP::SimpleVector P(new SimpleVector(3));
    P->zero();

    SP::NewtonEulerR relation1(new KneeJointR(beam1, P));

    SP::SimpleVector G20(new SimpleVector(3));
    P->zero();
    P->setValue(0, L1 / 2);
    SP::NewtonEulerR relation2(new KneeJointR(beam1, beam2, P));
    P->zero();
    P->setValue(0, -L1 / 2);
    SP::NewtonEulerR relation3(new KneeJointR(beam2, beam3, P));


    //relation Prismatic

    SP::SiconosMatrix H4(new SimpleMatrix(PrismaticJointR::_sNbEqualities, qDim));
    H4->zero();
    SP::NonSmoothLaw nslaw4(new EqualityConditionNSL(PrismaticJointR::_sNbEqualities));
    SP::SimpleVector axe1(new SimpleVector(3));
    axe1->zero();
    axe1->setValue(2, 1);
    SP::NewtonEulerR relation4(new PrismaticJointR(beam3, axe1));
    allDS4.insert(beam3);

    relation1->setJachq(H1);
    relation2->setJachq(H2);
    relation3->setJachq(H3);
    relation4->setJachq(H4);

    SP::Interaction inter1(new Interaction("axis-beam1", allDS1, 0, KneeJointR::_sNbEqualities, nslaw1, relation1));
    allInteractions.insert(inter1);
    SP::Interaction inter2(new Interaction("axis-beam2", allDS2, 1, KneeJointR::_sNbEqualities, nslaw2, relation2));
    allInteractions.insert(inter2);
    SP::Interaction inter3(new Interaction("axis-beam3", allDS3, 1, KneeJointR::_sNbEqualities, nslaw3, relation3));
    allInteractions.insert(inter3);
    SP::Interaction inter4(new Interaction("axis-beam4", allDS4, 1, PrismaticJointR::_sNbEqualities, nslaw4, relation4));
    allInteractions.insert(inter4);
    SP::Interaction interFloor(new Interaction("floor-ball", allDS_With_Floor , 0, 1, nslaw0, relation0));
    allInteractions.insert(interFloor);




    // -------------
    // --- Model ---
    // -------------
    SP::Model myModel(new Model(t0, T, allDS2, allInteractions));
    // add the dynamical system in the non smooth dynamical system
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam1);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam2);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam3);
    // link the interaction and the dynamical system
    myModel->nonSmoothDynamicalSystem()->link(inter1, beam1);

    myModel->nonSmoothDynamicalSystem()->link(inter2, beam1);
    myModel->nonSmoothDynamicalSystem()->link(inter2, beam2);

    myModel->nonSmoothDynamicalSystem()->link(inter3, beam2);
    myModel->nonSmoothDynamicalSystem()->link(inter3, beam3);

    myModel->nonSmoothDynamicalSystem()->link(inter4, beam3);

    myModel->nonSmoothDynamicalSystem()->link(interFloor, beam3);
    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::Moreau OSI1(new Moreau(beam1, theta));
    SP::Moreau OSI2(new Moreau(beam2, theta));
    SP::Moreau OSI3(new Moreau(beam3, theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new MLCP(SICONOS_MLCP_PATH));

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI1, osnspb));
    s->insertIntegrator(OSI2);
    s->insertIntegrator(OSI3);
    //    s->setComputeResiduY(true);
    //  s->setUseRelativeConvergenceCriteron(false);



    // =========================== End of model definition ===========================

    // ================================= Computation =================================

    // --- Simulation initialization ---

    cout << "====> Initialisation ..." << endl << endl;
    myModel->initialize(s);
    int N = (int)((T - t0) / h); // Number of time steps


    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    unsigned int outputSize = 15 + 7;
    SimpleMatrix dataPlot(N + 1, outputSize);

    SP::SiconosVector q1 = beam1->q();
    SP::SiconosVector q2 = beam2->q();
    SP::SiconosVector q3 = beam3->q();
    std::cout << "computeH1\n";
    relation1->computeh(0.);
    std::cout << "computeH2\n";
    relation2->computeh(0.);
    std::cout << "computeH3\n";
    relation3->computeh(0.);

    std::cout << "computeH0\n";
    relation0->computeh(0.);
    //    std::cout<<"computeH4\n";
    //    relation4->computeh(0.);

    // --- Time loop ---
    cout << "====> Start computation ... " << endl << endl;
    // ==== Simulation loop - Writing without explicit event handling =====
    int k = 0;
    boost::progress_display show_progress(N);

    boost::timer time;
    time.restart();
    SP::SimpleVector yAux(new SimpleVector(3));
    yAux->setValue(0, 1);
    SP::SimpleMatrix Jaux(new SimpleMatrix(3, 3));
    int NewtonIt = 0;
    Index dimIndex(2);
    Index startIndex(4);
    fprintf(pFile, "double T[%d*%d]={", N + 1, outputSize);
    while (s->nextTime() < T)
    {
      // solve ...
      s->newtonSolve(1e-4, 50);

      //   relation0->computeH(0.);
      //      yAux->setValue(0,inter->y(0)->getValue(0));
      //     yAux->setValue(1,inter->y(0)->getValue(1));
      //     yAux->setValue(2,inter->y(0)->getValue(2));

      //  NewtonIt=0;
      //  double aux = yAux->norm2();
      //  if (aux>0)
      //  while (NewtonIt <1 && yAux->norm2()>1e-3){
      //    relation0->computeJacQH(0.);
      //    Jaux->zero();
      //    startIndex[0]=0;
      //    startIndex[1]=0;
      //    startIndex[2]=0;
      //    startIndex[3]=0;
      //    dimIndex[0]=3;
      //    dimIndex[1]=3;
      //    setBlock(H,Jaux,dimIndex,startIndex);
      //    //    Jaux->display();
      //    Jaux->PLUForwardBackwardInPlace(*yAux);
      //    q->setValue(0,-yAux->getValue(0)+(*q)(0));
      //    q->setValue(1,-yAux->getValue(1)+(*q)(1));
      //    q->setValue(2,-yAux->getValue(2)+(*q)(2));
      //    NewtonIt++;
      //    relation0->computeH(0.);
      //    yAux->setValue(0,inter->y(0)->getValue(0));
      //    yAux->setValue(1,inter->y(0)->getValue(1));
      //    yAux->setValue(2,inter->y(0)->getValue(2));
      //    aux = yAux->norm2();
      //  }


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
      printf("reaction0:%lf", interFloor->lambda(0)->getValue(0));
      printf("reaction1:%lf \n", interFloor->lambda(1)->getValue(0));

      for (int jj = 0; jj < outputSize; jj++)
      {
        if ((k || jj))
          fprintf(pFile, ",");
        fprintf(pFile, "%f", dataPlot(k, jj));
      }
      fprintf(pFile, "\n");
      s->nextStep();
      ++show_progress;
      k++;
    }
    fprintf(pFile, "};");
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    // --- Output files ---
    //cout<<"====> Output file writing ..."<<endl;
    //ioMatrix io("result.dat", "ascii");
    //io.write(dataPlot,"noDim");
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

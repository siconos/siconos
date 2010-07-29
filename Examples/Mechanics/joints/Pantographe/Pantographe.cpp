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

/*COMPILATION WITH SICONOS:
siconos --opt -I/usr/include/qt4/Qt --opt -I/usr/include/qt4 --opt -I/usr/include/qt4/QtOpenGL --opt -I/usr/include/qt4/QtXml --opt -I/usr/include/qt4/QtCore --opt -I/usr/include/qt4/QtGui -L/usr/X11R6/lib -L/usr/lib -L/usr/lib -lQGLViewer -lpthread -lGLU -lGL -lQtXml -lQtOpenGL -lQtGui -lQtCore -v -g Pantographe.cpp
*/



#include "SiconosKernel.hpp"
#include "KneeJoint.hpp"
#include "PivotJoint.hpp"
#include "PrismaticJoint.hpp"
#include <QGLViewer/qglviewer.h>
using namespace std;
#include "qgl.h"
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
    double T = 3.0;                  // final computation time
    double h = 0.0005;                // time step
    int Freq = 20;

    double theta = 1.0;              // theta for Moreau integrator
    double g = 9.81; // Gravity
    double m = 1.;
    double wx = 0.0;
    double wz = 0.0;
    double wy = 0.0;
    double Angle = atan(0.5);
    double L = 6 / cos(Angle);
    double L1 = 2 * sqrt(2);
    double L2 = sqrt(2);
    // -------------------------
    // --- Dynamical systems ---
    // -------------------------




    cout << "====> Model loading ..." << endl << endl;
    DynamicalSystemsSet allDS, allDS14, allDS15, allDS23, allDS26, allDS56, allDS34;
    DynamicalSystemsSet allDS1; // the list of DS
    DynamicalSystemsSet allDS2; // the list of DS
    DynamicalSystemsSet allDS3; // the list of DS
    DynamicalSystemsSet allDS4; // the list of DS


    // -- Initial positions and velocities --

    //First DS
    SP::SimpleVector q10(new SimpleVector(qDim));
    SP::SimpleVector v10(new SimpleVector(nDim));
    SP::SimpleMatrix I1(new SimpleMatrix(3, 3));
    v10->zero();
    I1->eye();
    I1->setValue(0, 0, 0.1);
    (*q10)(0) = -1;
    (*q10)(1) = 0.5;
    (*q10)(2) = 0;
    (*q10)(3) = cos(Angle / 2.0);
    (*q10)(4) = 0.;
    (*q10)(5) = 0.;
    (*q10)(6) = -sin(Angle / 2.0);

    // -- The dynamical system --
    SP::NewtonEulerDS beam1(new NewtonEulerDS(q10, v10, m, I1));

    // -- Set external forces (weight) --
    SP::SimpleVector weight(new SimpleVector(nDof));
    (*weight)(1) = -m * g;
    beam1->setFExtPtr(weight);


    //second DS
    SP::SimpleVector q20(new SimpleVector(qDim));
    SP::SimpleVector v20(new SimpleVector(nDim));
    SP::SimpleMatrix I2(new SimpleMatrix(3, 3));
    v20->zero();
    I2->eye();
    I2->setValue(0, 0, 0.1);
    (*q20)(0) = -1;
    (*q20)(1) = -0.5;
    (*q20)(2) = 0;
    (*q20)(3) = cos(Angle / 2.0);
    (*q20)(4) = 0.;
    (*q20)(5) = 0.;
    (*q20)(6) = sin(Angle / 2.0);
    SP::NewtonEulerDS beam2(new NewtonEulerDS(q20, v20, m, I2));
    // -- Set external forces (weight) --
    SP::SimpleVector weight2(new SimpleVector(nDof));
    (*weight2)(1) = -m * g;
    beam2->setFExtPtr(weight2);

    //DS 3
    SP::SimpleVector q30(new SimpleVector(qDim));
    SP::SimpleVector v30(new SimpleVector(nDim));
    SP::SimpleMatrix I3(new SimpleMatrix(3, 3));
    v30->zero();
    I3->eye();
    I3->setValue(0, 0, 0.1);
    q30->zero();
    (*q30)(0) = 2.5;
    (*q30)(1) = 0.5;
    (*q30)(2) = 0;
    (*q30)(3) = cos(M_PI / 8.0);
    (*q30)(4) = 0.;
    (*q30)(5) = 0.;
    (*q30)(6) = -sin(M_PI / 8.0);
    SP::NewtonEulerDS beam3(new NewtonEulerDS(q30, v30, m, I3));
    // -- Set external forces (weight) --
    SP::SimpleVector weight3(new SimpleVector(nDof));
    (*weight3)(1) = -m * g;
    beam3->setFExtPtr(weight3);


    //DS 4
    SP::SimpleVector q40(new SimpleVector(qDim));
    SP::SimpleVector v40(new SimpleVector(nDim));
    SP::SimpleMatrix I4(new SimpleMatrix(3, 3));
    v40->zero();
    I4->eye();
    I4->setValue(0, 0, 0.1);
    q40->zero();
    (*q40)(0) = 2.5;
    (*q40)(1) = -0.5;
    (*q40)(2) = 0;
    (*q40)(3) = cos(M_PI / 8.0);
    (*q40)(4) = 0.;
    (*q40)(5) = 0.;
    (*q40)(6) = sin(M_PI / 8.0);
    SP::NewtonEulerDS beam4(new NewtonEulerDS(q40, v40, m, I4));
    // -- Set external forces (weight) --
    SP::SimpleVector weight4(new SimpleVector(nDof));
    (*weight4)(1) = -m * g;
    beam4->setFExtPtr(weight4);


    //DS 5
    SP::SimpleVector q50(new SimpleVector(qDim));
    SP::SimpleVector v50(new SimpleVector(nDim));
    SP::SimpleMatrix I5(new SimpleMatrix(3, 3));
    v50->zero();
    I5->eye();
    I5->setValue(0, 0, 0.1);
    q50->zero();
    (*q50)(0) = -5;
    (*q50)(1) = 1;
    (*q50)(2) = 0;
    (*q50)(3) = cos(M_PI / 8.0);
    (*q50)(4) = 0.;
    (*q50)(5) = 0.;
    (*q50)(6) = sin(M_PI / 8.0);
    SP::NewtonEulerDS beam5(new NewtonEulerDS(q50, v50, m, I5));
    // -- Set external forces (weight) --
    SP::SimpleVector weight5(new SimpleVector(nDof));
    (*weight5)(1) = -m * g;
    beam5->setFExtPtr(weight5);


    //DS 6
    SP::SimpleVector q60(new SimpleVector(qDim));
    SP::SimpleVector v60(new SimpleVector(nDim));
    SP::SimpleMatrix I6(new SimpleMatrix(3, 3));
    v60->zero();
    I6->eye();
    I6->setValue(0, 0, 0.1);
    q60->zero();
    (*q60)(0) = -5;
    (*q60)(1) = -1;
    (*q60)(2) = 0;
    (*q60)(3) = cos(M_PI / 8.0);
    (*q60)(4) = 0.;
    (*q60)(5) = 0.;
    (*q60)(6) = -sin(M_PI / 8.0);
    SP::NewtonEulerDS beam6(new NewtonEulerDS(q60, v60, m, I6));
    // -- Set external forces (weight) --
    SP::SimpleVector weight6(new SimpleVector(nDof));
    (*weight6)(1) = -m * g;
    beam6->setFExtPtr(weight6);


    // --------------------
    // --- Interactions ---
    // --------------------

    InteractionsSet allInteractions;





    // Interactions
    //
    SP::SiconosMatrix H1(new SimpleMatrix(KneeJointR::_sNbEqualities, qDim));
    SP::SiconosMatrix H2(new SimpleMatrix(KneeJointR::_sNbEqualities, qDim));
    SP::SiconosMatrix H14(new SimpleMatrix(PivotJointR::_sNbEqualities, 2 * qDim));
    SP::SiconosMatrix H15(new SimpleMatrix(PivotJointR::_sNbEqualities, 2 * qDim));
    SP::SiconosMatrix H23(new SimpleMatrix(PivotJointR::_sNbEqualities, 2 * qDim));
    SP::SiconosMatrix H26(new SimpleMatrix(PivotJointR::_sNbEqualities, 2 * qDim));
    SP::SiconosMatrix H34(new SimpleMatrix(PivotJointR::_sNbEqualities, 2 * qDim));
    SP::SiconosMatrix H56(new SimpleMatrix(PivotJointR::_sNbEqualities, 2 * qDim));
    H1->zero();
    H2->zero();
    H23->zero();
    H14->zero();
    H26->zero();
    H15->zero();
    H34->zero();
    H56->zero();



    //SP::NonSmoothLaw nslaw3(new EqualityConditionNSL(3));
    SP::SimpleVector P1(new SimpleVector(3));
    SP::SimpleVector P2(new SimpleVector(3));
    SP::SimpleVector P14(new SimpleVector(3));
    SP::SimpleVector P15(new SimpleVector(3));
    SP::SimpleVector P23(new SimpleVector(3));
    SP::SimpleVector P26(new SimpleVector(3));
    SP::SimpleVector P56(new SimpleVector(3));
    SP::SimpleVector P34(new SimpleVector(3));
    SP::SimpleVector Axe(new SimpleVector(3));
    P1->zero();
    P2->zero();
    P14->zero();
    P15->zero();
    P23->zero();
    P26->zero();
    P34->zero();
    P56->zero();
    Axe->zero();
    P14->setValue(0, L / 2.0);
    P15->setValue(0, -L / 2.0);
    P23->setValue(0, L / 2.0);
    P26->setValue(0, -L / 2.0);
    P34->setValue(0, 0.5 * sqrt(2));
    P56->setValue(0, -sqrt(2));
    Axe->setValue(2, 1.0);
    //    P14->setValue(1,-1);
    SP::NewtonEulerR relation1(new KneeJointR(beam1, P1));
    SP::NewtonEulerR relation2(new KneeJointR(beam2, P2));
    SP::NewtonEulerR relation14(new PivotJointR(beam1, beam4, P14, Axe));
    SP::NewtonEulerR relation15(new PivotJointR(beam1, beam5, P15, Axe));
    SP::NewtonEulerR relation23(new PivotJointR(beam2, beam3, P23, Axe));
    SP::NewtonEulerR relation26(new PivotJointR(beam2, beam6, P26, Axe));
    SP::NewtonEulerR relation34(new PivotJointR(beam3, beam4, P34, Axe));
    SP::NewtonEulerR relation56(new PivotJointR(beam5, beam6, P56, Axe));

    SP::NonSmoothLaw nslaw1(new EqualityConditionNSL(KneeJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw2(new EqualityConditionNSL(KneeJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw14(new EqualityConditionNSL(PivotJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw15(new EqualityConditionNSL(PivotJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw23(new EqualityConditionNSL(PivotJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw26(new EqualityConditionNSL(PivotJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw34(new EqualityConditionNSL(PivotJointR::_sNbEqualities));
    SP::NonSmoothLaw nslaw56(new EqualityConditionNSL(PivotJointR::_sNbEqualities));

    relation1->setJachq(H1);
    relation2->setJachq(H2);
    relation14->setJachq(H14);
    relation15->setJachq(H15);
    relation23->setJachq(H23);
    relation26->setJachq(H26);
    relation34->setJachq(H34);
    relation56->setJachq(H56);


    allDS1.insert(beam1);
    SP::Interaction inter1(new Interaction("axis-beam1", allDS1, 0, KneeJointR::_sNbEqualities, nslaw1, relation1));
    allInteractions.insert(inter1);

    allDS2.insert(beam2);
    SP::Interaction inter2(new Interaction("axis-beam2", allDS2, 1, KneeJointR::_sNbEqualities, nslaw2, relation2));
    allInteractions.insert(inter2);

    allDS14.insert(beam1);
    allDS14.insert(beam4);
    SP::Interaction inter14(new Interaction("axis-beam14", allDS14, 2, PivotJointR::_sNbEqualities, nslaw14, relation14));
    allInteractions.insert(inter14);

    allDS15.insert(beam1);
    allDS15.insert(beam5);
    SP::Interaction inter15(new Interaction("axis-beam15", allDS15, 3, PivotJointR::_sNbEqualities, nslaw15, relation15));
    allInteractions.insert(inter15);

    allDS23.insert(beam2);
    allDS23.insert(beam3);
    SP::Interaction inter23(new Interaction("axis-beam23", allDS23, 4, PivotJointR::_sNbEqualities, nslaw23, relation23));
    allInteractions.insert(inter23);

    allDS26.insert(beam2);
    allDS26.insert(beam6);
    SP::Interaction inter26(new Interaction("axis-beam26", allDS26, 5, PivotJointR::_sNbEqualities, nslaw26, relation26));
    allInteractions.insert(inter26);

    allDS26.insert(beam3);
    allDS26.insert(beam4);
    SP::Interaction inter34(new Interaction("axis-beam34", allDS34, 6, PivotJointR::_sNbEqualities, nslaw34, relation34));
    allInteractions.insert(inter34);

    allDS26.insert(beam5);
    allDS26.insert(beam6);
    SP::Interaction inter56(new Interaction("axis-beam56", allDS56, 7, PivotJointR::_sNbEqualities, nslaw56, relation56));
    allInteractions.insert(inter56);

    allDS.insert(beam1);
    allDS.insert(beam2);
    allDS.insert(beam3);
    allDS.insert(beam4);
    allDS.insert(beam5);
    allDS.insert(beam6);



    // -------------
    // --- Model ---
    // -------------
    SP::Model myModel(new Model(t0, T, allDS, allInteractions));
    // add the dynamical system in the non smooth dynamical system
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam1);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam2);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam3);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam4);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam5);
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(beam6);
    // link the interaction and the dynamical system
    myModel->nonSmoothDynamicalSystem()->link(inter1, beam1);
    myModel->nonSmoothDynamicalSystem()->link(inter2, beam2);


    myModel->nonSmoothDynamicalSystem()->link(inter14, beam1);
    myModel->nonSmoothDynamicalSystem()->link(inter14, beam4);

    myModel->nonSmoothDynamicalSystem()->link(inter15, beam1);
    myModel->nonSmoothDynamicalSystem()->link(inter15, beam5);

    myModel->nonSmoothDynamicalSystem()->link(inter23, beam2);
    myModel->nonSmoothDynamicalSystem()->link(inter23, beam3);

    myModel->nonSmoothDynamicalSystem()->link(inter26, beam2);
    myModel->nonSmoothDynamicalSystem()->link(inter26, beam6);

    myModel->nonSmoothDynamicalSystem()->link(inter34, beam3);
    myModel->nonSmoothDynamicalSystem()->link(inter34, beam4);

    myModel->nonSmoothDynamicalSystem()->link(inter56, beam5);
    myModel->nonSmoothDynamicalSystem()->link(inter56, beam6);

    // ------------------
    // --- Simulation ---
    // ------------------

    // -- (1) OneStepIntegrators --
    SP::Moreau OSI1(new Moreau(beam1, theta));
    SP::Moreau OSI2(new Moreau(beam2, theta));
    SP::Moreau OSI3(new Moreau(beam3, theta));
    SP::Moreau OSI4(new Moreau(beam4, theta));
    SP::Moreau OSI5(new Moreau(beam5, theta));
    SP::Moreau OSI6(new Moreau(beam6, theta));

    // -- (2) Time discretisation --
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));

    // -- (3) one step non smooth problem
    SP::OneStepNSProblem osnspb(new MLCP("PATH"));

    // -- (4) Simulation setup with (1) (2) (3)
    SP::TimeStepping s(new TimeStepping(t, OSI1, osnspb));
    s->insertIntegrator(OSI2);
    s->insertIntegrator(OSI3);
    s->insertIntegrator(OSI4);
    s->insertIntegrator(OSI5);
    s->insertIntegrator(OSI6);
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


    SP::SiconosVector q[6];
    q[0] = beam1->q();
    q[1] = beam2->q();
    q[2] = beam3->q();
    q[3] = beam4->q();
    q[4] = beam5->q();
    q[5] = beam6->q();


    std::cout << "computeH1\n";
    relation1->computeh(0.);
    std::cout << "computeH2\n";
    relation2->computeh(0.);
    std::cout << "computeH14\n";
    relation14->computeh(0.);

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
    SimpleMatrix dataPlot(N + 2, 7 * 6 + 1);

    dataPlot(k, 0) =  0.;
    for (int num = 0; num < 6; num++)
    {
      for (int coord = 0; coord < 7; coord++)
        dataPlot(k, coord + 1 + num * 7) = (*(q[num]))(coord);
    }
    //    std::cout<<"t "<<dataPlot.getValue(k,0)<<" qo "<<dataPlot.getValue(k,1)<<" "<<dataPlot.getValue(k,2)<<" "<<dataPlot.getValue(k,3)<<" "<<dataPlot.getValue(k,4)<<" "<<dataPlot.getValue(k,5)<<" "<<dataPlot.getValue(k,6)<<" "<<dataPlot.getValue(k,7)<<" \n";
    //    std::cout<<"t "<<dataPlot.getValue(k,0)<<" q1 "<<dataPlot.getValue(k,8)<<" "<<dataPlot.getValue(k,9)<<" "<<dataPlot.getValue(k,10)<<" "<<dataPlot.getValue(k,11)<<" "<<dataPlot.getValue(k,12)<<" "<<dataPlot.getValue(k,13)<<" "<<dataPlot.getValue(k,14)<<" \n";
    k++;

    while (s->nextTime() < T)
    {
      // solve ...
      s->newtonSolve(1e-4, 50);

      dataPlot(k, 0) =  s->nextTime();
      for (int num = 0; num < 6; num++)
      {
        for (int coord = 0; coord < 7; coord++)
          dataPlot(k, coord + 1 + num * 7) = (*(q[num]))(coord);
      }
      //  std::cout<<"t "<<dataPlot.getValue(k,0)<<" qo "<<dataPlot.getValue(k,1)<<" "<<dataPlot.getValue(k,2)<<" "<<dataPlot.getValue(k,3)<<" "<<dataPlot.getValue(k,4)<<" "<<dataPlot.getValue(k,5)<<" "<<dataPlot.getValue(k,6)<<" "<<dataPlot.getValue(k,7)<<" \n";
      // std::cout<<"t "<<dataPlot.getValue(k,0)<<" q1 "<<dataPlot.getValue(k,8)<<" "<<dataPlot.getValue(k,9)<<" "<<dataPlot.getValue(k,10)<<" "<<dataPlot.getValue(k,11)<<" "<<dataPlot.getValue(k,12)<<" "<<dataPlot.getValue(k,13)<<" "<<dataPlot.getValue(k,14)<<" \n";
      //      --- Get values to be plotted ---



      s->nextStep();
      ++show_progress;
      k++;
    }

    //fprintf(pFile,"};");
    cout << endl << "End of computation - Number of iterations done: " << k - 1 << endl;
    cout << "Computation Time " << time.elapsed()  << endl;

    FILE * pFile;
    pFile = fopen("data.wrl", "w");

    fprintf(pFile, "DEF Animation Group {\n  children [\n");
    fprintf(pFile, "\n");
    for (int num = 0; num < 6; num++)
    {
      fprintf(pFile, "DEF Solid%dRotInterp OrientationInterpolator {\n", num + 1);
      fprintf(pFile, "    key [ 0 ");
      for (int cmp = 1; cmp < N; cmp++)
      {
        if (cmp % Freq)
          continue;
        double dcmp = cmp;
        double dN = N;
        fprintf(pFile, ",%e", dcmp / dN);
      }
      fprintf(pFile, "]\n");
      fprintf(pFile, "    keyValue [  ");
      for (int cmp = 0; cmp < N; cmp++)
      {
        if (cmp % Freq)
          continue;

        qglviewer::Quaternion Q1(dataPlot.getValue(cmp, 5 + 7 * num), dataPlot.getValue(cmp, 6 + 7 * num), dataPlot.getValue(cmp, 7 + 7 * num), dataPlot.getValue(cmp, 4 + 7 * num));

        qglviewer::Vec V;
        float A;
        Q1.getAxisAngle(V, A);
        //std::cout<<"V "<<V<<"\n";
        //std::cout<<"V "<<V[0]<<" "<<V[1]<<" "<<V[2]<<"\n";
        fprintf(pFile, "%e %e %e %e,\n", V[0], V[1], V[2], A);
      }
      fprintf(pFile, "  ]\n}\n");



      fprintf(pFile, "DEF Solid%dPosInterp PositionInterpolator {\n", num + 1);
      fprintf(pFile, "    key [ 0 ");
      for (int cmp = 1; cmp < N; cmp++)
      {
        if (cmp % Freq)
          continue;
        double dcmp = cmp;
        double dN = N;
        fprintf(pFile, ",%e", dcmp / dN);
      }
      fprintf(pFile, "]\n");
      fprintf(pFile, "    keyValue [  ");
      for (int cmp = 0; cmp < N; cmp++)
      {
        if (cmp % Freq)
          continue;
        fprintf(pFile, "%e %e %e,\n", dataPlot.getValue(cmp, 1 + 7 * num), dataPlot.getValue(cmp, 2 + 7 * num), dataPlot.getValue(cmp, 3 + 7 * num));
        //fprintf(pFile,"%e %e %e,\n",dataPlot.getValue(0,1+7*num),dataPlot.getValue(0,2+7*num),dataPlot.getValue(0,3+7*num));
      }
      fprintf(pFile, "  ]\n}\n");

    }

    fprintf(pFile, "DEF Animation_Time TimeSensor {\n");
    fprintf(pFile, "     cycleInterval 10.0000\n");
    fprintf(pFile, "     startTime 0\n");
    fprintf(pFile, "     stopTime 0\n");
    fprintf(pFile, "     loop TRUE\n");
    fprintf(pFile, "   }\n ]\n}");
    fclose(pFile);
    system("cat Pantographe1.wrl data.wrl Pantographe2.wrl > run.wrl");
    // --- Output files ---
    //cout<<"====> Output file writing ..."<<endl;
    //ioMatrix io("result.dat", "ascii");
    //io.write(dataPlot,"noDim");

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

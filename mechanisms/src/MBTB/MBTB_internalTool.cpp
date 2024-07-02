/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "MBTB_DATA.hpp"
// #include "MBTB_PYTHON_API.hpp"
// #include "SiconosKernel.hpp"
#include "CADMBTB_API.hpp"
#include "MBTB_Body.hpp"
#include "MBTB_Contact.hpp"
#include "MBTB_JointR.hpp"
#include "MBTB_internalTool.hpp"
#include "ace.h"
// #include "MBTB_TimeSteppingProj.hpp"
// #include "MBTB_TimeSteppingCombinedProj.hpp"
#include "Interaction.hpp"
#include "NewtonEuler1DR.hpp"
#include "NewtonEulerJointR.hpp"
#include "RotationQuaternion.hpp"
// #include "SolverOptions.h"   // for SolverOptions struct
// #include "SiconosAlgebraProd.hpp" // for prod

void siconos::mechanisms::mbtb::internal::MBTB_updateContactFromDS() {
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
#ifdef PRINT_FORCE_CONTACTS
    printf("....> contact force of %s :",
           siconos::mechanisms::mbtb::data::sContacts[numC]->contactName());
    mbtb::data::sContacts[numC]->relation()->contactForce()->display();
#endif
    int index1 = mbtb::data::sContacts[numC]->_indexBody1;
    int index2 = mbtb::data::sContacts[numC]->_indexBody2;

    CADMBTB_moveModelFromModel(mbtb::data::sContacts[numC]->_indexCAD1, index1);
    if (index2 != -1)
      CADMBTB_moveModelFromModel(mbtb::data::sContacts[numC]->_indexCAD2,
                                 index2);

#ifdef DRAW_CONTACT_FORCES

#endif
    CADMBTB_moveGraphicalModelFromModel(mbtb::data::sContacts[numC]->_indexCAD1,
                                        index1);
    if (index2 != -1)
      CADMBTB_moveGraphicalModelFromModel(
          mbtb::data::sContacts[numC]->_indexCAD2, index2);
  }
}

void siconos::mechanisms::mbtb::internal::MBTB_updateContactFromDS(int numDS) {
  for (int numC = 0; numC < (int)mbtb::data::sNbOfContacts; numC++) {
    if ((int)(mbtb::data::sContacts[numC]->_indexBody1) == numDS) {
      CADMBTB_moveModelFromModel(mbtb::data::sContacts[numC]->_indexCAD1,
                                 numDS);
      CADMBTB_moveGraphicalModelFromModel(
          mbtb::data::sContacts[numC]->_indexCAD1, numDS);
    }
    if (mbtb::data::sContacts[numC]->_indexBody2 == numDS) {
      CADMBTB_moveModelFromModel(mbtb::data::sContacts[numC]->_indexCAD2,
                                 numDS);
      CADMBTB_moveGraphicalModelFromModel(
          mbtb::data::sContacts[numC]->_indexCAD2, numDS);
    }
  }
}

void siconos::mechanisms::mbtb::internal::MBTB_DRAW_STEP() {
  /*delete previous*/
  if (mbtb::data::sDrawMode) {
    for (unsigned int nC = 0; nC < 3 * mbtb::data::sNbOfContacts; nC++) {
      CADMBTB_buildLineArtefactLine(nC, nullptr, nullptr, nullptr, nullptr,
                                    nullptr, nullptr);
    }
  }
  if (mbtb::data::sDrawMode & MBTBConst::FaceNormal1) {
    for (unsigned int nC = 0; nC < mbtb::data::sNbOfContacts; nC++) {
      double x1, x2, y1, y2, z1, z2;  //,nx,ny,nz,MinDist;
      // int index1=mbtb::data::sContacts[nC]->_indexBody1;
      // int index2=mbtb::data::sContacts[nC]->_indexBody2;
      // int normalFromFace1;

      x2 = mbtb::data::sContacts[nC]->relation()->pc2()->getValue(0);
      y2 = mbtb::data::sContacts[nC]->relation()->pc2()->getValue(1);
      z2 = mbtb::data::sContacts[nC]->relation()->pc2()->getValue(2);
      x1 = x2 + 1.1 * mbtb::data::sArtefactLength *
                    mbtb::data::sContacts[nC]->relation()->nc()->getValue(0);
      y1 = y2 + 1.1 * mbtb::data::sArtefactLength *
                    mbtb::data::sContacts[nC]->relation()->nc()->getValue(1);
      z1 = z2 + 1.1 * mbtb::data::sArtefactLength *
                    mbtb::data::sContacts[nC]->relation()->nc()->getValue(2);
      /*  CADMBTB_getMinDistance(nC,index1,index2,
          x1,y1,z1,
          x2,y2,z2,
          nx,ny,nz,
          normalFromFace1,
          MinDist);
          printf ("second point of contact : x2=%lf
         ,y2=%lf,z2=%lf\n",x2,y2,z2);*/
      //  CADMBTB_buildOrientedLineArtefactLine(nC+mbtb::data::sNbOfContacts,&x2,&y2,&z2,&x1,&y1,&z1);
    }
  }

  if (mbtb::data::sDrawMode & MBTBConst::ArtefactP1P2) {
    for (unsigned int nC = 0; nC < mbtb::data::sNbOfContacts; nC++) {
      double x1, x2, y1, y2, z1, z2;
      x1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(0);
      y1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(1);
      z1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(2);
      x2 = mbtb::data::sContacts[nC]->relation()->pc2()->getValue(0);
      y2 = mbtb::data::sContacts[nC]->relation()->pc2()->getValue(1);
      z2 = mbtb::data::sContacts[nC]->relation()->pc2()->getValue(2);
      // printf ("second point of contact : x2=%lf ,y2=%lf,z2=%lf\n",x2,y2,z2);
      CADMBTB_buildLineArtefactLine(nC, &x1, &y1, &z1, &x2, &y2, &z2);
    }
  }

  if (mbtb::data::sDrawMode & MBTBConst::ArtefactNormal) {
    for (unsigned int nC = 0; nC < mbtb::data::sNbOfContacts; nC++) {
      double x1, x2, y1, y2, z1, z2;
      x1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(0);
      y1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(1);
      z1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(2);
      x2 = x1 + 1.1 * mbtb::data::sArtefactLength *
                    mbtb::data::sContacts[nC]->relation()->nc()->getValue(0);
      y2 = y1 + 1.1 * mbtb::data::sArtefactLength *
                    mbtb::data::sContacts[nC]->relation()->nc()->getValue(1);
      z2 = z1 + 1.1 * mbtb::data::sArtefactLength *
                    mbtb::data::sContacts[nC]->relation()->nc()->getValue(2);
      CADMBTB_buildOrientedLineArtefactLine(nC + mbtb::data::sNbOfContacts, &x1,
                                            &y1, &z1, &x2, &y2, &z2);
    }
  }

  if (mbtb::data::sDrawMode & MBTBConst::ArtefactReaction) {
    // printf("MBTB_DRAW_STEP REACTION\n");
    SP::Topology topo = sSimu->nonSmoothDynamicalSystem()->topology();
    double h = sSimu->timeStep();
    SP::InteractionsGraph indexSet1 = topo->indexSet(1);

    double FMax = 0;
    double aux, normF;
    int nbUR = 0;
    int nbR = 0;
    InteractionsGraph::VIterator ui1, ui1end;
    boost::tie(ui1, ui1end) = indexSet1->vertices();

    if (mbtb::data::sNominalForce > 1e-12)
      FMax = mbtb::data::sNominalForce;
    else
      for (; ui1 != ui1end; ++ui1) {
        nbUR++;
        std::shared_ptr<siconos::modeling::Interaction> inter1 =
            indexSet1->bundle(*ui1);
        SP::Relation R = inter1->relation();
        for (unsigned int nC = 0; nC < mbtb::data::sNbOfContacts; nC++) {
          if (mbtb::data::sContacts[nC]->relation() == R) {
            nbR++;
            auto F = mbtb::data::sContacts[nC]->relation()->contactForce();
            aux = sqrt(F->getValue(0) * F->getValue(0) +
                       F->getValue(1) * F->getValue(1) +
                       F->getValue(2) * F->getValue(2)) /
                  h;
            if (aux > FMax) FMax = aux;
          }
        }
      }
    //    printf("MBTB_DRAW_STEP REACTION nb Ur=%i, nb
    //    Rcontact=%i.\n",nbUR,nbR);

    // printf("siconos::mechanisms::mbtb::internal::MBTB_DRAW_STEP FMAX=%e
    // \n",FMax);

    boost::tie(ui1, ui1end) = indexSet1->vertices();

    for (; ui1 != ui1end; ++ui1) {
      auto inter1 = indexSet1->bundle(*ui1);
      SP::Relation R = inter1->relation();
      for (unsigned int nC = 0; nC < mbtb::data::sNbOfContacts; nC++) {
        if (mbtb::data::sContacts[nC]->relation() == R) {
          auto F = mbtb::data::sContacts[nC]->relation()->contactForce();
          normF = sqrt(F->getValue(0) * F->getValue(0) +
                       F->getValue(1) * F->getValue(1) +
                       F->getValue(2) * F->getValue(2));
          aux = normF / h;
          aux = aux / FMax;
          if (aux > 1)
            aux = log(aux) + 1.0;
          else
            aux = 1.0 / (-log(aux) + 1.0);
          if (aux > 5) aux = 5;

          // printf("siconos::mechanisms::mbtb::internal::MBTB_DRAW_STEP aux=%e
          // \n",aux);
          if (aux > mbtb::data::sArtefactThreshold) {
            auto F = mbtb::data::sContacts[nC]->relation()->contactForce();
            double x1, x2, y1, y2, z1, z2;
            x1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(0);
            y1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(1);
            z1 = mbtb::data::sContacts[nC]->relation()->pc1()->getValue(2);

            x2 = x1 +
                 aux * (mbtb::data::sArtefactLength / normF) * F->getValue(0);
            y2 = y1 +
                 aux * (mbtb::data::sArtefactLength / normF) * F->getValue(1);
            z2 = z1 +
                 aux * (mbtb::data::sArtefactLength / normF) * F->getValue(2);
            double radius = 0.03 * aux * mbtb::data::sArtefactLength;
            CADMBTB_buildCylinderArtefactLine(
                nC + 2 * mbtb::data::sNbOfContacts, &x1, &y1, &z1, &x2, &y2,
                &z2, &radius);
          }
        }
      }
    }
  }
  CADMBTB_updateGraphic();
}
void siconos::mechanisms::mbtb::internal::MBTB_STEP() {
  MBTB_updateDSFromSiconos();
  siconos::mechanisms::mbtb::internal::MBTB_updateContactFromDS();
  ACE_times[ACE_TIMER_SICONOS].start();
  sSimu->setNewtonTolerance(mbtb::data::sDParams[2]);
  sSimu->setNewtonMaxIteration(mbtb::data::sDParams[3]);
  sSimu->advanceToEvent();

  double* dd = sSimu->oneStepNSProblem(0)->numericsSolverOptions()->dparam;
  int* ii = sSimu->oneStepNSProblem(0)->numericsSolverOptions()->iparam;

  std::cout << "        OSNS solver info :: reached accuracy = "
            << dd[SICONOS_DPARAM_RESIDU] << " < " << dd[SICONOS_DPARAM_TOL];
  std::cout << "     nb iterations =" << ii[SICONOS_IPARAM_ITER_DONE] << " < "
            << ii[SICONOS_IPARAM_MAX_ITER] << std::endl;
  std::cout << "        Newton loop info ::  iterations = "
            << sSimu->getNewtonNbIterations() << " < "
            << sSimu->newtonMaxIteration()
            << " residu = " << sSimu->newtonResiduDSMax() << " < "
            << sSimu->newtonTolerance() << std::endl;

  Type::Siconos simuType;
  simuType = Type::value(*sSimu);
  if (simuType == Type::TimeStepping) {
  } else if (simuType == Type::TimeSteppingDirectProjection) {
  } else if (simuType == Type::TimeSteppingCombinedProjection) {
    std::cout << "     Number of projection iterations = "
              << (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(
                      sSimu))
                     ->nbProjectionIteration()
              << std::endl;
    std::cout << "     Number of cumulated Newton iterations = "
              << (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(
                      sSimu))
                     ->cumulatedNewtonNbIterations()
              << std::endl;
    std::cout << "     Number of set  iterations = "
              << (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(
                      sSimu))
                     ->nbIndexSetsIteration()
              << std::endl;
    std::cout << "     Max violation unilateral = "
              << (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(
                      sSimu))
                     ->maxViolationUnilateral()
              << std::endl;
    std::cout << "     Max violation equality = "
              << (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(
                      sSimu))
                     ->maxViolationEquality()
              << std::endl;
  }

  // sSimu->oneStepNSProblem(0)->display();
  ACE_times[ACE_TIMER_SICONOS].stop();

  if (mbtb::data::sTimerCmp % mbtb::data::sFreqGraphic == 0) {
    siconos::mechanisms::mbtb::internal::MBTB_DRAW_STEP();
  }
  sSimu->nextStep();
}

void siconos::mechanisms::mbtb::internal::MBTB_displayStep() {
  // fprintf(fp,"siconos::mechanisms::mbtb::internal::MBTB_printStep ");
  // printf("%d\t",mbtb::data::sTimerCmp);

  // Bodies display output
  if (mbtb::data::sDisplayStepBodies) {
    printf("[mechanisms] Time step  Number %d\t", mbtb::data::sTimerCmp);
    for (unsigned int numDS = 0; numDS < mbtb::data::sNbOfBodies; numDS++) {
      printf("Body number %i\n", numDS);
      printf("Position of body %i\n", numDS);
      for (unsigned int ii = 0; ii < mbtb::data::sDS[numDS]->q()->size();
           ii++) {
        printf("%e", mbtb::data::sDS[numDS]->q()->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Velocity of body %i\n", numDS);
      for (unsigned int ii = 0; ii < mbtb::data::sDS[numDS]->twist()->size();
           ii++) {
        printf("%e", mbtb::data::sDS[numDS]->twist()->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Kinetic Energy of body %i\n", numDS);
      /*Ec of the DS*/
      //   printf("MBTB Ec computattiom masse matrix:\n");
      //    (mbtb::data::sDS[numDS]->M())->display();
      //    (mbtb::data::sDS[numDS]->twist())->display();
      siconos::algebra::SiconosVector res(6);
      prod(*(mbtb::data::sDS[numDS]->mass()),
           *(mbtb::data::sDS[numDS]->twist()), res);
      double ec = 0.0;
      for (int i = 0; i < 6; i++)
        ec += res.getValue(i) * mbtb::data::sDS[numDS]->twist()->getValue(i);
      printf("%e\t", ec * 0.5);
      printf("\n");
    }
  }

  // Joints display output
  if (mbtb::data::sDisplayStepJoints) {
    printf("STEP Number = %d\t", mbtb::data::sTimerCmp);
    for (unsigned int numJ = 0; numJ < mbtb::data::sNbOfJoints; numJ++) {
      printf("Joint number %i\n", numJ);
      printf("interactionjointR->display  %i\n", numJ);
      mbtb::data::sJointRelations[numJ]->_interaction->display();
      printf("\n");
      printf("Forces in Joint  %i\n", numJ);
      for (int ii = 0; ii < 3; ii++) {
        printf("%e", mbtb::data::sJointRelations[numJ]
                         ->_jointR->contactForce()
                         ->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Moments in Joint  %i\n", numJ);

      auto vaux(new siconos::algebra::SiconosVector(3));
      for (int ii = 3; ii < 6; ii++) {
        vaux->setValue(ii - 3, mbtb::data::sJointRelations[numJ]
                                   ->_jointR->contactForce()
                                   ->getValue(ii));
        printf("%e", vaux->getValue(ii - 3));
        printf("\t");
      }
      /*convert momentum in abs frame*/
      printf("\n");

      printf("Moments in Joint %i in absolute frame \n", numJ);
      siconos::geometry::changeFrameBodyToAbs(
          mbtb::data::sJointRelations[numJ]->_ds1->q(), vaux);
      for (int ii = 0; ii < 3; ii++) {
        printf("%e", vaux->getValue(ii));
        printf("\t");
      }
      printf("\n");
      printf("Equivalent Forces to moments in Joint %i in absolute frame \n",
             numJ);
      if (mbtb::data::sJointRelations[numJ]->_G0C1) {
        mbtb::data::sJointRelations[numJ]->computeEquivalentForces();
        for (int ii = 0; ii < 6; ii++) {
          printf("%e", mbtb::data::sJointRelations[numJ]->_F->getValue(ii));
          printf("\t");
        }
      } else {
        printf("N/A\t");
      }
      printf("\n");
    }
  }
  // Contacts display output
  if (mbtb::data::sDisplayStepContacts) {
    printf("STEP Number = %d\t", mbtb::data::sTimerCmp);
    for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
      printf("Contact number %i\n", numC);
      printf("Contact forces in contact  %i\n", numC);
      for (int ii = 0; ii < 3; ii++) {
        printf(
            "%e",
            mbtb::data::sContacts[numC]->relation()->contactForce()->getValue(
                ii));
        printf("\t");
      }

      auto vaux(new siconos::algebra::SiconosVector(3));
      for (int ii = 3; ii < 6; ii++) {
        vaux->setValue(
            ii - 3,
            mbtb::data::sContacts[numC]->relation()->contactForce()->getValue(
                ii));
      }
      /*convert momentum in abs frame*/
      siconos::geometry::changeFrameBodyToAbs(
          mbtb::data::sDS[mbtb::data::sContacts[numC]->_indexBody1]->q(), vaux);
      printf("\n");
      printf("Moments of contact forces in contact  %i in absolute frame \n",
             numC);
      for (int ii = 0; ii < 3; ii++) {
        printf("%e", vaux->getValue(ii));
        printf("\t");
      }

      printf("\n");

      SP::InteractionsGraph indexSet1 = myNsds->topology()->indexSet(1);
      InteractionsGraph::VIterator ui1, ui1end, v1next;
      boost::tie(ui1, ui1end) = indexSet1->vertices();
      // Remove interactions from the indexSet1
      int find = 0;
      for (v1next = ui1; ui1 != ui1end; ui1 = v1next) {
        ++v1next;
        auto urI = indexSet1->bundle(*ui1);
        if (&(*urI) == &(*(mbtb::data::sInterContacts[numC]))) {
          find = 1;
        }
      }

      printf("Contact status %d\n", find);
      printf("Coordinates of first contact point for contact %i", numC);

      printf("%e\t%e\t%e\t",
             mbtb::data::sContacts[numC]->relation()->pc1()->getValue(0),
             mbtb::data::sContacts[numC]->relation()->pc1()->getValue(1),
             mbtb::data::sContacts[numC]->relation()->pc1()->getValue(2));
      printf("\n");
      if (mbtb::data::sContacts[numC]->relation()->pc2()) {
        printf("Coordinates of second contact point for contact %i", numC);
        printf("%e\t%e\t%e\t",
               mbtb::data::sContacts[numC]->relation()->pc2()->getValue(0),
               mbtb::data::sContacts[numC]->relation()->pc2()->getValue(1),
               mbtb::data::sContacts[numC]->relation()->pc2()->getValue(2));
        printf("\n");
      }
      printf("Gap  for contact %i\t =\t ", numC);
      printf("%e\n",
             mbtb::data::sContacts[numC]->interaction()->y(0)->getValue(0));
      printf("vitess  for contact %i\t =\t ", numC);
      printf("%e\n",
             mbtb::data::sContacts[numC]->interaction()->y(1)->getValue(0));
    }
  }
  // printf("\n");
}
FILE* siconos::mechanisms::mbtb::internal::MBTB_open(std::string filename,
                                                     std::string args) {
  FILE* fp = fopen(filename.c_str(), args.c_str());
  return fp;
}

void siconos::mechanisms::mbtb::internal::MBTB_close(FILE* fp) { fclose(fp); }

void siconos::mechanisms::mbtb::internal::MBTB_printStep(FILE* fp) {
  // fprintf(fp,"siconos::mechanisms::mbtb::internal::MBTB_printStep ");
  fprintf(fp, "%d\t", mbtb::data::sTimerCmp);
  for (unsigned int numDS = 0; numDS < mbtb::data::sNbOfBodies; numDS++) {
    for (unsigned int ii = 0; ii < mbtb::data::sDS[numDS]->q()->size(); ii++) {
      fprintf(fp, "%e", mbtb::data::sDS[numDS]->q()->getValue(ii));
      fprintf(fp, "\t");
    }
    for (unsigned int ii = 0; ii < mbtb::data::sDS[numDS]->twist()->size();
         ii++) {
      fprintf(fp, "%e", mbtb::data::sDS[numDS]->twist()->getValue(ii));
      fprintf(fp, "\t");
    }
    /*Ec of the DS*/
    //   printf("MBTB Ec computattiom masse matrix:\n");
    //    (mbtb::data::sDS[numDS]->M())->display();
    //    (mbtb::data::sDS[numDS]->twist())->display();
    siconos::algebra::SiconosVector res(6);
    prod(*(mbtb::data::sDS[numDS]->mass()), *(mbtb::data::sDS[numDS]->twist()),
         res);
    double ec = 0.0;
    for (int i = 0; i < 6; i++)
      ec += res.getValue(i) * mbtb::data::sDS[numDS]->twist()->getValue(i);
    fprintf(fp, "%e\t", ec * 0.5);
    //    printf("ec=%e",ec*0.5);
  }
  for (unsigned int numJ = 0; numJ < mbtb::data::sNbOfJoints; numJ++) {
    for (unsigned int ii = 0; ii < 3; ii++) {
      fprintf(
          fp, "%e",
          mbtb::data::sJointRelations[numJ]->_jointR->contactForce()->getValue(
              ii));
      fprintf(fp, "\t");
    }
    auto vaux(new siconos::algebra::SiconosVector(3));
    for (unsigned int ii = 3; ii < 6; ii++) {
      vaux->setValue(
          ii - 3,
          mbtb::data::sJointRelations[numJ]->_jointR->contactForce()->getValue(
              ii));
    }
    /*convert momentum in abs frame*/
    siconos::geometry::changeFrameBodyToAbs(
        mbtb::data::sJointRelations[numJ]->_ds1->q(), vaux);
    for (int ii = 0; ii < 3; ii++) {
      fprintf(fp, "%e", vaux->getValue(ii));
      fprintf(fp, "\t");
    }

    if (mbtb::data::sJointRelations[numJ]->_G0C1) {
      mbtb::data::sJointRelations[numJ]->computeEquivalentForces();
      for (int ii = 0; ii < 6; ii++) {
        fprintf(fp, "%e", mbtb::data::sJointRelations[numJ]->_F->getValue(ii));
        fprintf(fp, "\t");
      }
    }
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    for (unsigned int ii = 0; ii < 3; ii++) {
      fprintf(fp, "%e",
              mbtb::data::sContacts[numC]->relation()->contactForce()->getValue(
                  ii));
      fprintf(fp, "\t");
    }
    auto vaux(new siconos::algebra::SiconosVector(3));
    for (unsigned int ii = 3; ii < 6; ii++) {
      vaux->setValue(
          ii - 3,
          mbtb::data::sContacts[numC]->relation()->contactForce()->getValue(
              ii));
    }
    /*convert momentum in abs frame*/
    siconos::geometry::changeFrameBodyToAbs(
        mbtb::data::sDS[mbtb::data::sContacts[numC]->_indexBody1]->q(), vaux);
    for (int ii = 0; ii < 3; ii++) {
      fprintf(fp, "%e", vaux->getValue(ii));
      fprintf(fp, "\t");
    }
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    SP::InteractionsGraph indexSet1 = myNsds->topology()->indexSet(0);
    InteractionsGraph::VIterator ui1, ui1end, v1next;
    boost::tie(ui1, ui1end) = indexSet1->vertices();
    // Remove interactions from the indexSet1
    int find = 0;
    for (v1next = ui1; ui1 != ui1end; ui1 = v1next) {
      ++v1next;
      auto urI = indexSet1->bundle(*ui1);
      if (&(*urI) == &(*(mbtb::data::sInterContacts[numC]))) {
        find = 1;
      }
    }
    fprintf(fp, "%d\t", find);
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    fprintf(fp, "%e\t%e\t%e\t",
            mbtb::data::sContacts[numC]->relation()->pc1()->getValue(0),
            mbtb::data::sContacts[numC]->relation()->pc1()->getValue(1),
            mbtb::data::sContacts[numC]->relation()->pc1()->getValue(2));
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    unsigned int sizeY =
        mbtb::data::sContacts[numC]->interaction()->y(0)->size();
    if (sizeY == 1) {
      fprintf(fp, "%e\t0.\t0.\t",
              mbtb::data::sContacts[numC]->interaction()->y(0)->getValue(0));
    } else {
      fprintf(fp, "%e\t%e\t%e\t",
              mbtb::data::sContacts[numC]->interaction()->y(0)->getValue(0),
              mbtb::data::sContacts[numC]->interaction()->y(0)->getValue(1),
              mbtb::data::sContacts[numC]->interaction()->y(0)->getValue(2));
    }
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    unsigned int sizeY =
        mbtb::data::sContacts[numC]->interaction()->y(0)->size();
    if (sizeY == 1) {
      fprintf(fp, "%e\t0.\t0.\t",
              mbtb::data::sContacts[numC]->interaction()->y(1)->getValue(0));
    } else {
      fprintf(fp, "%e\t%e\t%e\t",
              mbtb::data::sContacts[numC]->interaction()->y(1)->getValue(0),
              mbtb::data::sContacts[numC]->interaction()->y(1)->getValue(1),
              mbtb::data::sContacts[numC]->interaction()->y(1)->getValue(2));
    }
  }

  fprintf(fp, "\n");
}
void siconos::mechanisms::mbtb::internal::MBTB_printHeader(FILE* fp) {
  unsigned int cmp = 1;
  // fprintf(fp,"siconos::mechanisms::mbtb::internal::MBTB_printStep\t");
  fprintf(fp, "stepNum1\t");
  cmp++;
  for (unsigned int numDS = 0; numDS < mbtb::data::sNbOfBodies; numDS++) {
    for (unsigned int icmp = cmp; icmp <= cmp + 7 - 1; ++icmp) {
      fprintf(fp, "position%d_ds_%d\t", icmp, numDS);
    }
    cmp += 7;
    for (unsigned int icmp = cmp; icmp <= cmp + 6 - 1; ++icmp) {
      fprintf(fp, "velocity%d_ds_%d\t", icmp, numDS);
    }
    cmp += 6;
    fprintf(fp, "EC%d_ds_%d\t", cmp, numDS);
    cmp++;
  }
  for (unsigned int numJ = 0; numJ < mbtb::data::sNbOfJoints; numJ++) {
    for (unsigned int icmp = cmp; icmp <= cmp + 6 - 1; ++icmp) {
      fprintf(fp, "jointF%d_%d\t", icmp, numJ);
    }
    cmp += 6;
    if (mbtb::data::sJointRelations[numJ]->_G0C1) {
      for (unsigned int icmp = cmp; icmp <= cmp + 6 - 1; ++icmp) {
        fprintf(fp, "jointEquiF%d_%d\t", icmp, numJ);
      }
      cmp += 6;
    }
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    for (unsigned int icmp = cmp; icmp <= cmp + 6 - 1; ++icmp) {
      fprintf(fp, "ContactForce_%s_%d,%d\t",
              mbtb::data::sContacts[numC]->contactName(), icmp, numC);
    }
    cmp += 6;
  }
  /*  for(int numC=0; numC<mbtb::data::sNbOfContacts; numC++)
  {
    fprintf(fp,"MomentsOfContactForce_%s_%d_%d\t",mbtb::data::sContacts[numC]->contactName(),cmp,cmp+6-1,numC);
    cmp+=6;
    }*/
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    fprintf(fp, "ContactState_%s_%d\t",
            mbtb::data::sContacts[numC]->contactName(), cmp);
    cmp++;
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    for (unsigned int icmp = cmp; icmp <= cmp + 2; ++icmp) {
      fprintf(fp, "ContactPoint_%s_%d\t",
              mbtb::data::sContacts[numC]->contactName(), icmp);
    }
    cmp += 3;
  }
  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    for (unsigned int icmp = cmp; icmp <= cmp + 2; ++icmp) {
      fprintf(fp, "ContactGap_%s_%d\t",
              mbtb::data::sContacts[numC]->contactName(), icmp);
    }
    cmp += 3;
  }

  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    for (unsigned int icmp = cmp; icmp <= cmp + 2; ++icmp) {
      fprintf(fp, "ContactVitess_%s_%d\t",
              mbtb::data::sContacts[numC]->contactName(), icmp);
    }
    cmp += 3;
  }
  fprintf(fp, "\n");
  siconos::mechanisms::mbtb::internal::MBTB_printStep(fp);
}

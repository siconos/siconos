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

#include "MBTB_PYTHON_API.hpp"

#include <cassert>

#include "CADMBTB_API.hpp"  // for CADMBTB_init ...
#include "MBTB_Body.hpp"
#include "MBTB_DATA.hpp"
#include "MBTB_internalTool.hpp"  // For MBTB_updateContactFromDS
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #include "CADMBTB_PYTHON_API.hpp"
#include <boost/math/quaternion.hpp>
// #include "KneeJointR.hpp"
#include "BoundaryCondition.hpp"
#include "EqualityConditionNSL.hpp"
#include "Interaction.hpp"
#include "MBTB_Contact.hpp"
#include "MBTB_ContactRelation.hpp"
#include "MBTB_FC3DContactRelation.hpp"
#include "MBTB_JointR.hpp"
#include "MBTB_TimeStepping.hpp"
#include "MBTB_TimeSteppingCombinedProj.hpp"
#include "MBTB_TimeSteppingProj.hpp"
#include "NonSmoothDynamicalSystem.hpp"
#include "OneStepNSProblem.hpp"
#include "PivotJointR.hpp"
#include "PrismaticJointR.hpp"
#include "TimeStepping.hpp"
#include "ace.h"
// #include <BRepTools.hxx>
#include "GenericMechanical.hpp"
#include "MLCPProjectOnConstraints.hpp"
#include "MoreauJeanCombinedProjectionOSI.hpp"
#include "MoreauJeanDirectProjectionOSI.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "SolverOptions.h"  // for SolverOptions struct
#include "TimeDiscretisation.hpp"
#include "Topology.hpp"
// #define MBTB_MOREAU_YES
//  #define DEBUG_STDOUT
//  #define DEBUG_NOCOLOR
//  #define DEBUG_MESSAGES
#include "siconos_debug.h"
#ifdef MBTB_MOREAU_YES
#include "MBTB_MoreauJeanOSI.hpp"
#endif

#define MBTB_LOAD_CONTACT

void siconos::mechanisms::MBTB_init(unsigned int NumOfBodies, unsigned int NumOfJoints,
                                    unsigned int NumOfContacts) {
  assert(NumOfBodies < MBTB_MAX_BODIES_NUMBER && "MBTB_init NumOfBodies out of range");
  assert(NumOfJoints < MBTB_MAX_JOINTS_NUMBER && "MBTB_init NumOfJoints out of range");
  assert(NumOfContacts < MBTB_MAX_CONTACTS_NUMBER && "MBTB_init NumOfContacts out of range");
  ACE_INIT_TIME();
  mbtb::data::sNbOfBodies = NumOfBodies;
  mbtb::data::sNbOfJoints = NumOfJoints;
  mbtb::data::sNbOfContacts = NumOfContacts;
  CADMBTB_init(mbtb::data::sNbOfBodies + 2 * NumOfContacts, NumOfContacts);
  CADMBTB_setNbOfArtefacts(4 * NumOfContacts); /** P1P2, NORMAL, REACTION */

  // -------------
  // --- Model ---
  // -------------
  mbtb::data::myNsds = std::make_shared<siconos::modeling::NonSmoothDynamicalSystem>(
      mbtb::data::myt0, mbtb::data::myTf);
}

/*get the quaternion from siconos and 1787update the CADs model*/
void siconos::mechanisms::MBTB_updateDSFromSiconos() {
  // ACE_times[ACE_TIMER_UPDATE_POS].start();
  for (unsigned int numDS = 0; numDS < mbtb::data::sNbOfBodies; numDS++) {
    auto q = mbtb::data::sDS[numDS]->q();
    // printf("step %d siconos %s ->q:\n",mTimerCmp,sPieceName[numDS]);
    // siconos::algebra::print(*q);
    double x = (*q)(0);
    double y = (*q)(1);
    double z = (*q)(2);
    double q1 = (*q)(3);
    double q2 = (*q)(4);
    double q3 = (*q)(5);
    double q4 = (*q)(6);
    ACE_times[ACE_TIMER_UPDATE_POS].start();
    CADMBTB_moveObjectFromQ(numDS, x, y, z, q1, q2, q3, q4);
    ACE_times[ACE_TIMER_UPDATE_POS].stop();
    int res = mbtb::data::sTimerCmp % FREQ_UPDATE_GRAPHIC;
    ACE_times[ACE_TIMER_GRAPHIC].start();
    if (!res) {
      /*THIS CODE REBUILD THE GRAPHICAL MODEL
      getContext()->Erase(data::spAISToposDS[numDS]);
      data::spAISToposDS[numDS] = new AIS_Shape( data::sTopoDSPiece[numDS] );
      getContext()->Display( data::spAISToposDS[numDS], false );*/

      // data::spAISToposDS[numDS]->SetTransformation(&(data::sGeomTrsf[numDS]),true,false);//new
      // Geom_Transformation(sTrsfPiece[numDS]),true);

      CADMBTB_moveGraphicalModelFromModel(numDS, numDS);

      // data::spAISToposDS[numDS]->SetTransformation(&(data::sGeomTrsf[numDS])
      //				     ,false,true);
      //       getContext()->Display( data::spAISToposDS[numDS], false );
    }
    ACE_times[ACE_TIMER_GRAPHIC].stop();
  }
}
void siconos::mechanisms::MBTB_BodyLoadCADFile(unsigned int numDS, const std::string& CADFile,
                                               unsigned int withGraphicModel) {
  assert(mbtb::data::sNbOfBodies > numDS && "MBTB_BodyLoadCADFile numDS out of range.");
  /*1) load the CAD model*/
  char* data = (char*)calloc((CADFile.length() + 1), sizeof(char));
  //  memset((void*)data,0,sizeof(*data));
  strcpy(data, CADFile.c_str());
  CADMBTB_loadCADFile(numDS, data);
  free(data);
  if (withGraphicModel) CADMBTB_buildGraphicalModel(numDS);
}

void siconos::mechanisms::MBTB_ContactLoadCADFile(unsigned int contactId,
                                                  const std::string& CADFile1,
                                                  const std::string& CADFile2,
                                                  unsigned int withGraphicModel1,
                                                  unsigned int withGraphicModel2) {
  assert(mbtb::data::sNbOfContacts > contactId &&
         "MBTB_ContactLoadCADFile contactId out of range.");
  char* data = (char*)calloc((CADFile1.length() + 1), sizeof(char));
  unsigned int IdInCAD = mbtb::data::sNbOfBodies + 2 * contactId;

  strcpy(data, CADFile1.c_str());
  CADMBTB_loadCADFile(IdInCAD, data);  // CADFile1.c_str());
  free(data);
  data = (char*)calloc((CADFile2.length() + 1), sizeof(char));

  strcpy(data, CADFile2.c_str());
  CADMBTB_loadCADFile(IdInCAD + 1, data);  // CADFile2.c_str());
  free(data);
  if (withGraphicModel1) CADMBTB_buildGraphicalModel(IdInCAD);
  if (withGraphicModel2) CADMBTB_buildGraphicalModel(IdInCAD + 1);
  CADMBTB_computeUVBounds(IdInCAD);
  CADMBTB_computeUVBounds(IdInCAD + 1);
  CADMBTB_initContact(contactId);
#ifdef MBTB_LOAD_CONTACT
  double U1, U2, V1, V2;
  CADMBTB_getUVBounds(IdInCAD, U1, U2, V1, V2);
  printf("MBTB_LOAD_CONTACT UVBOUNDS idContact1=%d,U1=%e,U2=%e,V1=%e,V2=%e\n", contactId, U1,
         U2, V1, V2);
  CADMBTB_getUVBounds(IdInCAD + 1, U1, U2, V1, V2);
  printf("MBTB_LOAD_CONTACT UVBOUNDS idContact2=%d,U1=%e,U2=%e,V1=%e,V2=%e\n", contactId, U1,
         U2, V1, V2);
#endif
}

namespace siconos::mechanisms::mbtb::internal {  // Local use only
void MBTB_BodyBuildComputeInitPosition(
    unsigned int numDS, double mass, std::shared_ptr<siconos::algebra::SiconosVector> initPos,
    std::shared_ptr<siconos::algebra::SiconosVector> modelCenterMass,
    std::shared_ptr<siconos::algebra::SiconosMatrix> inertialMatrix,
    std::shared_ptr<siconos::algebra::SiconosVector>& q10,
    std::shared_ptr<siconos::algebra::SiconosVector>& v10) {
  assert(mbtb::data::sNbOfBodies > numDS && "MBTB_BodyBuild numDS out of range.");
  /*2)  move the cad model to the initial position*/
  /*It consists in going to the position (x,y,z,q1,q2,q3,q4) starting from
    (0,0,0,1,0,0,0). Endeed, after loading the CAD, the cad model must be moved
    to the initial position of the simulation. This position is not q0 of the
    siconos::DS because siconos work in the frame of G, and G is not necessary
    at the origin.*/
  double q1 = cos(0.5 * (*initPos)(6));
  double q2 = (*initPos)(3) * sin(0.5 * (*initPos)(6));
  double q3 = (*initPos)(4) * sin(0.5 * (*initPos)(6));
  double q4 = (*initPos)(5) * sin(0.5 * (*initPos)(6));
  double x = (*initPos)(0);
  double y = (*initPos)(1);
  double z = (*initPos)(2);

  CADMBTB_moveObjectFromQ(numDS, x, y, z, q1, q2, q3, q4);
  internal::MBTB_updateContactFromDS(numDS);
  /*3) compute the q0 of Siconos, that is the coordinate of G at the initial
   * position*/
  // unsigned int qDim=7;
  // int nDof = 3;
  // unsigned int nDim = 6;
  // auto q10(new SiconosVector(qDim));
  // auto v10(new SiconosVector(nDim));
  q10->setZero();
  v10->setZero();

  /*From the siconos point of view, the dynamic equation are written at the
   * center of gravity.*/
  /*q10 is the coordinate of G in the initial pos:
    --> The initial orientation is still computed.
    --> The translation must be updated because of G.
   */
  ::boost::math::quaternion<double> quattrf(q1, q2, q3, q4);

  ::boost::math::quaternion<double> quatOG(0, (*modelCenterMass)(0),
                                           (*modelCenterMass)(1),
                                           (*modelCenterMass)(2));
  ::boost::math::quaternion<double> quatRes(0, 0, 0, 0);
  quatRes = quattrf * quatOG / quattrf;

  (*q10)(0) = quatRes.R_component_2() + (*initPos)(0);
  (*q10)(1) = quatRes.R_component_3() + (*initPos)(1);
  (*q10)(2) = quatRes.R_component_4() + (*initPos)(2);
  // In current version, the initial orientation is (1,0,0,0)
  (*q10)(3) = q1;
  (*q10)(4) = q2;
  (*q10)(5) = q3;
  (*q10)(6) = q4;
  // siconos::algebra::print(*sq10[numDS]);
  // gp_Ax3 aux=GetPosition(data::sTopoDSPiece[numDS]);
  // printf("and sould be : %e, %e,
  // %e\n",aux.Location().X(),aux.Location().Y(),aux.Location().Z());

  // set the translation of the CAD model.
  double q10x = (*q10)(0);
  double q10y = (*q10)(1);
  double q10z = (*q10)(2);
  CADMBTB_setLocation(numDS, q10x, q10y, q10z);

  // sStartPiece[numDS]=Ax3Aux2;
  CADMBTB_moveGraphicalModelFromModel(numDS, numDS);

  // //In current version I = Id3
  // sI[numDS] = std::make_shared<siconos::algebra::SiconosMatrix>(3,3));
  // sI[numDS]->setZero();
  // //sI[numDS]->setValue(0,0,sMass[numDS]);sI[numDS]->setValue(1,1,sMass[numDS]);sI[numDS]->setValue(2,2,sMass[numDS]);
  // sI[numDS]->setValue(0,0,sMassMatrix[9*numDS+0]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(1,0,sMassMatrix[9*numDS+1]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(2,0,sMassMatrix[9*numDS+2]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(0,1,sMassMatrix[9*numDS+3]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(1,1,sMassMatrix[9*numDS+4]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(2,1,sMassMatrix[9*numDS+5]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(0,2,sMassMatrix[9*numDS+6]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(1,2,sMassMatrix[9*numDS+7]*sMassMatrixScale[numDS]);
  // sI[numDS]->setValue(2,2,sMassMatrix[9*numDS+8]*sMassMatrixScale[numDS]);
  // MBTB_Body * p =new MBTB_Body(q10,v10,mass,inertialMatrix,
  //			       BodyName, CADFile,
  //			       pluginLib, plunginFct);
  // NewtonEulerDS * p1 =new NewtonEulerDS(q10,v10,mass,inertialMatrix);
}

}  // namespace siconos::mechanisms::mbtb::internal

/*Build the MBTB_body and set to the initial postion.*/
void siconos::mechanisms::MBTB_BodyBuild(
    unsigned int numDS, const std::string& BodyName, double mass,
    std::shared_ptr<siconos::algebra::SiconosVector> initPos,
    std::shared_ptr<siconos::algebra::SiconosVector> modelCenterMass,
    std::shared_ptr<siconos::algebra::SiconosMatrix> inertialMatrix,
    const siconos::modeling::func_prototypes::FunctionS_V& boundaryCondition_func,
    const siconos::modeling::BoundaryCondition::Indices& boundaryConditionIndex) {
  assert(mbtb::data::sNbOfBodies > numDS && "MBTB_BodyBuild numDS out of range.");
  unsigned int qDim = 7;
  // int nDof = 3;
  unsigned int nDim = 6;

  auto q10 = std::make_shared<siconos::algebra::SiconosVector>(qDim);
  auto v10 = std::make_shared<siconos::algebra::SiconosVector>(nDim);
  mbtb::internal::MBTB_BodyBuildComputeInitPosition(numDS, mass, initPos, modelCenterMass,
                                                    inertialMatrix, q10, v10);
  MBTB_Body* p = new MBTB_Body(Eigen::Ref<siconos::algebra::SiconosVector>(*q10),
                               Eigen::Ref<siconos::algebra::SiconosVector>(*v10), mass,
                               Eigen::Ref<siconos::algebra::SiconosMatrix>(*inertialMatrix),
                               Eigen::Ref<siconos::algebra::SiconosVector>(*modelCenterMass),
                               BodyName, BodyName);

  // Note FP: review process to set fext, fint ... with user-defined plugins for this MBTB_Body

  // We fix a ds number just to be able to use postprocessing based on hdf5 file
  p->setNumber(numDS + 1);
  // set external forces plugin
  //   p->setComputeFextFunction(fext_func);

  //   p->setComputeMextFunction(mext_func);

  //   p->setComputeFintFunction(fint_func);
  //   if (pluginFintJacqFct == "FiniteDifference") {
  //     std::cout << "setComputeJacobianFintOver_q_byFD(true)" << std::endl;
  //     p->setComputeJacobianFintOver_q_byFD(true);
  //   } else {
  //     p->setComputeJacobianFIntqFunction(pluginFintJacqLib, pluginFintJacqFct);
  //   }
  // }
  // if (pluginFintJacvFct.length() > 1) {
  //   if (pluginFintJacvFct == "FiniteDifference") {
  //     std::cout << "setComputeJacobianFintOver_twist_byFD(true)" << std::endl;
  //     p->setComputeJacobianFintOver_twist_byFD(true);
  //   } else {
  //     p->setComputeJacobianFIntvFunction(pluginFintJacvLib, pluginFintJacvFct);
  //   }
  // }
  // }

  // if (pluginMintFct.length() > 1) {
  //   p->setComputeMIntFunction(pluginMintLib, pluginMintFct);

  //   if (pluginMintJacqFct.length() > 1) {
  //     if (pluginMintJacqFct == "FiniteDifference") {
  //       std::cout << "setComputeJacobianMintOver_q_byFD(true)" << std::endl;
  //       p->setComputeJacobianMintOver_q_byFD(true);
  //     } else {
  //       p->setComputeJacobianMIntqFunction(pluginMintJacqLib, pluginMintJacqFct);
  //     }
  //   }
  //   if (pluginMintJacvFct.length() > 1) {
  //     if (pluginMintJacvFct == "FiniteDifference") {
  //       std::cout << "setComputeJacobianMintOver_twist_byFD(true)" << std::endl;
  //       p->setComputeJacobianMintOver_twist_byFD(true);
  //     } else {
  //       p->setComputeJacobianMIntvFunction(pluginMintJacvLib, pluginMintJacvFct);
  //     }
  //   }
  // }
  // set boundary condition
  DEBUG_PRINT("################################################################\n");

  DEBUG_PRINT("###\n");

  DEBUG_PRINT("###\n");

  DEBUG_PRINT("###\n");

  DEBUG_PRINTF("Set boundary Condition for body numDs = %i\n", numDS);
  DEBUG_EXPR(tools::print("bc indices ", *boundaryConditionIndex));
  auto bd = std::make_shared<siconos::modeling::BoundaryCondition>(boundaryConditionIndex);
  bd->setComputePrescribedVelocityFunction(boundaryCondition_func);
  p->setBoundaryConditions(bd);

  mbtb::data::sDS[numDS].reset(p);
  // sAllDS.insert(mbtb::data::sDS[numDS]);
  //  myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(mbtb::data::sDS[numDS]);
}

void siconos::mechanisms::MBTB_JointBuild(
    unsigned int numJ, const std::string& JointName, JointsType jointType,
    unsigned int indexDS1, unsigned int indexDS2,
    std::shared_ptr<siconos::algebra::SiconosVector> jointPosition) {
  assert(mbtb::data::sNbOfJoints > numJ && "MBTB_JointBuild numJ >=mbtb::data::sNbOfJoints.");
  if (numJ >= mbtb::data::sNbOfJoints) {
    printf("MBTB_JointBuild  numJoint >mbtb::data::sNbOfJoints\n");
    return;
  }
  unsigned int lNbEq = 0;
  unsigned int nbDS = 1;
  unsigned int qDim = 7;
  /*Data to build the graph*/
  mbtb::data::sJointType[numJ] = jointType;
  mbtb::data::sJointIndexDS[2 * numJ] = indexDS1;
  mbtb::data::sJointIndexDS[2 * numJ + 1] = indexDS2;
  /*BUILD H SiconosMatrix and NSLAW*/
  if (jointType == JointsType::Pivot0 || jointType == JointsType::Pivot1) {
    nbDS = 1;
    if (jointType == JointsType::Pivot1) {
      nbDS = 2;
    }
  } else if (jointType == JointsType::Prismatic0) {
    nbDS = 1;
  }

  siconos::algebra::SiconosVector3 P;
  siconos::algebra::SiconosVector3 A;
  auto ds1CenterOfMass = mbtb::data::sDS[indexDS1]->centerOfMass();
  P(0) = (*jointPosition)(3) - (*ds1CenterOfMass)(0);
  P(1) = (*jointPosition)(4) - (*ds1CenterOfMass)(1);
  P(2) = (*jointPosition)(5) - (*ds1CenterOfMass)(2);
  A(0) = (*jointPosition)(0);
  A(1) = (*jointPosition)(1);
  A(2) = (*jointPosition)(2);
  mbtb::data::sJointRelations[numJ] = new MBTB_JointR();
  if (jointType == JointsType::Pivot1) {
    mbtb::data::sJointRelations[numJ]->_jointR =
        std::make_shared<siconos::joints::PivotJointR>(P, A, false, mbtb::data::sDS[indexDS1],
                                                       mbtb::data::sDS[indexDS2]);
    mbtb::data::sJointRelations[numJ]->_ds1 = mbtb::data::sDS[indexDS1];
    // sAllDSByInter[numJ].insert(mbtb::data::sDS[indexDS1]);
    // sAllDSByInter[numJ].insert(mbtb::data::sDS[indexDS2]);
  } else if (jointType == JointsType::Pivot0) {
    mbtb::data::sJointRelations[numJ]->_jointR =
        std::make_shared<siconos::joints::PivotJointR>(P, A, false, mbtb::data::sDS[indexDS1]);
    mbtb::data::sJointRelations[numJ]->_ds1 = mbtb::data::sDS[indexDS1];
    // sAllDSByInter[numJ].insert(mbtb::data::sDS[indexDS1]);
  } else if (jointType == JointsType::Prismatic0) {
    mbtb::data::sJointRelations[numJ]->_jointR =
        std::make_shared<siconos::joints::PrismaticJointR>(A, false,
                                                           mbtb::data::sDS[indexDS1]);
    mbtb::data::sJointRelations[numJ]->_ds1 = mbtb::data::sDS[indexDS1];
    // sAllDSByInter[numJ].insert(mbtb::data::sDS[indexDS1]);
  }

  lNbEq = mbtb::data::sJointRelations[numJ]->_jointR->numberOfConstraints();

  auto lH = std::make_shared<siconos::algebra::SiconosMatrix>(lNbEq, nbDS * qDim);
  lH->setZero();
  auto lNSL = std::make_shared<siconos::modeling::EqualityConditionNSL>(lNbEq);

  mbtb::data::sJointRelations[numJ]->_jointR->setConstantH_NE(*lH);
  //  mbtb::data::sInterJoints[numJ].reset(new Interaction(JointName,
  //  sAllDSByInter[numJ],
  //                                           numJ, lNbEq , lNSL,
  //                                           mbtb::data::sJointRelations[numJ]->_jointR));
  mbtb::data::sInterJoints[numJ] = std::make_shared<siconos::modeling::Interaction>(
      lNSL, mbtb::data::sJointRelations[numJ]->_jointR);
  mbtb::data::sJointRelations[numJ]->_interaction = mbtb::data::sInterJoints[numJ];
  // myModel->nonSmoothDynamicalSystem()->link(mbtb::data::sInterJoints[numJ],
  //                                           mbtb::data::sDS[indexDS1]);
  // if(mbtb::data::sJointType[numJ]==JointsType::Pivot1)
  //   myModel->nonSmoothDynamicalSystem()->link(mbtb::data::sInterJoints[numJ],
  //                                             mbtb::data::sDS[indexDS2]);
}

void siconos::mechanisms::MBTB_ContactBuild(unsigned int numContact,
                                            const std::string& ContactName,
                                            unsigned int indexBody1, int indexBody2,
                                            unsigned int withFriction, double mu, double en,
                                            double et) {
  assert(mbtb::data::sNbOfContacts > numContact &&
         "MBTB_ContactBuild contactId out of range.");
  mbtb::data::sContacts[numContact] = std::make_shared<MBTB_Contact>(
      numContact, ContactName, indexBody1, indexBody2,
      mbtb::data::sNbOfBodies + 2 * numContact, mbtb::data::sNbOfBodies + 2 * numContact + 1,
      withFriction);
  mbtb::data::sContacts[numContact]->set_en(en);

  if (withFriction) {
    auto relation =
        std::make_shared<MBTB_FC3DContactRelation>(mbtb::data::sContacts[numContact]);
    mbtb::data::sContacts[numContact]->relation() = relation;
    mbtb::data::sContacts[numContact]->set_et(et);
    auto nslaw0 = std::make_shared<siconos::modeling::NewtonImpactFrictionNSL>(en, et, mu, 3);
    mbtb::data::sInterContacts[numContact] = std::make_shared<siconos::modeling::Interaction>(
        nslaw0, mbtb::data::sContacts[numContact]->relation());
    // MB : contactName is already in MBTB_Contact!
    // mbtb::data::sInterContacts[numContact]->setId(ContactName);
  } else {
    auto relation = std::make_shared<MBTB_ContactRelation>(mbtb::data::sContacts[numContact]);
    mbtb::data::sContacts[numContact]->relation() = relation;

    auto lNSL = std::make_shared<siconos::modeling::NewtonImpactNSL>(
        mbtb::data::sContacts[numContact]->en());
    mbtb::data::sInterContacts[numContact] = std::make_shared<siconos::modeling::Interaction>(
        lNSL, mbtb::data::sContacts[numContact]->relation());
    //    mbtb::data::sInterContacts[numContact]->setId(ContactName);
  }

  mbtb::data::sContacts[numContact]->setInteraction(mbtb::data::sInterContacts[numContact]);

  // myModel->nonSmoothDynamicalSystem()->link(mbtb::data::sInterContacts[numContact],
  //                                           mbtb::data::sDS[mbtb::data::sContacts[numContact]->_indexBody1]);
  // std::cout << "MBTB_ContactBuild() insert "<<
  // mbtb::data::sContacts[numContact]->_indexBody1 <<std::endl;

  // if(mbtb::data::sContacts[numContact]->_indexBody2!=-1)
  //   myModel->nonSmoothDynamicalSystem()->link(mbtb::data::sInterContacts[numContact],
  //                                             mbtb::data::sDS[mbtb::data::sContacts[numContact]->_indexBody2]);
}
void siconos::mechanisms::MBTB_setSolverIOption(int i, int value) {
  mbtb::data::sSimu->oneStepNSProblem(0)->numericsSolverOptions()->iparam[i] = value;
}
void siconos::mechanisms::MBTB_setSolverDOption(int i, double value) {
  mbtb::data::sSimu->oneStepNSProblem(0)->numericsSolverOptions()->dparam[i] = value;
}
void siconos::mechanisms::MBTB_initSimu(double hTS, int withProj) {
  for (unsigned int numDS = 0; numDS < mbtb::data::sNbOfBodies; numDS++)
    mbtb::data::myNsds->insertDynamicalSystem(mbtb::data::sDS[numDS]);
  for (unsigned int numJ = 0; numJ < mbtb::data::sNbOfJoints; numJ++) {
    if (mbtb::data::sJointType[numJ] == JointsType::Pivot0)
      mbtb::data::myNsds->link(mbtb::data::sInterJoints[numJ],
                               mbtb::data::sDS[mbtb::data::sJointIndexDS[2 * numJ]]);
    if (mbtb::data::sJointType[numJ] == JointsType::Pivot1)
      mbtb::data::myNsds->link(mbtb::data::sInterJoints[numJ],
                               mbtb::data::sDS[mbtb::data::sJointIndexDS[2 * numJ]],
                               mbtb::data::sDS[mbtb::data::sJointIndexDS[2 * numJ + 1]]);
  }

  for (unsigned int numC = 0; numC < mbtb::data::sNbOfContacts; numC++) {
    if (mbtb::data::sContacts[numC]->indexBody2() != -1) {
      DEBUG_PRINT(
          "MBTB_initSimu(double hTS, int withProj). Link contact with two "
          "bodies\n");
      mbtb::data::myNsds->link(mbtb::data::sInterContacts[numC],
                               mbtb::data::sDS[mbtb::data::sContacts[numC]->indexBody1()],
                               mbtb::data::sDS[mbtb::data::sContacts[numC]->indexBody2()]);
      // mbtb::data::sInterContacts[numC]->insert(
      // mbtb::data::sDS[mbtb::data::sContacts[numC]->_indexBody2]  );

    } else {
      DEBUG_PRINT(
          "MBTB_initSimu(double hTS, int withProj). Link contact with one "
          "body\n");
      mbtb::data::myNsds->link(mbtb::data::sInterContacts[numC],
                               mbtb::data::sDS[mbtb::data::sContacts[numC]->indexBody1()]);

      // std::cout <<"link(mbtb::data::sInterContacts[numC],
      // mbtb::data::sDS[mbtb::data::sContacts[numC]->_indexBody1]); " <<
      // std::endl; std::cout <<
      // "============"<<   mbtb::data::sInterContacts[numC] <<std::endl;
    }
  }

  // -- (2) Time discretisation --
  auto t = std::make_shared<siconos::simulation::TimeDiscretisation>(mbtb::data::myt0, hTS);

  // -- (3) one step non smooth problem
  // osnspb.reset(new Equality());
  // osnspb.reset(new MLCP(SICONOS_MLCP_PATH));
  auto osnspb = std::make_shared<siconos::nonsmooth_formulations::GenericMechanical>(
      SICONOS_FRICTION_3D_ONECONTACT_NSN_GP);

  osnspb->setKeepLambdaAndYState(true);
  // osnspb->numericsSolverOptions()->iparam[1]=0;
  osnspb->numericsSolverOptions()->dWork = (double*)malloc(512 * sizeof(double));
  // osnspb->setNumericsVerboseMode(true);

  // osnspb->numericsSolverOptions()->iparam[1]=0;
  // osnspb->numericsSolverOptions()->dparam[0]=1e-5;
  std::shared_ptr<siconos::nonsmooth_formulations::MLCPProjectOnConstraints> osnspb_pos;

  if (withProj) {
    osnspb_pos = std::make_shared<siconos::nonsmooth_formulations::MLCPProjectOnConstraints>(
        SICONOS_MLCP_ENUM);
    // osnspb_pos->setNumericsVerboseMode(1);
  }

  // -- (4) Simulation setup with (1) (2) (3)
#ifdef MBTB_MOREAU_YES
  std::shared_ptr<siconos::mechanisms::MBTB_MoreauJeanOSI> pOSI1;
  std::shared_ptr<siconos::integrators::MoreauJeanCombinedProjectionOSI> pOSI2;
  if (withProj == 0 or withProj == 1) {
    pOSI1 = std::make_shared<siconos::mechanisms::MBTB_MoreauJeanOSI>(mbtb::data::sDParams[0]);
    pOS11->insertDynamicalSystem(mbtb::data::sDS[0]);
    pOSI1->_deactivateYPosThreshold = mbtb::data::sDParams[4];
    pOSI1->_deactivateYVelThreshold = mbtb::data::sDParams[5];
    pOSI1->_activateYPosThreshold = mbtb::data::sDParams[6];
    pOSI1->_activateYVelThreshold = mbtb::data::sDParams[7];
  }
#else
  std::shared_ptr<siconos::integrators::MoreauJeanOSI> pOSI0;
  std::shared_ptr<siconos::integrators::MoreauJeanDirectProjectionOSI> pOSI1;
  std::shared_ptr<siconos::integrators::MoreauJeanCombinedProjectionOSI> pOSI2;
  if (withProj == 0) {
    pOSI0 = std::make_shared<siconos::integrators::MoreauJeanOSI>(mbtb::data::sDParams[0]);
    // pOSI0->insertDynamicalSystem(mbtb::data::sDS[0]);
  } else if (withProj == 1) {
    pOSI1 = std::make_shared<siconos::integrators::MoreauJeanDirectProjectionOSI>(
        mbtb::data::sDParams[0]);
    // pOSI1->insertDynamicalSystem(mbtb::data::sDS[0]);
    pOSI1->setDeactivateYPosThreshold(mbtb::data::sDParams[4]);
    pOSI1->setDeactivateYVelThreshold(mbtb::data::sDParams[5]);
    pOSI1->setActivateYPosThreshold(mbtb::data::sDParams[6]);
    pOSI1->setActivateYVelThreshold(mbtb::data::sDParams[7]);
  }
#endif
  else if (withProj == 2) {
    pOSI2 = std::make_shared<siconos::integrators::MoreauJeanCombinedProjectionOSI>(
        mbtb::data::sDParams[0]);
    // pOSI2->insertDynamicalSystem(mbtb::data::sDS[0]);
  }

  if (withProj == 0) {
    mbtb::data::sSimu = std::make_shared<siconos::mechanisms::MBTB_TimeStepping>(
        mbtb::data::myNsds, t, pOSI0, osnspb);
    auto spSimu = (std::static_pointer_cast<MBTB_TimeStepping>(mbtb::data::sSimu));
  } else if (withProj == 1) {
    mbtb::data::sSimu = std::make_shared<siconos::mechanisms::MBTB_TimeSteppingProj>(
        mbtb::data::myNsds, t, pOSI1, osnspb, osnspb_pos, mbtb::data::sDParams[11]);
    (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))
        ->setProjectionMaxIteration(mbtb::data::sDParams[8]);
    (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))
        ->setConstraintTol(mbtb::data::sDParams[9]);
    (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))
        ->setConstraintTolUnilateral(mbtb::data::sDParams[10]);
  } else if (withProj == 2) {
    mbtb::data::sSimu = std::make_shared<siconos::mechanisms::MBTB_TimeSteppingCombinedProj>(
        mbtb::data::myNsds, t, pOSI2, osnspb, osnspb_pos, 2);
    (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(mbtb::data::sSimu))
        ->setProjectionMaxIteration(mbtb::data::sDParams[8]);
    (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(mbtb::data::sSimu))
        ->setConstraintTol(mbtb::data::sDParams[9]);
    (std::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(mbtb::data::sSimu))
        ->setConstraintTolUnilateral(mbtb::data::sDParams[10]);
  }

  // --  OneStepIntegrators --

  /** \warning VA 3/12/2011
   *  I do not understand why pOSI is multiply reset to another pointer
   *  Is it justified ?
   *  9/06/2015 :  VA commented the use of multiple OSI
   */

  for (unsigned int numDS = 1; numDS < mbtb::data::sNbOfBodies; numDS++) {
#ifdef MBTB_MOREAU_YES
    if (withProj == 0 or withProj == 1) {
      // pOSI1.reset(new MBTB_MoreauJeanOSI(mbtb::data::sDParams[0]));
      // pOSI1->insertDynamicalSystem(mbtb::data::sDS[numDS]);
      //  pOSI1->_deactivateYPosThreshold= mbtb::data::sDParams[4];
      //  pOSI1->_deactivateYVelThreshold= mbtb::data::sDParams[5];
      //  pOSI1->_activateYPosThreshold= mbtb::data::sDParams[6];
      //  pOSI1->_activateYVelThreshold= mbtb::data::sDParams[7];
      //  mbtb::data::sSimu->insertIntegrator(pOSI1);
    }
#else
    if (withProj == 0) {
      // pOSI0.reset(new MoreauJeanOSI(mbtb::data::sDParams[0]));
      // pOSI0->insertDynamicalSystem(mbtb::data::sDS[numDS]);
      //  mbtb::data::sSimu->insertIntegrator(pOSI0);
    } else if (withProj == 1) {
      // pOSI1.reset(new
      // MoreauJeanDirectProjectionOSI(mbtb::data::sDParams[0]));
      // pOSI1->insertDynamicalSystem(mbtb::data::sDS[numDS]);
      //  pOSI1->setDeactivateYPosThreshold(mbtb::data::sDParams[4]);
      //  pOSI1->setDeactivateYVelThreshold(mbtb::data::sDParams[5]);
      //  pOSI1->setActivateYPosThreshold(mbtb::data::sDParams[6];
      //  pOSI1->setActivateYVelThreshold(mbtb::data::sDParams[7]);
      //  mbtb::data::sSimu->insertIntegrator(pOSI1);
    }
#endif
    else if (withProj == 2) {
      // pOSI2.reset(new
      // MoreauJeanCombinedProjectionOSI(mbtb::data::sDParams[0]));
      // pOSI2->insertDynamicalSystem(mbtb::data::sDS[numDS]);
      //  mbtb::data::sSimu->insertIntegrator(pOSI2);
    }
  }

  printf("====> COMPUTE H OF INTERATIONS: (just for display)\n");
  auto indexSet0 = mbtb::data::myNsds->topology()->indexSet0();
  for (unsigned int numJ = 0; numJ < mbtb::data::sNbOfJoints; numJ++) {
    printf("-->compute h of %d \n", numJ);
    auto inter = mbtb::data::sJointRelations[numJ]->_interaction;
    auto& y = *(inter->y(0));
    mbtb::data::sJointRelations[numJ]->_jointR->computeOutput(0., *inter, 0);
    siconos::algebra::print(y);
  }
  printf("====> COMPUTE H OF INTERATION END)\n");

  std::ofstream myfile("simulation_results.dat");
  mbtb::internal::MBTB_printHeader(myfile);
  myfile.close();
  cout << "====> end of initialisation\n\n";
}
std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> siconos::mechanisms::MBTB_nsds() {
  return mbtb::data::myNsds;
}
std::shared_ptr<siconos::simulation::Simulation> siconos::mechanisms::MBTB_simulation() {
  return mbtb::data::sSimu;
}

void siconos::mechanisms::MBTB_doProj(unsigned int v) {
  (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))->setDoProj(v);
}
void siconos::mechanisms::MBTB_doOnlyProj(unsigned int v) {
  (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))->setDoOnlyProj(v);
}
void siconos::mechanisms::MBTB_projectionMaxIteration(unsigned int v) {
  (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))
      ->setProjectionMaxIteration(v);
}
void siconos::mechanisms::MBTB_constraintTol(double v) {
  (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))->setConstraintTol(v);
}

void siconos::mechanisms::MBTB_constraintTolUnilateral(double v) {
  (std::static_pointer_cast<MBTB_TimeSteppingProj>(mbtb::data::sSimu))
      ->setConstraintTolUnilateral(v);
}

void siconos::mechanisms::MBTB_run(int NbSteps) {
  std::ofstream fp("simulation_results.dat", std::ios::app);
  int currentTimerCmp = mbtb::data::sTimerCmp;
  for (int ii = 0; ii < NbSteps; ii++) {
    // while (true){
    mbtb::data::sTimerCmp++;
    if (mbtb::data::sTimerCmp % mbtb::data::sFreqOutput == 0) {
      printf("STEP Number = %d < %d.\n", mbtb::data::sTimerCmp, NbSteps + currentTimerCmp);
    }
    /*NB: first step not useful*/
    mbtb::internal::MBTB_STEP();
    mbtb::internal::MBTB_displayStep();
    if (mbtb::data::sTimerCmp % mbtb::data::sFreqOutput == 0) {
      mbtb::internal::MBTB_printStep(fp);
    }

    // if (mbtb::data::sTimerCmp%mbtb::data::sFreqGraphic==0){
    //    CADMBTB_DumpGraphic();
    // }
    // break;
  }
  fp.close();
  ACE_PRINT_TIME();
  // updateDSFromSiconos();
  // updateContactFromDS();
  // QList<MDIWindow*>::iterator i;
}

void siconos::mechanisms::MBTB_moveBodyToPosWithSpeed(
    unsigned int numDS, std::shared_ptr<siconos::algebra::SiconosVector> aPos,
    std::shared_ptr<siconos::algebra::SiconosVector> aVel) {
  auto q = mbtb::data::sDS[numDS]->q();
  *q = *aPos;
  auto v = mbtb::data::sDS[numDS]->twist();
  *v = *aVel;
  MBTB_updateDSFromSiconos();
  mbtb::internal::MBTB_updateContactFromDS();
  mbtb::data::sDS[numDS]->swapInMemory();
}
void siconos::mechanisms::MBTB_setGraphicFreq(unsigned int freq) {
  mbtb::data::sFreqGraphic = freq;
}
void siconos::mechanisms::MBTB_setOutputFreq(unsigned int freq) {
  mbtb::data::sFreqOutput = freq;
}
void siconos::mechanisms::MBTB_setJointPoints(
    unsigned int numJ, std::shared_ptr<siconos::algebra::SiconosVector> G0C1,
    std::shared_ptr<siconos::algebra::SiconosVector> G0C2) {
  mbtb::data::sJointRelations[numJ]->_G0C1 = G0C1;
  mbtb::data::sJointRelations[numJ]->_G0C2 = G0C2;
}

void siconos::mechanisms::MBTB_ContactSetDParam(unsigned int paramId, unsigned int contactId,
                                                unsigned int idShape, double v) {
  assert(mbtb::data::sNbOfContacts > contactId &&
         "MBTB_ContactLoadCADFile contactId out of range.");
  // unsigned int IdInCAD=mbtb::data::sNbOfBodies+2*contactId;
  switch (paramId) {
    case 1:
      mbtb::data::sContacts[contactId]->set_offset(v);
      break;
    case 2:
      mbtb::data::sArtefactLength = v;
      break;
    case 3:
      mbtb::data::sArtefactThreshold = v;
      break;
    case 4:
      mbtb::data::sNominalForce = v;
      break;
    default:
      printf("Error: MBTB_ContactSetDParam paramId out of range \n");
  }
}
void siconos::mechanisms::MBTB_ContactSetIParam(unsigned int paramId, unsigned int contactId,
                                                unsigned int idShape, bool v) {
  switch (paramId) {
    case 0:
      mbtb::data::sContacts[contactId]->set_offset_to_P1(v);
      break;
    case 1:
      mbtb::data::sContacts[contactId]->set_normal_from_face1(v);
      break;
    case 2:
      mbtb::data::sDrawMode = v;
      break;
    default:
      printf("Error: MBTB_ContactSetIParam paramId out of range \n");
  }
}
void siconos::mechanisms::MBTB_BodySetDParam(unsigned int paramId, unsigned int bodyId,
                                             double v) {
  printf("MBTB_BodySetDParam not yet implemented\n");
}

void siconos::mechanisms::MBTB_BodySetIParam(unsigned int paramId, unsigned int bodyId,
                                             int v) {
  printf("MBTB_BodySetIParam not yet implemented\n");
}
// void siconos::mechanisms::MBTB_BodySetVelocity(
//     unsigned int numDS, std::shared_ptr<siconos::algebra::SiconosVector> aVel) {
//   auto v = mbtb::data::sDS[numDS]->twist();
//   *v = *aVel;
//   auto v0 = mbtb::data::sDS[numDS]->twist0();
//   *v0 = *aVel;
// }
void siconos::mechanisms::MBTB_SetDParam(unsigned int paramId, double v) {
  mbtb::data::sDParams[paramId] = v;
}
void siconos::mechanisms::MBTB_print_dist(unsigned int v) { mbtb::data::sPrintDist = v; }

void siconos::mechanisms::MBTB_displayStep_bodies(unsigned int v) {
  mbtb::data::sDisplayStepBodies = v;
}
void siconos::mechanisms::MBTB_displayStep_joints(unsigned int v) {
  mbtb::data::sDisplayStepJoints = v;
}
void siconos::mechanisms::MBTB_displayStep_contacts(unsigned int v) {
  mbtb::data::sDisplayStepContacts = v;
}

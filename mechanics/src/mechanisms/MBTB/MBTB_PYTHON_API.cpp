#include "MBTB_PYTHON_API.hpp"
#include "MBTB_DATA.hpp"
#include "MBTB_internalTool.hpp"
#include "CADMBTB_API.hpp"
#include "CADMBTB_PYTHON_API.hpp"
#include <boost/math/quaternion.hpp>
#include "KneeJointR.hpp"
#include "PivotJointR.hpp"
#include "PrismaticJointR.hpp"
#include "ace.h"
#include "MBTB_FC3DContactRelation.hpp"
#include "MBTB_TimeStepping.hpp"
#include "MBTB_TimeSteppingProj.hpp"
#include "MBTB_TimeSteppingCombinedProj.hpp"
#include <BRepTools.hxx>

//#define MBTB_MOREAU_YES
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "debug.h"



#ifdef MBTB_MOREAU_YES
#include "MBTB_MoreauJeanOSI.hpp"
#endif

#define MBTB_LOAD_CONTACT

void MBTB_init(unsigned int NumOfBodies, unsigned int NumOfJoints, unsigned int NumOfContacts)
{
  assert(NumOfBodies<MBTB_MAX_BODIES_NUMBER &&"MBTB_init NumOfBodies out of range");
  assert(NumOfJoints<MBTB_MAX_JOINTS_NUMBER &&"MBTB_init NumOfJoints out of range");
  assert(NumOfContacts<MBTB_MAX_CONTACTS_NUMBER &&"MBTB_init NumOfContacts out of range");
  ACE_INIT_TIME();
  sNbOfBodies=NumOfBodies;
  sNbOfJoints=NumOfJoints;
  sNbOfContacts=NumOfContacts;
  //sDS = (SP::MBTB_Body *) malloc(sNbOfBody*sizeof(SP::MBTB_Body));
  //sPieceDraw = (bool*) malloc(sNbOfBody*sizeof(bool));
  CADMBTB_init(sNbOfBodies + 2*NumOfContacts,NumOfContacts);
  CADMBTB_setNbOfArtefacts(4*NumOfContacts); /** P1P2, NORMAL, REACTION */

  myt0 = 0;
  myTf = std::numeric_limits<double>::max();

  // -------------
  // --- Model ---
  // -------------
  myModel.reset(new Model(myt0, myTf));
}

/*get the quaternion from siconos and 1787update the CADs model*/
void MBTB_updateDSFromSiconos()
{
  //ACE_times[ACE_TIMER_UPDATE_POS].start();
  for(unsigned int numDS=0; numDS<sNbOfBodies; numDS++)
  {
    SP::SiconosVector q = sDS[numDS]->q();
    //printf("step %d siconos %s ->q:\n",mTimerCmp,sPieceName[numDS]);
    //q->display();
    double x=q->getValue(0);
    double y=q->getValue(1);
    double z=q->getValue(2);
    double q1=q->getValue(3);
    double q2=q->getValue(4);
    double q3=q->getValue(5);
    double q4=q->getValue(6);
    ACE_times[ACE_TIMER_UPDATE_POS].start();
    CADMBTB_moveObjectFromQ(numDS,x,y,z,q1,q2,q3,q4);
    ACE_times[ACE_TIMER_UPDATE_POS].stop();
    int res = sTimerCmp%FREQ_UPDATE_GRAPHIC;
    ACE_times[ACE_TIMER_GRAPHIC].start();
    if(!res)
    {
      /*THIS CODE REBUILD THE GRAPHICAL MODEL
      getContext()->Erase(spAISToposDS[numDS]);
      spAISToposDS[numDS] = new AIS_Shape( sTopoDSPiece[numDS] );
      getContext()->Display( spAISToposDS[numDS], false );*/

      //spAISToposDS[numDS]->SetTransformation(&(sGeomTrsf[numDS]),true,false);//new Geom_Transformation(sTrsfPiece[numDS]),true);

      CADMBTB_moveGraphicalModelFromModel(numDS,numDS);

      //spAISToposDS[numDS]->SetTransformation(&(sGeomTrsf[numDS])
      //				     ,false,true);
      //      getContext()->Display( spAISToposDS[numDS], false );
    }
    ACE_times[ACE_TIMER_GRAPHIC].stop();
  }

}
void MBTB_BodyLoadCADFile(unsigned int numDS,const std::string& CADFile,unsigned int withGraphicModel)
{
  assert(sNbOfBodies > numDS &&"MBTB_BodyLoadCADFile numDS out of range.");
  /*1) load the CAD model*/
  char * data = (char*)calloc((CADFile.length()+1) , sizeof(char));
  //  memset((void*)data,0,sizeof(*data));
  strcpy(data,CADFile.c_str());
  CADMBTB_loadCADFile(numDS, data);
  free(data);
  if(withGraphicModel)
    CADMBTB_buildGraphicalModel(numDS);
}

void MBTB_ContactLoadCADFile(unsigned int contactId,const std::string& CADFile1,const std::string& CADFile2,unsigned int withGraphicModel1,unsigned int withGraphicModel2)
{
  assert(sNbOfContacts > contactId &&"MBTB_ContactLoadCADFile contactId out of range.");
  char * data = (char*)calloc((CADFile1.length()+1) , sizeof(char));
  unsigned int IdInCAD=sNbOfBodies+2*contactId;

  strcpy(data,CADFile1.c_str());
  CADMBTB_loadCADFile(IdInCAD, data);//CADFile1.c_str());
  free(data);
  data = (char*)calloc((CADFile2.length()+1) , sizeof(char));

  strcpy(data,CADFile2.c_str());
  CADMBTB_loadCADFile(IdInCAD+1, data);//CADFile2.c_str());
  free(data);
  if(withGraphicModel1)
    CADMBTB_buildGraphicalModel(IdInCAD);
  if(withGraphicModel2)
    CADMBTB_buildGraphicalModel(IdInCAD+1);
  CADMBTB_computeUVBounds(IdInCAD);
  CADMBTB_computeUVBounds(IdInCAD+1);
  CADMBTB_initContact(contactId);
#ifdef MBTB_LOAD_CONTACT
  double U1,U2,V1,V2;
  CADMBTB_getUVBounds(IdInCAD,U1,U2,V1,V2);
  printf("MBTB_LOAD_CONTACT UVBOUNDS idContact1=%d,U1=%e,U2=%e,V1=%e,V2=%e\n",contactId,U1,U2,V1,V2);
  CADMBTB_getUVBounds(IdInCAD+1,U1,U2,V1,V2);
  printf("MBTB_LOAD_CONTACT UVBOUNDS idContact2=%d,U1=%e,U2=%e,V1=%e,V2=%e\n",contactId,U1,U2,V1,V2);
#endif
}
void _MBTB_BodyBuildComputeInitPosition(unsigned int numDS,   double mass,
                                        SP::SiconosVector initPos, SP::SiconosVector modelCenterMass,SP::SimpleMatrix inertialMatrix, SP::SiconosVector& q10,SP::SiconosVector& v10)
{
  assert(sNbOfBodies > numDS &&"MBTB_BodyBuild numDS out of range.");
  /*2)  move the cad model to the initial position*/
  /*It consists in going to the position (x,y,z,q1,q2,q3,q4) starting from (0,0,0,1,0,0,0).
    Endeed, after loading the CAD, the cad model must be moved to the initial position of the simulation.
    This position is not q0 of the siconos::DS because siconos work in the frame of G, and G is not necessary at the origin.*/
  double q1=cos(0.5*initPos->getValue(6));
  double q2=initPos->getValue(3)*sin(0.5*initPos->getValue(6));
  double q3=initPos->getValue(4)*sin(0.5*initPos->getValue(6));
  double q4=initPos->getValue(5)*sin(0.5*initPos->getValue(6));
  double x=initPos->getValue(0);
  double y=initPos->getValue(1);
  double z=initPos->getValue(2);

  CADMBTB_moveObjectFromQ(numDS,
                          x,
                          y,
                          z,
                          q1,
                          q2,
                          q3,
                          q4);
  _MBTB_updateContactFromDS(numDS);
  /*3) compute the q0 of Siconos, that is the coordinate of G at the initial position*/
  //unsigned int qDim=7;
  //unsigned int nDof = 3;
  //unsigned int nDim = 6;
  //SP::SiconosVector q10(new SiconosVector(qDim));
  //SP::SiconosVector v10(new SiconosVector(nDim));
  q10->zero();
  v10->zero();


  /*From the siconos point of view, the dynamic equation are written at the center of gravity.*/
  /*q10 is the coordinate of G in the initial pos:
    --> The initial orientation is still computed.
    --> The translation must be updated because of G.
   */
  ::boost::math::quaternion<double>    quattrf(q1,q2,q3,q4);

  ::boost::math::quaternion<double>    quatOG(0,
      modelCenterMass->getValue(0),
      modelCenterMass->getValue(1),
      modelCenterMass->getValue(2));
  ::boost::math::quaternion<double>    quatRes(0,0,0,0);
  quatRes=quattrf*quatOG/quattrf;

  q10->setValue(0,quatRes.R_component_2()+initPos->getValue(0));
  q10->setValue(1,quatRes.R_component_3()+initPos->getValue(1));
  q10->setValue(2,quatRes.R_component_4()+initPos->getValue(2));
  //In current version, the initial orientation is (1,0,0,0)
  q10->setValue(3,q1);
  q10->setValue(4,q2);
  q10->setValue(5,q3);
  q10->setValue(6,q4);
  //sq10[numDS]->display();
  //gp_Ax3 aux=GetPosition(sTopoDSPiece[numDS]);
  //printf("and sould be : %e, %e, %e\n",aux.Location().X(),aux.Location().Y(),aux.Location().Z());

  //set the translation of the CAD model.
  double q10x=q10->getValue(0);
  double q10y=q10->getValue(1);
  double q10z=q10->getValue(2);
  CADMBTB_setLocation(numDS,q10x,q10y,q10z);

  // sStartPiece[numDS]=Ax3Aux2;
  CADMBTB_moveGraphicalModelFromModel(numDS,numDS);

  // //In current version I = Id3
  // sI[numDS].reset(new SimpleMatrix(3,3));
  // sI[numDS]->zero();
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
/*Build the MBTB_body and set to the initial postion.*/
void MBTB_BodyBuild(unsigned int numDS, const std::string& BodyName,  double mass,
                    SP::SiconosVector initPos, SP::SiconosVector modelCenterMass,
                    SP::SimpleMatrix inertialMatrix,
                    const std::string& pluginFextLib,  const std::string& pluginFextFct,
                    const std::string& pluginMextLib,  const std::string& pluginMextFct,
                    const std::string& pluginFintLib,  const std::string& pluginFintFct,
                    const std::string& pluginMintLib,  const std::string& pluginMintFct,
                    const std::string& pluginFintJacqLib,  const std::string& pluginFintJacqFct,
                    const std::string& pluginMintJacqLib,  const std::string& pluginMintJacqFct,
                    const std::string& pluginFintJacvLib,  const std::string& pluginFintJacvFct,
                    const std::string& pluginMintJacvLib,  const std::string& pluginMintJacvFct,
                    const std::string& pluginBoundaryConditionLib,  const std::string& pluginBoundaryConditionFct,
                    SP::IndexInt boundaryConditionIndex)
{
  assert(sNbOfBodies > numDS &&"MBTB_BodyBuild numDS out of range.");
  unsigned int qDim=7;
  //unsigned int nDof = 3;
  unsigned int nDim = 6;

  SP::SiconosVector q10(new SiconosVector(qDim));
  SP::SiconosVector v10(new SiconosVector(nDim));
  _MBTB_BodyBuildComputeInitPosition(numDS,mass,initPos,modelCenterMass,inertialMatrix,q10,v10);
  MBTB_Body * p=0;
  p =new MBTB_Body(q10,v10,mass,inertialMatrix,modelCenterMass,
                   BodyName, BodyName);

  // set external forces plugin
  if(pluginFextFct.length()>1)
  {
    p->setComputeFExtFunction(pluginFextLib,pluginFextFct);
  }
  if(pluginMextFct.length()>1)
  {
    p->setComputeMExtFunction(pluginMextLib,pluginMextFct);
  }
  // set internal forces plugin
  if(pluginFintFct.length()>1)
  {
    p->setComputeFIntFunction(pluginFintLib,pluginFintFct);

    if(pluginFintJacqFct.length()>1)
    {
      if (pluginFintJacqFct == "FiniteDifference")
      {
        std::cout <<"setComputeJacobianFIntqByFD(true)" <<std::endl;
        p->setComputeJacobianFIntqByFD(true);
      }
      else
      {
        p->setComputeJacobianFIntqFunction(pluginFintJacqLib,pluginFintJacqFct);
      }
    }
    if(pluginFintJacvFct.length()>1)
    {
      if (pluginFintJacvFct == "FiniteDifference")
      {
        std::cout <<"setComputeJacobianFIntvByFD(true)" <<std::endl;
        p->setComputeJacobianFIntvByFD(true);
      }
      else
      {
        p->setComputeJacobianFIntvFunction(pluginFintJacvLib,pluginFintJacvFct);
      }
    }
  }

  if(pluginMintFct.length()>1)
  {
    p->setComputeMIntFunction(pluginMintLib,pluginMintFct);

    if(pluginMintJacqFct.length()>1)
    {
      if (pluginMintJacqFct == "FiniteDifference")
      {
        std::cout <<"setComputeJacobianMIntqByFD(true)" <<std::endl;
        p->setComputeJacobianMIntqByFD(true);
      }
      else
      {
        p->setComputeJacobianMIntqFunction(pluginMintJacqLib,pluginMintJacqFct);
      }
    }
    if(pluginMintJacvFct.length()>1)
    {
      if (pluginMintJacvFct == "FiniteDifference")
      {
        std::cout <<"setComputeJacobianMIntvByFD(true)" <<std::endl;
        p->setComputeJacobianMIntvByFD(true);
      }
      else
      {
        p->setComputeJacobianMIntvFunction(pluginMintJacvLib,pluginMintJacvFct);
      }
    }
  }
  // set boundary condition
  if (pluginBoundaryConditionFct.length() >1)
  {
    //SP::IndexInt bdindex(new IndexInt(1));
    //(*bdindex)[0] = 4;
    DEBUG_PRINT("################################################################\n");

    DEBUG_PRINT("###\n");

    DEBUG_PRINT("###\n");

    DEBUG_PRINT("###\n");

    DEBUG_PRINTF("Set boundary Condition for body numDs = %i\n", numDS);
    DEBUG_EXPR(
      for (std::vector<unsigned int>::iterator  itindex = boundaryConditionIndex->begin() ;
           itindex != boundaryConditionIndex->end();
           ++itindex)
      {std::cout << *itindex <<std::endl;};
          );

    SP::BoundaryCondition bd(new BoundaryCondition(boundaryConditionIndex));
    bd->setComputePrescribedVelocityFunction(pluginBoundaryConditionLib, pluginBoundaryConditionFct);
    p->setBoundaryConditions(bd);
  }


  sDS[numDS].reset(p);
  sAllDS.insert(sDS[numDS]);
  // std::cout << "MBTB_BodyBuild()" <<std::endl;
  // sDS[numDS]->display();
  // myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(sDS[numDS]);
}

void MBTB_JointBuild(unsigned int numJ,const std::string& JointName,
                     unsigned int jointType,
                     unsigned int indexDS1, unsigned int indexDS2,
                     SP::SiconosVector jointPosition)
{
  assert(sNbOfJoints > numJ &&"MBTB_JointBuild numJ >=sNbOfJoints.");
  if(numJ >= sNbOfJoints)
  {
    printf("MBTB_JointBuild  numJoint >sNbOfJoints\n");
    return;
  }
  unsigned int lNbEq=0;
  unsigned int nbDS=1;
  unsigned int qDim=7;
  /*Data to build the graph*/
  sJointType[numJ]=jointType;
  sJointIndexDS[2*numJ]=indexDS1;
  sJointIndexDS[2*numJ+1]=indexDS2;
  /*BUILD H SimpleMatrix and NSLAW*/
  if(jointType == PIVOT_0 || jointType == PIVOT_1)
  {
    lNbEq = PivotJointR::numberOfConstraints();
    nbDS=1;
    if(jointType == PIVOT_1)
    {
      nbDS=2;
    }
  }
  else if(jointType == PRISMATIC_0)
  {
    nbDS=1;
    lNbEq =PrismaticJointR::numberOfConstraints();
  }
  SP::SimpleMatrix lH(new SimpleMatrix(lNbEq ,nbDS*qDim));
  lH->zero();
  SP::EqualityConditionNSL lNSL(new EqualityConditionNSL(lNbEq));


  SP::SiconosVector P(new SiconosVector(3));
  SP::SiconosVector A(new SiconosVector(3));
  SP::SiconosVector ds1CenterOfMass = sDS[indexDS1]->centerOfMass();
  //DynamicalSystemsSet lallDS;
  P->setValue(0,jointPosition->getValue(3)-ds1CenterOfMass->getValue(0));
  P->setValue(1,jointPosition->getValue(4)-ds1CenterOfMass->getValue(1));
  P->setValue(2,jointPosition->getValue(5)-ds1CenterOfMass->getValue(2));
  A->setValue(0,jointPosition->getValue(0));
  A->setValue(1,jointPosition->getValue(1));
  A->setValue(2,jointPosition->getValue(2));
  sJointRelations[numJ]= new MBTB_JointR();
  if(jointType == PIVOT_1)
  {
    sJointRelations[numJ]->_jointR.reset(new PivotJointR(sDS[ indexDS1],sDS[indexDS2],P,A));
    sJointRelations[numJ]->_ds1 = sDS[indexDS1];
    sAllDSByInter[numJ].insert(sDS[indexDS1]);
    sAllDSByInter[numJ].insert(sDS[indexDS2]);
  }
  else if(jointType == PIVOT_0)
  {
    sJointRelations[numJ]->_jointR.reset(new PivotJointR(sDS[indexDS1],P,A,false));
    sJointRelations[numJ]->_ds1 = sDS[indexDS1];
    sAllDSByInter[numJ].insert(sDS[indexDS1]);
  }
  else if(jointType == PRISMATIC_0)
  {
    sJointRelations[numJ]->_jointR.reset(new PrismaticJointR(sDS[indexDS1],A));
    sJointRelations[numJ]->_ds1 = sDS[indexDS1];
    sAllDSByInter[numJ].insert(sDS[indexDS1]);
  }
  sJointRelations[numJ]->_jointR->setJachq(lH);
//  sInterJoints[numJ].reset(new Interaction(JointName, sAllDSByInter[numJ],
//                                           numJ, lNbEq , lNSL,
//                                           sJointRelations[numJ]->_jointR));
  sInterJoints[numJ].reset(new Interaction(lNbEq , lNSL,
                                           sJointRelations[numJ]->_jointR));
  sJointRelations[numJ]->_interaction = sInterJoints[numJ];
  // myModel->nonSmoothDynamicalSystem()->link(sInterJoints[numJ],
  //                                           sDS[indexDS1]);
  // if(sJointType[numJ]==PIVOT_1)
  //   myModel->nonSmoothDynamicalSystem()->link(sInterJoints[numJ],
  //                                             sDS[indexDS2]);


}

void MBTB_ContactBuild(unsigned int numContact, const std::string& ContactName,
                       unsigned int indexBody1, int indexBody2,
                       unsigned int withFriction, double mu, double en, double et)
{
  assert(sNbOfContacts > numContact &&"MBTB_ContactBuild contactId out of range.");
  sContacts[numContact]=new MBTB_Contact(numContact,ContactName,
                                         indexBody1, indexBody2,
                                         sNbOfBodies+2*numContact,
                                         sNbOfBodies+2*numContact+1,
                                         withFriction);
  sContacts[numContact]->_en=en;

  if(withFriction)
  {
    sContacts[numContact]->_et=et;
    SP::NonSmoothLaw nslaw0(new NewtonImpactFrictionNSL(en,et,mu,3));
    sInterContacts[numContact].reset(new Interaction(3,nslaw0,sContacts[numContact]->relation(),numContact));
    // MB : contactName is already in MBTB_Contact!
    // sInterContacts[numContact]->setId(ContactName);
  }
  else
  {
    SP::NewtonImpactNSL lNSL(new NewtonImpactNSL(sContacts[numContact]->_en));
    sInterContacts[numContact].reset(new Interaction(1,lNSL,
                                                     sContacts[numContact]->relation(),numContact));
//    sInterContacts[numContact]->setId(ContactName);
  }

  sContacts[numContact]->setInteraction(sInterContacts[numContact]);

  // myModel->nonSmoothDynamicalSystem()->link(sInterContacts[numContact],
  //                                           sDS[sContacts[numContact]->_indexBody1]);
  // std::cout << "MBTB_ContactBuild() insert "<< sContacts[numContact]->_indexBody1 <<std::endl;

  // sDS[sContacts[numContact]->_indexBody1]->display();
  // if(sContacts[numContact]->_indexBody2!=-1)
  //   myModel->nonSmoothDynamicalSystem()->link(sInterContacts[numContact],
  //                                             sDS[sContacts[numContact]->_indexBody2]);


}
void MBTB_setSolverIOption(int i,int value)
{
  myModel->simulation()->oneStepNSProblem(0)->numericsSolverOptions()->iparam[i]=value;
}
void MBTB_setSolverDOption(int i,double value)
{
  myModel->simulation()->oneStepNSProblem(0)->numericsSolverOptions()->dparam[i]=value;
}
void  MBTB_initSimu(double hTS, int withProj)
{

  for(unsigned int numDS =0; numDS<sNbOfBodies; numDS++)
    myModel->nonSmoothDynamicalSystem()->insertDynamicalSystem(sDS[numDS]);
  for(unsigned int numJ=0; numJ<sNbOfJoints; numJ++)
  {
    if (sJointType[numJ]==PIVOT_0)
      myModel->nonSmoothDynamicalSystem()->link(sInterJoints[numJ],
                                                sDS[sJointIndexDS[2*numJ]]);
    if (sJointType[numJ]==PIVOT_1)
      myModel->nonSmoothDynamicalSystem()->link(sInterJoints[numJ],
                                                sDS[sJointIndexDS[2*numJ]],
                                                sDS[sJointIndexDS[2*numJ+1]]);
  }

  for(unsigned int numC=0; numC<sNbOfContacts; numC++)
  {

    if(sContacts[numC]->_indexBody2!=-1)
    {
      DEBUG_PRINT("MBTB_initSimu(double hTS, int withProj). Link contact with two bodies\n");
      myModel->nonSmoothDynamicalSystem()->link(sInterContacts[numC],
                                                sDS[sContacts[numC]->_indexBody1],
                                                sDS[sContacts[numC]->_indexBody2]);
      // sInterContacts[numC]->insert(   sDS[sContacts[numC]->_indexBody2]  );

    }
    else
    {
      DEBUG_PRINT("MBTB_initSimu(double hTS, int withProj). Link contact with one body\n");
      myModel->nonSmoothDynamicalSystem()->link(sInterContacts[numC],
                                                sDS[sContacts[numC]->_indexBody1]);


      // std::cout <<"link(sInterContacts[numC],       sDS[sContacts[numC]->_indexBody1]); " << std::endl;
      // std::cout <<  "============"<<   sInterContacts[numC] <<std::endl;
      // sInterContacts[numC]->display();
      // std::cout << sDS[sContacts[numC]->_indexBody1] << std::endl;
      // sDS[sContacts[numC]->_indexBody1]->display();
      // sInterContacts[numC]->insert(   sDS[sContacts[numC]->_indexBody1]  );
      // sInterContacts[numC]->display();
      //    sInterContacts[numC]->dynamicalSystem(0)->display();
    }
  }



  // -- (2) Time discretisation --
  SP::TimeDiscretisation t(new TimeDiscretisation(myt0,hTS));

  // -- (3) one step non smooth problem
  //osnspb.reset(new Equality());
  //osnspb.reset(new MLCP(SICONOS_MLCP_PATH));
  SP::LinearOSNS osnspb(new GenericMechanical(SICONOS_FRICTION_3D_ONECONTACT_NSN_AC_GP));
  //SP::LinearOSNS osnspb(new GenericMechanical());
  osnspb->setKeepLambdaAndYState(true);
  //osnspb->numericsSolverOptions()->iparam[1]=0;
  osnspb->numericsSolverOptions()->dWork=(double*) malloc(512*sizeof(double));
  //osnspb->setNumericsVerboseMode(true);

  //osnspb->numericsSolverOptions()->iparam[1]=0;
  //osnspb->numericsSolverOptions()->dparam[0]=1e-5;
  SP::MLCPProjectOnConstraints osnspb_pos;

  if(withProj)
  {
    osnspb_pos.reset(new MLCPProjectOnConstraints(SICONOS_MLCP_ENUM));
    //osnspb_pos->setNumericsVerboseMode(1);
  }




  // -- (4) Simulation setup with (1) (2) (3)
#ifdef MBTB_MOREAU_YES
  SP::MBTB_MoreauJeanOSI pOSI1;
  SP::MoreauJeanCombinedProjectionOSI pOSI2;
  if (withProj==0 or withProj==1)
  {
    pOSI1.reset(new MBTB_MoreauJeanOSI(sDParams[0]));
    pOS11->insertDynamicalSystem(sDS[0]);
    pOSI1->_deactivateYPosThreshold= sDParams[4];
    pOSI1->_deactivateYVelThreshold= sDParams[5];
    pOSI1->_activateYPosThreshold= sDParams[6];
    pOSI1->_activateYVelThreshold= sDParams[7];
  }
#else
  SP::MoreauJeanOSI pOSI0;
  SP::MoreauJeanDirectProjectionOSI pOSI1;
  SP::MoreauJeanCombinedProjectionOSI pOSI2;
  if (withProj==0)
  {
    pOSI0.reset(new MoreauJeanOSI(sDParams[0]));
    //pOSI0->insertDynamicalSystem(sDS[0]);
  }
  else if(withProj==1)
  {
    pOSI1.reset(new MoreauJeanDirectProjectionOSI(sDParams[0]));
    //pOSI1->insertDynamicalSystem(sDS[0]);
    pOSI1->setDeactivateYPosThreshold(sDParams[4]);
    pOSI1->setDeactivateYVelThreshold(sDParams[5]);
    pOSI1->setActivateYPosThreshold(sDParams[6]);
    pOSI1->setActivateYVelThreshold(sDParams[7]);
  }
#endif
  else if  (withProj==2)
  {
    pOSI2.reset(new MoreauJeanCombinedProjectionOSI(sDParams[0]));
    //pOSI2->insertDynamicalSystem(sDS[0]);
  }


  if(withProj==0)
  {
    sSimu.reset(new MBTB_TimeStepping(t,pOSI0,osnspb));
    SP::MBTB_TimeStepping spSimu = (std11::static_pointer_cast<MBTB_TimeStepping>(sSimu));
  }
  else if (withProj==1)
  {
    sSimu.reset(new MBTB_TimeSteppingProj(t,pOSI1,osnspb,osnspb_pos,sDParams[11]));
    (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setProjectionMaxIteration(sDParams[8]);
    (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setConstraintTol(sDParams[9]);
    (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setConstraintTolUnilateral(sDParams[10]);
  }
  else if (withProj==2)
  {
    sSimu.reset(new MBTB_TimeSteppingCombinedProj(t,pOSI2,osnspb,osnspb_pos,2));
    (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->setProjectionMaxIteration(sDParams[8]);
    (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->setConstraintTol(sDParams[9]);
    (std11::static_pointer_cast<MBTB_TimeSteppingCombinedProj>(sSimu))->setConstraintTolUnilateral(sDParams[10]);
  }

  // --  OneStepIntegrators --

  /** \warning VA 3/12/2011
   *  I do not understand why pOSI is multiply reset to another pointer
   *  Is it justified ?
   *  9/06/2015 :  VA commented the use of multiple OSI
   */

  for(unsigned int numDS =1; numDS<sNbOfBodies; numDS++)
  {
#ifdef MBTB_MOREAU_YES
    if (withProj==0 or withProj==1)
     {
       //pOSI1.reset(new MBTB_MoreauJeanOSI(sDParams[0]));
       //pOSI1->insertDynamicalSystem(sDS[numDS]);
       // pOSI1->_deactivateYPosThreshold= sDParams[4];
       // pOSI1->_deactivateYVelThreshold= sDParams[5];
       // pOSI1->_activateYPosThreshold= sDParams[6];
       // pOSI1->_activateYVelThreshold= sDParams[7];
       // sSimu->insertIntegrator(pOSI1);
     }
#else
    if (withProj==0)
    {
      //pOSI0.reset(new MoreauJeanOSI(sDParams[0]));
      //pOSI0->insertDynamicalSystem(sDS[numDS]);
      // sSimu->insertIntegrator(pOSI0);
    }
    else if (withProj==1)
    {
      //pOSI1.reset(new MoreauJeanDirectProjectionOSI(sDParams[0]));
      //pOSI1->insertDynamicalSystem(sDS[numDS]);
      // pOSI1->setDeactivateYPosThreshold(sDParams[4]);
      // pOSI1->setDeactivateYVelThreshold(sDParams[5]);
      // pOSI1->setActivateYPosThreshold(sDParams[6];
      // pOSI1->setActivateYVelThreshold(sDParams[7]);
      // sSimu->insertIntegrator(pOSI1);
    }
#endif
    else if  (withProj==2)
    {
      //pOSI2.reset(new MoreauJeanCombinedProjectionOSI(sDParams[0]));
      //pOSI2->insertDynamicalSystem(sDS[numDS]);
      // sSimu->insertIntegrator(pOSI2);
    }
  }

  // --- Simulation initialization ---
  cout <<"\n====> Initialisation ..." <<endl<<endl;
  myModel->setSimulation(sSimu);
  myModel->initialize();

  printf("====> COMPUTE H OF INTERATIONS: (just for display)\n");
  SP::InteractionsGraph indexSet0 = myModel->nonSmoothDynamicalSystem()->topology()->indexSet0();
  for(unsigned int numJ=0; numJ<sNbOfJoints; numJ++)
  {
    printf("-->compute h of %d \n",numJ);
    SP::Interaction inter = sJointRelations[numJ]->_interaction;
    InteractionsGraph::VDescriptor ui = indexSet0->descriptor(inter);
    SiconosVector& y = *(inter->y(0));
    VectorOfBlockVectors& DSlink = *(indexSet0->properties(ui)).DSlink;

    sJointRelations[numJ]->_jointR->computeh(0., *DSlink[NewtonEulerR::q0], y);
  }
  printf("====> COMPUTE H OF INTERATION END)\n");

  FILE *fp;
  fp = fopen("simulation_results.dat", "w");
  _MBTB_printHeader(fp);
  fclose(fp) ;
  NumericsOptions global_options;
  global_options.verboseMode=0;
  setNumericsOptions(&global_options);
  cout <<"====> end of initialisation" <<endl<<endl;
}
SP::Model MBTB_model()
{
  return myModel;
}

void MBTB_doProj(unsigned int v)
{
  (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setDoProj(v);
}
void MBTB_doOnlyProj(unsigned int v)
{
  (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setDoOnlyProj(v);
}
void MBTB_projectionMaxIteration(unsigned int v)
{
  (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setProjectionMaxIteration(v);
}
void MBTB_constraintTol(double v)
{
  (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setConstraintTol(v);
}

void MBTB_constraintTolUnilateral(double v)
{
  (std11::static_pointer_cast<MBTB_TimeSteppingProj>(sSimu))->setConstraintTolUnilateral(v);
}

void MBTB_run(int NbSteps)
{

  FILE *fp;
  fp = fopen("simulation_results.dat", "a");
  int currentTimerCmp = sTimerCmp;
  for(int ii=0; ii<NbSteps; ii++)
  {
    //while (true){
    sTimerCmp++;
 if(sTimerCmp%sFreqOutput==0)
  {
    printf("STEP Number = %d < %d.\n",sTimerCmp,NbSteps+currentTimerCmp);
  }
  /*NB: first step not useful*/
    _MBTB_STEP();
    _MBTB_displayStep();
  if(sTimerCmp%sFreqOutput==0)
  {
    _MBTB_printStep(fp);
  }



    //if (sTimerCmp%sFreqGraphic==0){
    //   CADMBTB_DumpGraphic();
    //}
    //break;
  }
  fclose(fp) ;
  ACE_PRINT_TIME();
  //updateDSFromSiconos();
  //updateContactFromDS();
  //QList<MDIWindow*>::iterator i;
}

void MBTB_moveBodyToPosWithSpeed(unsigned int numDS, SP::SiconosVector aPos, SP::SiconosVector aVel)
{
  SP::SiconosVector q = sDS[numDS]->q();
  *q= *aPos;
  SP::SiconosVector v = sDS[numDS]->velocity();
  *v= *aVel;
  MBTB_updateDSFromSiconos();
  _MBTB_updateContactFromDS();
  sDS[numDS]->swapInMemory();
}
void MBTB_setGraphicFreq(unsigned int freq)
{
  sFreqGraphic=freq;
}
void MBTB_setOutputFreq(unsigned int freq)
{
  sFreqOutput=freq;
}
void MBTB_setJointPoints(unsigned int numJ, SP::SiconosVector G0C1,SP::SiconosVector G0C2)
{
  sJointRelations[numJ]->_G0C1=G0C1;
  sJointRelations[numJ]->_G0C2=G0C2;
}

void MBTB_ContactSetDParam(unsigned int paramId,unsigned int contactId,unsigned int idShape,double v)
{
  assert(sNbOfContacts > contactId &&"MBTB_ContactLoadCADFile contactId out of range.");
  // unsigned int IdInCAD=sNbOfBodies+2*contactId;
  switch(paramId)
  {
  case 1:
    sContacts[contactId]->_Offset=v;
    break;
  case 2:
    sArtefactLength=v;
    break;
  case 3:
    sArtefactThreshold=v;
    break;
  case 4:
    sNominalForce=v;
    break;
  default:
    printf("Error: MBTB_ContactSetDParam paramId out of range \n");
  }
}
void MBTB_ContactSetIParam(unsigned int paramId,unsigned int contactId,unsigned int idShape,  int v)
{
  switch(paramId)
  {
  case 0:
    sContacts[contactId]->_OffsetP1=v;
    break;
  case 1:
    sContacts[contactId]->_normalFromFace1=v;
    break;
  case 2:
    sDrawMode=(unsigned int)v;
    break;
  default:
    printf("Error: MBTB_ContactSetIParam paramId out of range \n");
  }
}
void MBTB_BodySetDParam(unsigned int paramId,unsigned int bodyId,double v)
{
  printf("MBTB_BodySetDParam not yet implemented\n");
}

void MBTB_BodySetIParam(unsigned int paramId,unsigned int bodyId,int v)
{
  printf("MBTB_BodySetIParam not yet implemented\n");
}
void MBTB_BodySetVelocity(unsigned int numDS, SP::SiconosVector aVel)
{
  SP::SiconosVector v = sDS[numDS]->velocity();
  *v= *aVel;
  SP::SiconosVector v0 = sDS[numDS]->v0();
  *v0=*aVel;
}
void MBTB_SetDParam(unsigned int paramId,double v)
{
  sDParams[paramId]=v;
}
void MBTB_print_dist(unsigned int v)
{
  sPrintDist=v;
}

void MBTB_displayStep_bodies(unsigned int v)
{
  sDisplayStepBodies=v;
}
void MBTB_displayStep_joints(unsigned int v)
{
  sDisplayStepJoints=v;
}
void MBTB_displayStep_contacts(unsigned int v)
{
  sDisplayStepContacts=v;
}

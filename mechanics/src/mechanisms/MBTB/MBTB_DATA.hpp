/*! \addtogroup MBTB_DATA
   *  \brief This file contains the static memory of the MBTB module.
   *
   * The memory allocation is done using MBTB_MAX_BODIES_NUMBER, MBTB_MAX_JOINTS_NUMBER and MBTB_MAX_CONTACTS_NUMBER.
   *  @{
   */
#ifndef MBTB_DATA
#define MBTB_DATA
#include "MBTB_Body.hpp"
#include "MBTB_ContactRelation.hpp"
#include "MBTB_JointR.hpp"
//!Must be 1. It is the update frequency setting the transformation of the graphical object.
#define FREQ_UPDATE_GRAPHIC 1
//! The maximal number of bodies.
#define MBTB_MAX_BODIES_NUMBER 100
//! The maximal number of joints.
#define MBTB_MAX_JOINTS_NUMBER 100
//! The maximal number of contacts.
#define MBTB_MAX_CONTACTS_NUMBER 100
//!The dynamical bodies.
extern SP::MBTB_Body sDS[];
//!The joint relations.
extern MBTB_JointR * sJointRelations[];
//!The contacts.
extern MBTB_Contact * sContacts[];
//!The number of bodies.
extern unsigned int sNbOfBodies;
//!The number of joints.
extern unsigned int sNbOfJoints;
//!The number of contacts.
extern unsigned int sNbOfContacts;
//!The counter of step of simulation.
extern unsigned int sTimerCmp;
//!The graphical frequency.
extern unsigned int sFreqGraphic;
//!The output frequency.
extern unsigned int sFreqOutput;
//!The siconos joint interactions.
extern SP::Interaction sInterJoints[];
//!The siconos contact interactions.
extern SP::Interaction sInterContacts[];
//!siconos model.
extern SP::Model myModel;
//!siconos model t0.
extern double myt0;
//!siconos model Tf.
extern double myTf;
//!use the gravity vector
extern unsigned int sUseGravity;

//!for the graph building.
//! The dynamical systems involved in the joint 'numJ' have indices sJointIndexDS[2*numJ] and sJointIndexDS[2*numJ+1].
extern int sJointIndexDS[];
//!The type of joint see JOINTS_TYPE.
extern int sJointType[];
//!The siconos simulation.
extern SP::TimeStepping sSimu;
//!The draw mode of the artefacts (forces, normals). Used with bit to bit test with MBTB_CST.
extern unsigned int sDrawMode;
//!The verbose mode for print_dist
extern unsigned int sPrintDist;
//!The verbose mode for displayStep_bodies
extern unsigned int sDisplayStepBodies;
//!The verbose mode for displayStep_joints
extern unsigned int sDisplayStepJoints;
//!The verbose mode for displayStep_contacts
extern unsigned int sDisplayStepContacts;
//!The nominal length of an artefact.
extern double sArtefactLength;
//!The minimal length drawing.
extern double sArtefactThreshold;
//!The nominal forces.
extern double sNominalForce;
extern double sDParams[20];
#endif
/*! @} */

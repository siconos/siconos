#include "MBTB_DATA.hpp"

unsigned int sNbOfBodies=0;
unsigned int sNbOfJoints=0;
unsigned int sNbOfContacts=0;
unsigned int sTimerCmp=0;
unsigned int sFreqGraphic=100;
unsigned int sFreqOutput =100;
SP::MBTB_Body sDS[MBTB_MAX_BODIES_NUMBER];
MBTB_JointR* sJointRelations[MBTB_MAX_JOINTS_NUMBER];
MBTB_Contact* sContacts[MBTB_MAX_CONTACTS_NUMBER];
SP::Interaction sInterJoints[MBTB_MAX_JOINTS_NUMBER];
SP::Interaction sInterContacts[MBTB_MAX_CONTACTS_NUMBER];
SP::Model myModel;
double myt0;
double myTf;
int sJointIndexDS[2*MBTB_MAX_JOINTS_NUMBER];
int sJointType[MBTB_MAX_JOINTS_NUMBER];
SP::TimeStepping sSimu;
unsigned int sDrawMode=0;
unsigned int sPrintDist=0;
unsigned int sDisplayStepBodies=0;
unsigned int sDisplayStepJoints=0;
unsigned int sDisplayStepContacts=0;
double sArtefactLength=1.0;
double sArtefactThreshold=1e-7;
double sNominalForce=0;
double sDParams[20];

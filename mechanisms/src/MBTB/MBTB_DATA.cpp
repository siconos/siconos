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

unsigned int siconos::mechanisms::mbtb::data::sNbOfBodies = 0;
unsigned int siconos::mechanisms::mbtb::data::sNbOfJoints = 0;
unsigned int siconos::mechanisms::mbtb::data::sNbOfContacts = 0;
unsigned int siconos::mechanisms::mbtb::data::sTimerCmp = 0;
unsigned int siconos::mechanisms::mbtb::data::sFreqGraphic = 100;
unsigned int siconos::mechanisms::mbtb::data::sFreqOutput = 100;
std::shared_ptr<siconos::mechanisms::MBTB_Body>
    siconos::mechanisms::mbtb::data::sDS[MBTB_MAX_BODIES_NUMBER];
siconos::mechanisms::MBTB_JointR*
    siconos::mechanisms::mbtb::data::sJointRelations[MBTB_MAX_JOINTS_NUMBER];
// std::vector<siconos::mechanisms::MBTB_Contact>
//     siconos::mechanisms::mbtb::data::sContacts(MBTB_MAX_CONTACTS_NUMBER);
std::shared_ptr<siconos::modeling::Interaction>
    siconos::mechanisms::mbtb::data::sInterJoints[MBTB_MAX_JOINTS_NUMBER];
std::shared_ptr<siconos::modeling::Interaction>
    siconos::mechanisms::mbtb::data::sInterContacts[MBTB_MAX_CONTACTS_NUMBER];
std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem>
    siconos::mechanisms::mbtb::data::myNsds;
double siconos::mechanisms::mbtb::data::myt0;
double siconos::mechanisms::mbtb::data::myTf;
int siconos::mechanisms::mbtb::data::sJointIndexDS[2 * MBTB_MAX_JOINTS_NUMBER];
siconos::mechanisms::JointsType
    siconos::mechanisms::mbtb::data::sJointType[MBTB_MAX_JOINTS_NUMBER];
std::shared_ptr<siconos::simulation::TimeStepping>
    siconos::mechanisms::mbtb::data::sSimu;
bool siconos::mechanisms::mbtb::data::sDrawMode = false;
unsigned int siconos::mechanisms::mbtb::data::sPrintDist = 0;
unsigned int siconos::mechanisms::mbtb::data::sDisplayStepBodies = 0;
unsigned int siconos::mechanisms::mbtb::data::sDisplayStepJoints = 0;
unsigned int siconos::mechanisms::mbtb::data::sDisplayStepContacts = 0;
double siconos::mechanisms::mbtb::data::sArtefactLength = 1.0;
double siconos::mechanisms::mbtb::data::sArtefactThreshold = 1e-7;
double siconos::mechanisms::mbtb::data::sNominalForce = 0;
double siconos::mechanisms::mbtb::data::sDParams[20];

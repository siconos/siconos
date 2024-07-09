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

/*! \addtogroup MBTB_DATA
 *  \brief This file contains the static memory of the MBTB module.
 *
 * The memory allocation is done using MBTB_MAX_BODIES_NUMBER,
 * MBTB_MAX_JOINTS_NUMBER and MBTB_MAX_CONTACTS_NUMBER.
 *  @{
 */
#ifndef MBTB_DATA
#define MBTB_DATA

//! Must be 1. It is the update frequency setting the transformation of the
//! graphical object.
#define FREQ_UPDATE_GRAPHIC 1
//! The maximal number of bodies.
#define MBTB_MAX_BODIES_NUMBER 100
//! The maximal number of joints.
#define MBTB_MAX_JOINTS_NUMBER 100
//! The maximal number of contacts.
#define MBTB_MAX_CONTACTS_NUMBER 100

#include <limits>
#include <memory>
#include <vector>

namespace siconos::modeling {
class NonSmoothDynamicalSystem;
class Interaction;

}  // namespace siconos::modeling

namespace siconos::simulation {
class Simulation;
class TimeStepping;

}  // namespace siconos::simulation

namespace siconos::mechanisms {

class MBTB_Body;
class MBTB_JointR;
class MBTB_Contact;

/** The joint type.

PIVOT_0, involves one ds.
PIVOT_1, involves two ds.
PRISMATIC, not used, the connection to siconos is not done.
*/
enum class JointsType {
  Pivot0 = 0,
  Pivot1 = 1,
  Prismatic0 = 2,
  Prismatic1 = 3
};

// /** Artefact constants.

// This constant are used to display info at contact :
// <ul>
// <li> ArtefactP1P2 draws the detected contact points and a line that
// links them</li> <li> ArtefactReaction draws the reaction forces</li>
// <li> ArtefactNormal draws the unit normal vector at contact</li>
// <li> FaceNormal1 draws the unit normal of surface of contact</li>
// </ul>
// Use with bit to bit test.
// */
//   enum class MBTBConstant : int {
//   ArtefactP1P2 = 1,
//   ArtefactReaction = 2,
//   ArtefactNormal = 4,
//   FaceNormal1 = 5
// };

namespace mbtb::draw {

inline constexpr bool ArtefactP1P2{true};
inline constexpr bool ArtefactReaction{true};
inline constexpr bool ArtefactNormal{true};
inline constexpr bool FaceNormal1{true};

}  // namespace mbtb::draw

namespace mbtb::data {

//! The dynamical bodies.
inline std::shared_ptr<siconos::mechanisms::MBTB_Body>
    sDS[MBTB_MAX_BODIES_NUMBER];
//! The joint relations.
inline MBTB_JointR* sJointRelations[MBTB_MAX_JOINTS_NUMBER];
//! The contacts.
inline std::vector<std::shared_ptr<siconos::mechanisms::MBTB_Contact>>
    sContacts(MBTB_MAX_CONTACTS_NUMBER);
//! The number of bodies.
inline unsigned int sNbOfBodies{0};
//! The number of joints.
inline unsigned int sNbOfJoints{0};
//! The number of contacts.
inline unsigned int sNbOfContacts{0};
//! The counter of step of simulation.
inline unsigned int sTimerCmp{0};
//! The graphical frequency.
inline unsigned int sFreqGraphic{100};
//! The output frequency.
inline unsigned int sFreqOutput{100};
//! The siconos joint interactions.
inline std::shared_ptr<siconos::modeling::Interaction>
    sInterJoints[MBTB_MAX_JOINTS_NUMBER];
//! The siconos contact interactions.
inline std::shared_ptr<siconos::modeling::Interaction>
    sInterContacts[MBTB_MAX_CONTACTS_NUMBER];
//! siconos model.
inline std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> myNsds;
//! siconos model t0.
inline double myt0{0.};
//! siconos model Tf.
inline double myTf{std::numeric_limits<double>::max()};
//! use the gravity vector
inline unsigned int sUseGravity{0};

//! for the graph building.
//!  The dynamical systems involved in the joint 'numJ' have indices
//!  sJointIndexDS[2*numJ] and sJointIndexDS[2*numJ+1].
inline int sJointIndexDS[2 * MBTB_MAX_JOINTS_NUMBER];
//! The type of joint see JOINTS_TYPE.
inline siconos::mechanisms::JointsType sJointType[MBTB_MAX_JOINTS_NUMBER];
//! The siconos simulation.
inline std::shared_ptr<siconos::simulation::TimeStepping> sSimu;
//! The draw mode of the artefacts (forces, normals). Used with bit to bit test
//! with MBTB_CST.
inline bool sDrawMode{false};
//! The verbose mode for print_dist
inline unsigned int sPrintDist{0};
//! The verbose mode for displayStep_bodies
inline unsigned int sDisplayStepBodies{0};
//! The verbose mode for displayStep_joints
inline unsigned int sDisplayStepJoints{0};
//! The verbose mode for displayStep_contacts
inline unsigned int sDisplayStepContacts{0};
//! The nominal length of an artefact.
inline double sArtefactLength{1.0};
//! The minimal length drawing.
inline double sArtefactThreshold{1.e-7};
//! The nominal forces.
inline double sNominalForce{0};
inline double sDParams[20];
}  // namespace mbtb::data
}  // namespace siconos::mechanisms
#endif
/*! @} */

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

/*! \addtogroup MBTB_PYTHON_API
   *  \brief This file provides the python API of the MBTB module.
   *
   It provides :
   * - the functions to build the multi bobies system
   * - the functions to set internal parameters
   * - the functions to run the simulation.
   *
   *  @{
   */

#ifndef MBTB_PYTHON_API
#define MBTB_PYTHON_API

#include <memory>
#include <string>

#include "BoundaryCondition.hpp"
#include "MBTB_DATA.hpp"  // for JointsType enum
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

namespace siconos::modeling {
class NonSmoothDynamicalSystem;
}

namespace siconos::simulation {
class Simulation;
}

namespace siconos::mechanisms {
using FunctionT_V = std::function<void(double, Eigen::Ref<siconos::algebra::MapVectorType>)>;

/**  To initialize the MBTB library (no yet dynamical memory).

     \param [in] NumOfBodies : unsigned int NumOfBodies, the number of bodies.
     \param [in] NumOfJoints : unsigned int , the number of joints.
     \param [in] NumberOfContacts : unsigned int , the number of contacts.
*/
void MBTB_init(unsigned int NumOfBodies, unsigned int NumOfJoints,
               unsigned int NumberOfContacts);

/** To update the CADs model */
void MBTB_updateDSFromSiconos();

/** Loads the cad model of a body

    \param [in] numDS unsigned int , identifier of Body.
    \param [in] CADFile const std::string& , the cad file.
    \param [in] withGraphicModel if 0 the graphic model is not built
 */
void MBTB_BodyLoadCADFile(unsigned int numDS, const std::string& CADFile,
                          unsigned int withGraphicModel);

/** Build the MBTB_body and set to the initial postion

    This function build a mechanical body in the simulator

    \param [in] numDS an identifier of body.
    \param [in] BodyName the body name.
    \param [in] mass the mass
    \param [in] initPos  a R^7 vector representing the position (translation
    in R^3, vector in R^3, angle in R) that must be appyed after the load to
    get the initial position of the object

    \param [in] initCenterMass coordinate of the mass center in the  model
    \param [in] inertialMatrix matrix in R^{3,3}
    \param [in] pluginFextLib the path to the plugin library
    \param [in] pluginFextFct the name of the pluged fonction
    \param [in] pluginMextLib the path to the plugin library
    \param [in] pluginMextFct the name of the pluged fonction
    \param [in] pluginFintLib the path to the plugin library
    \param [in] pluginFintFct the name of the pluged fonction
    \param [in] pluginMintLib the path to the plugin library
    \param [in] pluginMintFct the name of the pluged fonction
    \param [in] pluginFintJacqLib the path to the plugin library
    \param [in] pluginFintJacqFct the name of the pluged fonction
    \param [in] pluginMintJacqLib path to the plugin library
    \param [in] pluginMintJacqFct name of the pluged fonction
    \param [in] pluginFintJacvLib path to the plugin library
    \param [in] pluginFintJacvFct name of the pluged fonction
    \param [in] pluginMintJacvLib path to the plugin library
    \param [in] pluginMintJacvFct name of the pluged fonction
    \param [in] pluginBoundaryConditionLib path to the plugin library
    \param [in] pluginBoundaryConditionFct name of the pluged fonction
    \param [in] boundaryConditionIndex the indices of the velocities
   prescribed by the boundary condition
 */
void MBTB_BodyBuild(
    unsigned int numDS, const std::string& BodyName, double mass,
    std::shared_ptr<siconos::algebra::SiconosVector> initPos,
    std::shared_ptr<siconos::algebra::SiconosVector> modelCenterMass,
    std::shared_ptr<siconos::algebra::SiconosMatrix> inertialMatrix,
    // const siconos::modeling::newton_euler::FunctionT_V& fext_func,
    // const siconos::modeling::newton_euler::FunctionT_V& mext_func,
    // const siconos::modeling::newton_euler::FunctionVQT_V& fint_func,
    // const siconos::modeling::newton_euler::FunctionVQT_V& mint_func,
    // const siconos::modeling::newton_euler::FunctionVQT_V& jacobianfint_qfunc,
    // const siconos::modeling::newton_euler::FunctionVQT_V& jacobianfint_twistfunc,
    // const siconos::modeling::newton_euler::FunctionVQT_V& jacobianmint_qfunc,
    // const siconos::modeling::newton_euler::FunctionVQT_V& jacobianmint_twistfunc,
    const std::string& pluginBoundaryConditionLib,
    const std::string& pluginBoundaryConditionFct,
    const siconos::modeling::BoundaryCondition::Indices& boundaryConditionIndex);

/** To build a joint.
 *
 *  \param [in] numJ  an identifier
 *  \param [in] JointName a string for the joint name
 *  \param [in] jointType see enum JOINTS_TYPE
 *  \param [in] indexDS1 index of the first body attached to the joint
 *  \param [in] indexDS2 index of the second body  attached
 *    to the joint(ignored for JOINTS_TYPE _0)
 *  \param [in] jointPosition a R^7 vector representing the position
 * (translation in R^3, vector in R^3, angle in R) with respect to the
 * indexDS1 body frame( So, it must be applied after the load of the first
 * body).
 *
 */
void MBTB_JointBuild(unsigned int numJ, const std::string& JointName, JointsType jointType,
                     unsigned int indexDS1, unsigned int indexDS2,
                     std::shared_ptr<siconos::algebra::SiconosVector> xsujointPosition);

/** To set the location where is computed the equivalente forces.

    It consists in defining the location of the points in view of compute the
    equivalente forces of the joint

    \param [in] numJ an identifier
    \param [in] G0C1 vector where C1 Contact points were the
    forces of joint must be computed
    \param [in] G0C2 vector , where C2 Contact
    points were the forces of joint must be computed.
 */
void MBTB_setJointPoints(unsigned int numJ,
                         std::shared_ptr<siconos::algebra::SiconosVector> G0C1,
                         std::shared_ptr<siconos::algebra::SiconosVector> G0C2);

/** Loads the cad model of a contact
 *
 *  It is the shape involved in the contact. Each shape is includes in a
 *   container. It contains either some faces or some edges
 *
 *  \param [in]  contactId identifier of Contact
 *  \param [in]  CADFile1 cad file 1, it must contain one or two faces
 *  \param [in]  CADFile2 cad file 2, it must contains one or two faces,
 *    or one or two edges
 *  \param [in]  withGraphicModel1 1 to draw the corresponding object else 0
 *  \param [in]  withGraphicModel2 1 to draw the corresponding object else 0
 */
void MBTB_ContactLoadCADFile(unsigned int contactId, const std::string& CADFile1,
                             const std::string& CADFile2, unsigned int withGraphicModel1,
                             unsigned int withGraphicModel2);

/** To set a double parameter.(extendable, without modifying the API)

    This type of function has been chosen to easely set any parameters without
    modify the module API

    \param  [in] paramId : identifier of the param.<br>
    1 for offset.<br>
    2 for artefact lenghth.<br>
    3 for artefact thershold.<br>
    4 for the nominal force.<br>
    \param  [in] contactId : identifier of the contact.
    \param  [in] idShape : identifier of the shape of the contact (0 or 1).
    \param  [in] v : value.
 */
void MBTB_ContactSetDParam(unsigned int paramId, unsigned int contactId, unsigned int idShape,
                           double v);

/** To set a integer parameter.(extendable, without modifying the API)

    This type of function has been chosen to easely set any parameters without
    modify the module API

    \param paramId : identifier of the param.<br>
    0 for translate offset P1 parameters.<br>
    1 normal from face 1.<br>
    2 Artefact lenght.<br>
    \param contactId : identifier of the contact.
    \param idShape : identifier of the shape of the contact (0 or 1).
    \param v : value.
 */
void MBTB_ContactSetIParam(unsigned int paramId, unsigned int contactId, unsigned int idShape,
                           bool v);

/** To build a contact.

    It builds a relation of the convenient type doing the connection between
   the CAD model and the simulator

    \param [in] numContact the id of the contact
    \param [in] ContactName const name of the contact
    \param [in] indexBody1 id of the body carrying the first contact shape
    \param [in] indexBody2 id of the body carrying the second contact shape
    \param [in] withFriction true to set friction
    \param [in] mu
    \param [in] en
    \param [in] et (not used).
 */
void MBTB_ContactBuild(unsigned int numContact, const std::string& ContactName,
                       unsigned int indexBody1, int indexBody2, unsigned int withFriction,
                       double mu, double en, double et);

/** initializes the simulation.

    \param [in] hTS time step size.
    \param [in] withProj if 0 the projection in done.
*/
void MBTB_initSimu(double hTS, int withProj);

/** \return the NSDS */
std::shared_ptr<siconos::modeling::NonSmoothDynamicalSystem> MBTB_nsds();

/** \retrun the simulation */
std::shared_ptr<siconos::simulation::Simulation> MBTB_simulation();

/** Runs the simulation.

  \param [in] nbSteps int , the number of steps to be run
*/
void MBTB_run(int nbSteps);

/** Warm start.

  It sets the siconos state.

  \param [in] numDS the id of the ds.
  \param [in] aPos the target position.
  \param [in] aVel the target velocity.
 */
void MBTB_moveBodyToPosWithSpeed(unsigned int numDS,
                                 std::shared_ptr<siconos::algebra::SiconosVector> aPos,
                                 std::shared_ptr<siconos::algebra::SiconosVector> aVel);

// /** Set the velocity.

//     \param [in] numDS the id of the ds.
//     \param [in] aVel the targeted velocity.
//  */
// void MBTB_BodySetVelocity(unsigned int numDS,
//                           std::shared_ptr<siconos::algebra::SiconosVector> aVel);

/** Defines the graphic frequency.

    \param [in] freq the frequency.
*/
void MBTB_setGraphicFreq(unsigned int freq);

/** Defines the graphic frequency.

    \param [in] freq the frequency.
*/
void MBTB_setOutputFreq(unsigned int freq);

/** Sets an integer value of the solver's parameters.

  \param [in] i index of the parameter.
  \param [in] value  the new value.
*/
void MBTB_setSolverIOption(int i, int value);

/** Sets a double value of the solver's parameters.

    \param [in]  i index of the parameter.
    \param [in]  value the new value.
*/
void MBTB_setSolverDOption(int i, double value);

/** To enable/disable the projection algorithm.

    This function is usefull only when proj has been activated in MBTB_init.

    \param [in] v if 0 then the projection is disabled.

*/
void MBTB_doProj(unsigned int v);

/** To perform only the projection algorithm.

  This function is usefull only when proj has been activated in MBTB_init.

  \param [in] v if 0 then only the projection algorithm is
  done, the mechanical equations are not simulated.

 */
void MBTB_doOnlyProj(unsigned int v);

/** To set the max iteration of the projection algorithm.

  This function is usefull only when proj has been activated in MBTB_init.

  \param [in]  v the max number of iteration

 */
void MBTB_projectionMaxIteration(unsigned int v);

/** To set the tolerance on constraints (joints) of the projection

    This function is usefull only when proj has been activated in MBTB_init.

    \param [in]  tol tolerance

 */
void MBTB_constraintTol(double tol);

/** To set the tolerance on unilateral constraints (contact) of the projection
   algorithm.

    This function is usefull only when proj has been activated in MBTB_init.

    \param [in]  tol v tolerance

 */
void MBTB_constraintTolUnilateral(double tol);

/** Verbose bodies mode

    \param [in] v if 0 no verbose

 */
void MBTB_displayStep_bodies(unsigned int v);

/** Verbose joints mode

    \param [in] v if 0 no verbose

 */
void MBTB_displayStep_joints(unsigned int v);

/** Verbose contacts mode

    \param [in] v if 0 no verbose

*/
void MBTB_displayStep_contacts(unsigned int v);

/** Verbose distance and contact mode

    \param [in] v if 0 no verbose

 */
void MBTB_print_dist(unsigned int v);

/** MBTB_BodySetDParam not yet used
 * \param paramId
 * \param bodyId
 * \param v
 */
void MBTB_BodySetDParam(unsigned int paramId, unsigned int bodyId, double v);

/** Must be call before MBTB_InitSimu
    It sets the mbtb::data::sDParams value:
    0 : MBTB_TimeStepping::_deactivateYPosThreshold;
    1 : MBTB_TimeStepping::_deactivateYVelThreshold;
    2 : MBTB_TimeStepping::_activateYPosThreshold;
    3 : MBTB_TimeStepping::_activateYVelThreshold;

    \param paramId
    \param v
 */
void MBTB_SetDParam(unsigned int paramId, double v);

/** MBTB_BodySetDParam not yet used
 *
 * \param paramId
 * \param bodyId
 * \param v
 */
void MBTB_BodySetIParam(unsigned int paramId, unsigned int bodyId, int v);

}  // namespace siconos::mechanisms
#endif
/*! @} */

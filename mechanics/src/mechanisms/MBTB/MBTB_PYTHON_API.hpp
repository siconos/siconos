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
#include "MBTB_Body.hpp"
//! To initialize the MBTB library. It calls the CADMBTB_init.
/*!
  To initialize the MBTB library (no yet dynamical memory).
  \param [in] NumOfBodies : unsigned int NumOfBodies, the number of bodies.
  \param [in] NumOfJoints : unsigned int , the number of joints.
  \param [in] NumberOfContacts : unsigned int , the number of contacts.
 */
void MBTB_init(unsigned int NumOfBodies, unsigned int NumOfJoints, unsigned int NumberOfContacts);

//! updating the CADs model
/*! MBTB_updateDSFromSiconos
  It consists in updating the CADs model form the siconos quaternion.
 */

void MBTB_updateDSFromSiconos();

//! Load the cad model of a body
/*!
 Load the cad model of body numDS
  \param [in] numDS unsigned int , identifier of Body.
  \param [in] CADFile const std::string& , the cad file.
  \param [in] withGraphicModel unsigned int , iff 0 the graphic model is not built
 */
void MBTB_BodyLoadCADFile(unsigned int numDS,const std::string& CADFile,unsigned int withGraphicModel);
//! Build the MBTB_body and set to the initial postion.
/*!
 This function build a mechanical body in the simulator.
  \param [in] numDS unsigned int , an identifier of body.
  \param [in] BodyName const std::string& , a string for the body name.
  \param [in] mass double , the mass.
  \param [in] initPos SP::SiconosVector , a R^7 vector representing the position (translation in R^3, vector in R^3, angle in R) that must be appyed after the load to get the initial position of the object.
  \param [in] initCenterMass SP::SiconosVector , coordinate of the mass center in the just loaded model
  \param [in] inertialMatrix SP::SimpleMatrix , matrix in R^{3,3}
  \param [in] pluginFextLib const std::string& , the path to the plugin library.
  \param [in] pluginFextFct const std::string& , the name of the pluged fonction.
  \param [in] pluginMextLib const std::string& , the path to the plugin library.
  \param [in] pluginMextFct const std::string& , the name of the pluged fonction.
  \param [in] pluginFintLib const std::string& , the path to the plugin library.
  \param [in] pluginFintFct const std::string& , the name of the pluged fonction.
  \param [in] pluginMintLib const std::string& , the path to the plugin library.
  \param [in] pluginMintFct const std::string& , the name of the pluged fonction.
  \param [in] pluginFintJacqLib const std::string& , the path to the plugin library.
  \param [in] pluginFintJacqFct const std::string& , the name of the pluged fonction.
  \param [in] pluginMintJacqLib const std::string& , the path to the plugin library.
  \param [in] pluginMintJacqFct const std::string& , the name of the pluged fonction.
  \param [in] pluginFintJacvLib const std::string& , the path to the plugin library.
  \param [in] pluginFintJacvFct const std::string& , the name of the pluged fonction.
  \param [in] pluginMintJacvLib const std::string& , the path to the plugin library.
  \param [in] pluginMintJacvFct const std::string& , the name of the pluged fonction.
  \param [in] pluginBoundaryConditionLib const std::string& , the path to the plugin library.
  \param [in] pluginBoundaryConditionFct const std::string& , the name of the pluged fonction.
  \param [in] boundaryConditionIndex SP::IndexInt, the indices of the velocities prescribed by the boundary condition
 */
void MBTB_BodyBuild(unsigned int numDS, const std::string& BodyName,double mass,
                    SP::SiconosVector initPos, SP::SiconosVector initCenterMass,
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
                    SP::IndexInt boundaryConditionIndex);
//! To build a joint.
/*!
 * It builds the joint in the simulator.
 * \param [in] numJ unsigned int , an identifier .
 * \param [in] JointName const std::string& , a string for the joint name.
 * \param [in] jointType unsigned int , see enum JOINTS_TYPE.
 * \param [in] indexDS1 unsigned int , index of the first body attached to the joint.
 * \param [in] indexDS2 unsigned int , index of the second body attached to the joint(ignored for JOINTS_TYPE _0).
 * \param [in] jointPosition SP::SiconosVector  a R^7 vector representing the position (translation in R^3, vector in R^3, angle in R) with respect to the indexDS1 body frame( So, it must be applied after the load of the first body).
 *
 */
void MBTB_JointBuild(unsigned int numJ, const std::string& JointName,unsigned int jointType,
                     unsigned int indexDS1, unsigned int indexDS2, SP::SiconosVector jointPosition);
//!  To set the location where is computed the equivalente forces.
/*!
  It consists in defining the location of the points in view of compute the equivalente forces of the joint.
  \param [in] numJ unsigned int , an identifier .
  \param [in] G0C1 vector , where C1 Contact points where the forces of joint must be computed.
  \param [in] G0C2 vector , where C2 Contact points where the forces of joint must be computed.
 */
void MBTB_setJointPoints(unsigned int numJ, SP::SiconosVector G0C1,SP::SiconosVector G0C2);

//!MBTB_ContactLoadCADFile load the cad model of a contact
/*!
 *It is the shape involved in the contact. Each is include in a container. It contains either some faces or some edges.
 * \param [in]  contactId unsigned int , identifier of Contact.
 * \param [in]  CADFile1 const std::string& , the cad file 1, it must contain one or two faces
 * \param [in]  CADFile2 const std::string& , the cad file 2, it must contains one or two faces, or one or two edges.
 * \param [in]  withGraphicModel1  unsigned int : 1 to draw the corresponding object else 0;
 * \param [in]  withGraphicModel2 unsigned int : 1 to draw the corresponding object else 0;
 */
void MBTB_ContactLoadCADFile(unsigned int contactId,const std::string& CADFile1,const std::string& CADFile2,unsigned int withGraphicModel1,unsigned int withGraphicModel2);

//!To set a double parameter.(extendable, without modifying the API)
/**
 This type of function has been chosen to easely set any parameters without modify the module API.
 *\param  [in] paramId : identifier of the param.<br>
 1 for offset.<br>
 2 for artefact lenghth.<br>
 3 for artefact thershold.<br>
 4 for the nominal force.<br>
 *\param  [in] contactId : identifier of the contact.
 *\param  [in] idShape : identifier of the shape of the contact (0 or 1).
 *\param  [in] v : value.
 */
void MBTB_ContactSetDParam(unsigned int paramId,unsigned int contactId,unsigned int idShape,double v);

//!To set a integer parameter.(extendable, without modifying the API)
/*!
   This type of function has been chosen to easely set any parameters without modify the module API.
 *\param paramId : identifier of the param.<br>
 0 for translate offset P1 parameters.<br>
 1 normal from face 1.<br>
 2 Artefact lenght.<br>
 *\param contactId : identifier of the contact.
 *\param idShape : identifier of the shape of the contact (0 or 1).
 *\param v : value.
 */
void MBTB_ContactSetIParam(unsigned int paramId,unsigned int contactId,unsigned int idShape, int v);

//! To build a contact.
/*!
 *It builds a relation of the convenient type doing the connection between the CAD model and the simulator.
 \param [in] numContact unsigned int , the id os the contact.
 \param [in] ContactName const std::string& , the name of the contact.
 \param [in] indexBody1 unsigned int , the id of the body carrying the first contact shape.
 \param [in] indexBody2 unsigned int , the id of the body carrying the second contact shape.
 \param [in] withFriction unsigned int , 0 or 1.
 \param [in] mu double.
 \param [in] en double.
 \param [in] et double  (not used).
 */
void MBTB_ContactBuild(unsigned int numContact, const std::string& ContactName, unsigned int indexBody1, int indexBody2, unsigned int withFriction, double mu, double en, double et);
//! It initializes the simulation.
/*!
 It consists in building the siconos simulation and al.
 \param [in] hTS double , time step size.
 \param [in] withProj int , iff 0 the projection in done.
*/
void  MBTB_initSimu(double hTS, int withProj);

//! Get Siconos model.
/*!
  The model may be used outside MBTB in Siconos Front-End.
  \return  SP::Model
*/
SP::NonSmoothDynamicalSystem MBTB_nsds();

//! It runs the simulation.
/*!
  It consists in running nbSteps simulation steps.
  \param [in] nbSteps int , the number of run step.
 */
void MBTB_run(int nbSteps);

//! It is a warm start.
/*!
  It sets the siconos state.
  \param [in] numDS unsigned int , the id of the ds.
  \param [in] aPos SP::SiconosVector , the target position.
  \param [in] aVel SP::SiconosVector , the target velocity.
 */
void MBTB_moveBodyToPosWithSpeed(unsigned int numDS, SP::SiconosVector aPos, SP::SiconosVector aVel);

//! Set the velocity.
/*!
  It sets the siconos state.
  \param [in] numDS unsigned int , the id of the ds.
  \param [in] aVel SP::SiconosVector , the target velocity.
 */
void MBTB_BodySetVelocity(unsigned int numDS, SP::SiconosVector aVel);

//!It defines the graphic frequency.
/*!
  \param [in] freq unsigned int, the frequency.
 */
void MBTB_setGraphicFreq(unsigned int freq);

//!It defines the graphic frequency.
/*!
  \param [in] freq unsigned int, the frequency.
 */
void MBTB_setOutputFreq(unsigned int freq);

//!It sets an integer value of the solver's parameters.
/*!
  \param [in] i index of the parameter.
  \param [in] value  the new value.
 */
void MBTB_setSolverIOption(int i,int value);

//!It sets a double value of the solver's parameters.
/*!
  \param [in]  i index of the parameter.
  \param [in]  value the new value.
 */
void MBTB_setSolverDOption(int i,double value);

//! It allows to enable/disable the projection algorithm.
/*!
  This function is usefull only when proj has been activated in MBTB_init.
  \param [in] v unsigned int  iff 0 then the projection is disabled.

 */
void MBTB_doProj(unsigned int v);

//! It allows to perform only the projection algorithm.
/*!
  This function is usefull only when proj has been activated in MBTB_init.
  \param [in] v unsigned int , iff 0 then only the projection algorithm is done, the mechanical equations are not simulated.

 */
void MBTB_doOnlyProj(unsigned int v);


//! It allows to set the max iteration of the projection algorithm.
/*!
  This function is usefull only when proj has been activated in MBTB_init.
  \param [in]  v unsigned int  the max number of iteration

 */
void MBTB_projectionMaxIteration(unsigned int v);

//! It allows to set the tolerance on constraints (joints) of the projection algorithm.
/*!
  This function is usefull only when proj has been activated in MBTB_init.
  \param [in]  tol double  tolerance

 */
void MBTB_constraintTol(double tol);

//! It allows to set the tolerance on unilateral constraints (contact) of the projection algorithm.
/*!
  This function is usefull only when proj has been activated in MBTB_init.
  \param [in]  tol double v tolerance

 */
void MBTB_constraintTolUnilateral(double tol);








//! It allows to verbose  bodies information
/*!
  \param [in] v unsigned int , if 0 no verbose

 */
void MBTB_displayStep_bodies(unsigned int v);

//! It allows to verbose  joints information
/*!
  \param [in] v unsigned int, if 0 no verbose

 */
void MBTB_displayStep_joints(unsigned int v);

//! It allows to verbose  contacts information
/*!
  \param [in] v unsigned int , if 0 no verbose

 */

void MBTB_displayStep_contacts(unsigned int v);

//! It allows to verbose distance and contact information
/*!
  \param [in] v unsigned int , if 0 no verbose

 */
void MBTB_print_dist(unsigned int v);



/** MBTB_BodySetDParam not yet used
 * \param paramId
 * \param bodyId
 * \param v
 */
void MBTB_BodySetDParam(unsigned int paramId,unsigned int bodyId ,double v);

//! MBTB_SetDParam
/*!
  It must be call before MBTB_InitSimu
  It sets the sDParams value:
  0 : MBTB_TimeStepping::_deactivateYPosThreshold;
  1 : MBTB_TimeStepping::_deactivateYVelThreshold;
  2 : MBTB_TimeStepping::_activateYPosThreshold;
  3 : MBTB_TimeStepping::_activateYVelThreshold;
  \param paramId
  \param v
 */
void MBTB_SetDParam(unsigned int paramId,double v);


/** MBTB_BodySetDParam not yet used
 *
 * \param paramId
 * \param bodyId
 * \param v
 */
void MBTB_BodySetIParam(unsigned int paramId, unsigned int bodyId, int v);

//! The joint type.
/*!
  PIVOT_0, involves one ds.
  PIVOT_1, involves two ds.
  PRISMATIC, not used, the connection to siconos is not done.
 */
enum JOINTS_TYPE
{
  PIVOT_0=0,
  PIVOT_1=1,
  PRISMATIC_0=2,
  PRISMATIC_1=3
};
//!Artefact constants.
/*!
  This constant are used to display info at contact :
  <ul>
  <li> MBTB_ARTEFACT_P1P2 draws the detected contact points and a line that links them</li>
  <li> MBTB_ARTEFACT_Reaction draws the reaction forces</li>
  <li> MBTB_ARTEFACT_NORMAL draws the unit normal vector at contact</li>
  <li> MBTB_FACE_NORMAL draws the unit normal of surface of contact</li>
  </ul>
  Use with bit to bit test.
 */
enum MBTB_CST
{
  MBTB_ARTEFACT_P1P2 = 1,
  MBTB_ARTEFACT_REACTION = 2,
  MBTB_ARTEFACT_NORMAL=4,
  MBTB_FACE_NORMAL1=5
};
#endif
/*! @} */

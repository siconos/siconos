//$Id: DynamicalSystem.h,v 1.61 2005/03/14 16:05:26 jbarbier Exp $
#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include <string>
#include <vector>

#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "SimpleVector.h"
#include "CompositeVector.h"
#include "SiconosMemory.h"
#include "SiconosSharedLibrary.h"
#include "NonSmoothDynamicalSystem.h"
#include "DSInputOutput.h"
#include "BoundaryCondition.h"
#include "DSXML.h"
#include "SiconosConst.h"
#include "RuntimeException.h"
#include "check.h"

//#include "XMLTagsName.h"

using namespace std;

class NonSmoothDynamicalSystem;
class BoundaryCondition;
class DSInputOutput;

class DSXML;

/** \class DynamicalSystem
 *  \brief  Super class of the dynamical systems
 *  \author V. ACARY,  JB CHARLETY
 *  \version 1.0
 *  \date (Creation) April 29, 2004
 *
 *  $Date: 2005/03/14 16:05:26 $
 *  $Revision: 1.61 $
 *  $Author: jbarbier $
 *  $Source: /CVS/Siconos/SICONOS/src/modelformalisation/DynamicalSystem.h,v $
 *
 * The class DynamicalSystem allows to define and compute a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * \dot x = f(x,t)+r,
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *
 *  The function \f$ f : R^{n} \times R  \mapsto  R^{n}   \f$ defines the VectorField.
 *
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial condition is given by
 *  * \f[
 *  x(t_0)=x_0
 * \f]
 * To define an boundary Value Problem, the pointer on  a BoundaryCondition must be set.
 *
 * One word on the bilateral constraint
 *
 *
 * Particular cases such as linear system (LinearSystemDS) or
 * Lagrangian Non Linear System (LagrangianNLDS)  are specialization of this class.
 *
 * \todo Add a pointer to an object Constraint .
 */
class DynamicalSystem
{
public:

  DynamicalSystem();

  /** \fn DynamicalSystem( DSXML* )
   *  \brief constructor with XML object of the DynamicalSystem
   *  \param DSXML* : the XML object corresponding to the DynamicalSystem
   */
  DynamicalSystem(DSXML*);

  virtual ~DynamicalSystem();

  /*getter et setter*/

  /** \fn NonSmoothDynamicalSystem* getNSDS(void);
   *  \brief allows to get the NonSmoothDynamicalSystem containing the DynamicalSystem
   *  \return NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem containing the DynamicalSystem
   */
  inline NonSmoothDynamicalSystem* getNSDS(void) const
  {
    return this->nsds;
  }

  /** \fn int getNumber(void);
   *  \brief allows to get the number of the DynamicalSystem
   *  \return the value of number
   */
  inline int getNumber(void) const
  {
    return this->number;
  }

  /** \fn string getId(void)
   *  \brief allows to get the id of the DynamicalSystem
   *  \return the value of ths id
   */
  inline string getId(void) const
  {
    return this->id;
  }

  /** \fn int getN(void);
   *  \brief allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline int getN(void) const
  {
    return this->n;
  }


  /** \fn SiconosVector* getX0(void)
   *  \brief allows to get x0, the state at the initial time t0 of the DynamicalSystem
   *  \return pointer on x0
   */
  inline SiconosVector* getX0(void) const
  {
    return this->x0;
  }

  /** \fn SiconosVector* getX(void)
   *  \brief allows to get the state vector x of the dynamical system
   *  \return pointer on vector x
   */
  inline SiconosVector* getX(void) const
  {
    return this->x;
  }

  /** \fn SiconosMemory* getXMemories(void)
   *  \brief allows to get all the values of the state vector x stored in memory
   *  \return the memory object which stores previous values of x
   */
  inline SiconosMemory* getXMemories(void)
  {
    return &this->xMemory;
  }

  /** \fn SiconosVector& getXDotPtr(void)
   *  \brief allow to get the memory adress of vector xDot
   *  \exception to be defined
   *  \return SiconosVector* xDot
   *  \warning such type of method dont have to be included in the final API
   */
  inline SiconosVector* getXDotPtr(void)
  {
    return &(this->xDot);
  }

  /** \fn SimpleVector getXDot(void)
   *  \brief allows to get the vector xDot
   *  \return a SimpleVector
   */
  inline SimpleVector getXDot(void) const
  {
    return this->xDot;
  }

  /** \fn SiconosMemory* getXDotMemories(void)
   *  \brief allows to get all the values of old xDot vectors
   *  \return the memory object which stores previous values of xDot
   */
  inline SiconosMemory* getXDotMemories(void)
  {
    return &this->xDotMemory;
  }


  /** \fn SiconosVector* getXFree(void)
   *  \brief allows to get the vector xFree
   *  \return pointer on SiconosVector xFree
   */
  inline SiconosVector* getXFree(void)
  {
    return this->xFree;
  }


  /** \fn SiconosVector* getRPtr(void)
   *  \brief allow to get the memory adress of vector r
   *  \exception to be defined
   *  \return SiconosVector* r
   *  \warning such type of method dont have to be included in the final API
   */
  inline SiconosVector* getRPtr(void)
  {
    return &(this->r);
  }


  /** \fn SimpleVector getR(void)
   *  \brief allow to get the vector r
   *  \return SimpleVector r
   */
  inline /*SiconosVector*/SimpleVector getR(void) const
  {
    return this->r;
  }

  /** \fn SiconosMemory* getRMemories(void)
   *  \brief allows to get all the values of old r vectors
   *  \return the memory object which stores previous values of r
   */
  inline SiconosMemory* getRMemories(void)
  {
    return &this->rMemory;
  }


  /** \fn SiconosMatrix getJacobianX(void)
   *  \brief allows to get the SiconosMatrix jacobianXMat
   *  \return the SiconosMatrix jacobianXMat
   */
  inline SiconosMatrix getJacobianX(void) const
  {
    return this->jacobianX;
  }

  /** \fn int getStepsInMemory(void)
   *  \brief allows to get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline int getStepsInMemory(void) const
  {
    return this->stepsInMemory;
  }

  /** \fn BoundaryCondition PeriodicBCXML(void)
   *  \brief allows to get the BoundaryCondition
   *  \return a pointer on the BoundaryCondition object
   */
  inline BoundaryCondition* getBoundaryCondition(void)
  {
    return this->BC;
  }



  /** \fn void setNSDS(NonSmoothDynamicalSystem*);
   *  \brief allows to set the NonSmoothDynamicalSystem containing the DynamicalSystem
   *  \param NonSmoothDynamicalSystem* : the NonSmoothDynamicalSystem containing the DynamicalSystem
   */
  inline void setNSDS(NonSmoothDynamicalSystem*nsds)
  {
    this->nsds = nsds;
  }

  /** \fn void setNumber(int)
   *  \brief allows to set the value of number
   *  \param int number : an integer to set the value of number
   */
  inline void setNumber(int number)
  {
    this->number = number;
  }

  /** \fn void setId(string)
   *  \brief allows to set the value of id
   *  \param string id : a string to set the value of id
   */
  inline void setId(string id)
  {
    this->id = id;
  }

  /** \fn void setN(int)
   *  \brief allows to set the value of n
   *  \param int n : an integer to set the value of n
   */
  inline void setN(int n)
  {
    this->n = n;
  }

  /** \fn void setX0(SiconosVector*)
   *  \brief allows to set the value of x0
   *  \param SiconosVector* x0 : a pointer on vector to set the value of x0
   */
  inline void setX0(SiconosVector *x0)
  {
    this->x0 = x0;
  }

  /** \fn void setX(SiconosVector*)
   *  \brief allows to set the value of x
   *  \param SiconosVector x : a vector to set the value of x
   */
  inline void setX(SiconosVector *x)
  {
    this->x = x;
  }

  /** \fn void setXMemories(SiconosMemory &)
   *  \brief allows to set the value of xMemory
   *  \param SiconosMemory & xMem : a memory object.
   */
  inline void setXMemories(SiconosMemory &xMem)
  {
    this->xMemory = xMem;
  }

  /** \fn void setXDot(SimpleVector &)
   *  \brief allows to set the value of xDot
   *  \param SimpleVector & xDot : a SimpleVector to set the value of xDot
   */
  inline void setXDot(SimpleVector &xDot)
  {
    this->xDot = xDot;
  }

  /** \fn void setXDotMemories(SiconosMemory &)
   *  \brief allows to set the value of xDotMemory
   *  \param SiconosMemory &xDotMem : a memory object.
   */
  inline void setXDotMemories(SiconosMemory &xDotMem)
  {
    this->xDotMemory = xDotMem;
  }

  /** \fn void setXFree(SiconosVector*)
   *  \brief allows to set the value of xFree
   *  \param SiconosVector *xFree : a pointer on SiconosVector to set the value of xFree
   */
  inline void setXFree(SiconosVector *xFree)
  {
    this->xFree = xFree;
  }

  /** \fn void setR(SimpleVector)
   *  \brief allows to set the value of r
   *  \param SiconosVector &r : a SimpleVector to set the value of r
   */
  inline void setR(SimpleVector &r)
  {
    this->r = r;
  }

  /** \fn void setRMemories(SiconosMemory &)
   *  \brief allows to set the value of rMemory
   *  \param SiconosMemory &rMem : a memory object
   */
  inline void setRMemories(SiconosMemory &rMem)
  {
    this->rMemory = rMem;
  }

  /** \fn void setJacobianX(SiconosMatrix)
   *  \brief allows to set the SiconosMatrix jacobianXMat
   *  \param SiconosMatrix jacobian : the matrix to set jacobianX
   */
  inline void setJacobianX(SiconosMatrix jacobian)
  {
    this->jacobianX = jacobian;
  }

  /** \fn void setStepsInMemory(int)
   *  \brief allows to set the value of stepsInMemory
   *  \param int steps : the value to set stepsInMemory
   */
  inline void setStepsInMemory(int steps)
  {
    this->stepsInMemory = steps;
  }

  /** \fn void setBoundaryCondition(BoundaryCondition*)
   *  \brief allows to set the BoundaryCondition
   *  \param BoundaryCondition *bc : the BoundaryCondition to set BC
   */
  inline void setBoundaryCondition(BoundaryCondition *bc)
  {
    this->BC = bc;
  }

  /** \fn inline DSXML* getDynamicalSystemXML()
   *  \brief allows to get the object DSXML of the DynamicalSystem
   *  \return a pointer on the DSXML of the DynamicalSystem
   */
  inline DSXML* getDynamicalSystemXML() const
  {
    return this->dsxml;
  }

  /** \fn inline void getDynamicalSystemXML(DSXML *dsxml)
   *  \brief allows to set the DSXML of the DynamicalSystem
   *  \param DSXML* dsxml : the address of theDSXML to set
   */
  inline void setDynamicalSystemXML(DSXML *dsxml)
  {
    this->dsxml = dsxml;
  }

  /** \fn void setVectorFieldFunction(string, string)
   *  \brief allow to set a specified function to compute vector field
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setVectorFieldFunction(string pluginPath, string functionName);

  /** \fn void setComputeJacobianXFunction(string, string)
   *  \brief allow to set a specified function to compute jacobianX
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianXFunction(string pluginPath, string functionName);


  /** \fn inline string getType()
   *  \brief allows to get the type of a DynamicalSystem
   *  \return string : the type of the DynamicalSystem
   */
  inline string getType() const
  {
    return this->DSType;
  }

  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the DynamicalSystem
   *  \return the vector of DSInputOutput
   */
  vector<DSInputOutput*> getDSInputOutputs(void);

  /** \fn DSInputOutput* getDSInputOutput(int)
   *  \brief allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
   *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the DynamicalSystem
   *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
   */
  DSInputOutput* getDSInputOutput(int);

  /** \fn void setDSInputOutputs(vector<DSInputOutput*>)
   *  \brief allows to set all the DSInputOutputs of the DynamicalSystem
   *  \param vector<DSInputOutput*> : the vector to set
   */
  void setDSInputOutputs(vector<DSInputOutput*>);

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the DynamicalSystem
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  void addDSInputOutput(DSInputOutput*);

  ////////////////////////////////////

  /** \fn void initMemory(int steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  virtual void initMemory(int steps) ;

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, xDot and r in the stored previous values
   *  xMemory, xDotMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  virtual void swapInMemory(void);

  /** \fn void vectorField (double time)
   * \brief Default function for computing the vector field \f$ f: (x,t) \in R^{n} \times R  \mapsto  R^{n}\f$
   * \param double time : the time for the computation
   *  \exception RuntimeException
   */
  virtual void vectorField(double time);

  /** \fn static void computeJacobianX (double time)
   *  \brief Default function for computing the gradient of the vector field with the respect
   *  to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : the time for the computation
   *  \exception RuntimeException
   */
  virtual void computeJacobianX(double time);

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** \fn void display()
   *  \brief print the data of the dynamical system on the standard output
   */
  void display() const;

  /** \fn void createDynamicalSystem(DSXML * nsdsXML, int number, int n,
      SiconosVector* x0, string vectorFieldPlugin, NonSmoothDynamicalSystem * nsds, BoundaryCondition* bc)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML* : the XML object for this DynamicalSystem
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial state of this DynamicalSystem
   *  \param string : the plugin name for vectorField of this DynamicalSystem
   *  \exception RuntimeException
   */
  void createDynamicalSystem(DSXML * dsXML, int number = -1, int n = -1,
                             SiconosVector* x0 = NULL, string vectorFieldPlugin = "BasicPlugin:vectorField");
  //, NonSmoothDynamicalSystem * nsds = NULL, BoundaryCondition* bc=NULL);

  /** \fn BoundaryCondition* createPeriodicBC()
   *  \brief create the Periodic Boundary Condition of this DynamicalSystem
   */
  BoundaryCondition* createPeriodicBC();

  /** \fn BoundaryCondition* createLinearBC( SiconosMatrix*, SiconosMatrix*, SiconosMatrix* )
   *  \brief create the Linear Boundary Condition of this DynamicalSystem
   *  \param SiconosVector* : the omega vector of this boundary condition
   *  \param SiconosMatrix* : the omega0 matrix of this boundary condition
   *  \param SiconosMatrix* : the omegaT matrix of this boundary condition
   */
  BoundaryCondition* createLinearBC(SiconosVector* omega = NULL,
                                    SiconosMatrix* omega0 = NULL, SiconosMatrix* omegaT = NULL);

  /** \fn BoundaryCondition* createNLinearBC()
   *  \brief create the NLinear Boundary Condition of this DynamicalSystem
   */
  BoundaryCondition* createNLinearBC();

  /** \var typedef void (*vfPtr) (int* sizeOfX, double* time, double* xPtr, double* xdotPtr);
   *  \brief signature of plugin function computing the vectorfield
   *  \param int* sizeOfX : the size of the vector X
   *  \param double* time : the time for the computation
   *  \param double* xPtr : the pointer to the first element of the vector X
   *  \param double* jacobianXPtr : the pointer to the first element of the matrix jacobianX (in-out parameter)
   */
  typedef void (*vfPtr)(int* sizeOfX, double* time, double* xPtr, double* xdotPtr);

  /** \fn vfPtr getVectorFieldPtr()
   *  \brief return the function adress of the plugin computing vectorfield
   */
  vfPtr getVectorFieldPtr()
  {
    return *vectorFieldPtr;
  }


protected:
  /** \fn void fillDSWithDSXML()
   *  \brief uses the DSXML of the DynamicalSystem to fill the field of this DynamicalSystem
   *  \exception RuntimeException
   */
  virtual void fillDSWithDSXML();

  /** \fn void linkDSXML()
   *  \brief makes the links between the BoundaryConditionXML of the DSXML of the DynamicalSystem and the BoundaryCondition
   */
  virtual void linkDSXML();

  /** the type of the DS : LagrangianNLDS, LagrangianTIDS, LinearSystemDS */
  string DSType;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** this number defines in a single way the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc.)*/
  string id;

  /** the dimension of the system (i.e. size of the state vector x, or the vector r, ...)*/
  int n;

  /** x0 the value of the state at the initial time t0 */
  SiconosVector *x0;

  /** the state vector of the system, \f$  x \in R^{n}\f$ */
  SiconosVector *x;

  /** the  previous state vectors stored in memory*/
  SiconosMemory xMemory; //vector< SiconosVector* > xMemory;

  /** the time derivative of the state x (the velocity) */
  /*SiconosVector*/
  SimpleVector xDot;

  /** the  previous vector of the time derivative of the state x  */
  SiconosMemory xDotMemory; //vector< SiconosVector* > xDotMemory;

  /** the  free state vector (the state vector for r=0) */
  SiconosVector *xFree;

  /** the  input vector due to the non-smooth law \f$  r \in R^{n}\f$ (multiplier, force, ...)*/
  /*SiconosVector*/
  SimpleVector r;

  /**  the previous vectors, r, stored in memory */
  SiconosMemory rMemory; //vector< SiconosVector* > rMemory;

  /** number of previous states stored in memory */
  int stepsInMemory;

  /* contains the name of the plugin used to compute the vectorField */
  string vectorFieldFunctionName;
  /* contains the name of the plugin used to compute the JacobianX */
  string computeJacobianXFunctionName;

  /** Gradient of the vectorfield \f$ f(x,t) \f$ with the respect to \f$ x\f$*/
  SiconosMatrix jacobianX;

  /** the boundary conditions defined if the DynamicalSystem has boundary conditions */
  BoundaryCondition *BC;

  /** the XML object linked to the DynamicalSystem to read xML data */
  DSXML *dsxml;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** adress of the plugin function computing vectorfield */
  vfPtr vectorFieldPtr;


  /** contains a link to the DSInputOutput of the DynamicalSystems */
  vector<DSInputOutput*> dsioVector;

  /** \fn void (*computeJacobianXPtr) (int* sizeOfX, double* time, double* xPtr, double* jacobianXPtr)
   *  \brief  Pointer on function for computing the gradient of the vector field with the respect to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param int* sizeOfX : the size of the vector X
   *  \param double* time : the time for the computation
   *  \param double* xPtr : the pointer to the first element of the vector X
   *  \param double* jacobianXPtr : the pointer to the first element of the matrix jacobianX (in-out parameter)
   */
  void (*computeJacobianXPtr)(int* sizeOfX, double* time, double* xPtr, double* jacobianXPtr);

  /** \fn void init()
   *  \brief initialise value of a Dynamical System
   */
  virtual void init();


private :

};

#endif // DYNAMICALSYSTEM_H


//$Log: DynamicalSystem.h,v $
//Revision 1.61  2005/03/14 16:05:26  jbarbier
//- manual creation of DSInputOutput saving OK
//
//- in progress for EqualityConstraint
//
//Revision 1.60  2005/03/08 12:41:36  jbarbier
//- constant variables files modified :
//Some constants added in SiconosConst
//
//all global tag of the modeling tools are in XMLTagsName, other tags are specific to an XML class
//
//Revision 1.59  2005/03/07 13:17:19  jbarbier
//- new test : Ball2D, with a ball moving in a 2D system
//
//- another constant variables moved/refactored in XMLTagsName
//- making uniform the name of the constant variables
//
//Revision 1.58  2005/02/11 17:35:55  charlety
//
//_ little "inspection of code"
//_ basic getters and setters passed inline
//_ getters functions passed const
//
//Revision 1.57  2005/02/10 10:35:18  jbarbier
//- new file regrouping all the const values of the model, modelingTools and numericalStrategy
//
//- new function in the LagrangianLinearR to get the H matrix corresponding to one of the 2 dynamical systems linked to the relation
//
//- new atribute of the OneStepNSProblem. A visibility table of the Interaction.
//
//Revision 1.56  2005/01/25 14:51:46  jbarbier
//- attributes id, type and XML object added to EqualityConstraint
//
//Revision 1.55  2005/01/25 13:56:03  jbarbier
//- link DynamicalSystem-DSInputOutput, NonSmoothDynamicalSystem-EqualityConstraint, EquaityConstraint-DSInputOutput and Relation-DSInputOutput available
//
//Revision 1.54  2005/01/20 14:44:48  jbarbier
//- NSDS class renamed NonSmoothDynamicalSystem
//
//- code reduce, some comments remove
//
//Revision 1.53  2004/09/30 08:35:01  jbarbier
//- fonction of the formalisation : fill..With...XML and link... are now
//"protected" and no more "public"
//
//Revision 1.52  2004/09/22 11:16:28  charlety
//
//_ revision of Doxygen comments in modelformalisation
//
//Revision 1.51  2004/09/16 11:35:24  jbarbier
//- save of the TimeDiscretisation in a XML file in manual creation of the
//platform which was forgotten is now available.
//
//- the save of the platform's data can be done when the platform is created with
//an XML input file and completed with dynmical systems, interactions, one-step
//non smooth problem and one-step integrator.
//
//Revision 1.50  2004/09/10 11:26:07  charlety
//
//_ Integration of the new version of the SiconosVector in the platform. the class simpleVector is used mostly to replace old SiconosVector. When a vector can be composite or simple, like the state of a dynamical system, a pointer on SiconosVector is used, and the vector is initialized simple or composite when the system is initialized.
//
//_ All the tests which worked with the previous version of the vector are OK with the new version.
//
//_ Example SICONOS and bouncingBall are OK
//
//_ some comments have still to be adapted to NewSiconosVector .
//
//_ This version of NewSiconosVector could be called 0.9. some details have to be fixed, it will be done before the end of September.
//
//Revision 1.49  2004/09/10 08:04:46  jbarbier
//- XML save available for BoundaryCondition and Interaction
//
//Revision 1.48  2004/09/03 14:41:41  jbarbier
//- new functions to create the boundary condition of the dynamical systems
//- new functions to add an interaction to a NonSmoothDynamicalSystem
//- new functions to create the relation and the non-smooth law of an interaciton
//
//Revision 1.47  2004/08/20 15:26:44  jbarbier
//- creation of a Model and save in the XML is ok
//- creation of a NonSmoothDynamicalSystem and save in the XML is ok
//- creation of a NonLinearSystemDS and save in the XML is OK
//
//Revision 1.46  2004/08/17 15:12:37  jbarbier
//- methods createDynamicalSystem, createBoundaryCondition, createInteraction,
//createRelation and createNSLaw completed with the required attributes
//
//Revision 1.45  2004/08/13 11:26:58  jbarbier
//- function createNSDS complete
//
//- function createDynamicalSystem and createLinearSystemDS complete
//
//- function  createLagrangianNLDS in progress
//
//Revision 1.44  2004/08/10 13:04:15  charlety
//
//_ try to pass a pointer of a C function of th eplugin to a numerical routine
//
//Revision 1.43  2004/08/05 12:44:40  jbarbier
//- loading XML file with no OneStepNSProblem succesfull
//
//- NonLinearSystemDS is now available
//
//Revision 1.42  2004/07/23 14:39:25  jbarbier
//- createModel, createNSDS, createDynamicalSystem, createBoundaryCondition OK
//it's the first step, it do the same thing that before, but the process is
//unified and it must simply add the traitment for the creation of the nodes in
//the DOM tree
//
//Revision 1.41  2004/07/09 14:18:55  jbarbier
//-management of the optional attributes of the DynamicalSystem
//-new node t for current time in the time of the NonSmoothDynamicalSystem
//
//Revision 1.40  2004/07/09 11:14:53  charlety
//
//_ Added a constructor by copy and an operator = in class SiconosMemory
//_ the getters on memory in DynamicalSystems return now some pointers
//
//Revision 1.39  2004/06/30 13:08:21  jbarbier
//DoxygenSiconos.cfg moved into config/
//ICDLLSharedLibrary renamed SiconosSharedLibrary
//
//Revision 1.38  2004/06/29 08:49:57  acary
//Ajout des commentaires Doxygen et des tages CVS
//

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

//using namespace std;

class NonSmoothDynamicalSystem;
class BoundaryCondition;
class DSInputOutput;

class DSXML;

/** \class DynamicalSystem
 *  \brief  Super class of the dynamical systems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) April 29, 2004
 *
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
 * Lagrangian Non Linear System (LagrangianDS)  are specialization of this class.
 *
 * \todo Add a pointer to an object Constraint .
 */
class DynamicalSystem
{
public:

  // --- Constructors ---

  /** \fn DynamicalSystem(DSXML * nsdsXML)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML* : the XML object for this DynamicalSystem
   *  \exception RuntimeException
   */

  DynamicalSystem(DSXML * dsXML);

  /** \fn DynamicalSystem(DSXML * nsdsXML, int number, int n,
      SiconosVector* x0, string vectorFieldPlugin, NonSmoothDynamicalSystem * nsds, BoundaryCondition* bc)
      *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
      *  \param int : the number for this DynamicalSystem
      *  \param int : the dimension of this DynamicalSystem
      *  \param SiconosVector* : the initial state of this DynamicalSystem
      *  \param string : the plugin name for vectorField of this DynamicalSystem
      *  \exception RuntimeException
      */

  DynamicalSystem(int number, int n,
                  SiconosVector* x0, string vectorFieldPlugin = "BasicPlugin:vectorField");

  // --- Destructor ---

  virtual ~DynamicalSystem();

  // ---  Getters and setters ---

  /** \fn NonSmoothDynamicalSystem* getNSDSPtr(void) const;
   *  \brief get the NonSmoothDynamicalSystem containing this DynamicalSystem
   *  \return NonSmoothDynamicalSystem*
   */
  inline NonSmoothDynamicalSystem* getNSDSPtr(void) const
  {
    return this->nsds;
  }

  /** \fn void setNSDSPtr(NonSmoothDynamicalSystem*);
   *  \brief set the NonSmoothDynamicalSystem containing the DynamicalSystem
   *  \param NonSmoothDynamicalSystem*
   */
  inline void setNSDSPtr(NonSmoothDynamicalSystem *newNsds)
  {
    nsds = newNsds;
  }

  /** \fn const int getNumber(void) const;
   *  \brief allows to get the number of the DynamicalSystem
   *  \return the value of number
   */
  inline const int getNumber(void) const
  {
    return this->number;
  }

  /** \fn void setNumber(const int&)
   *  \brief allows to set the value of number
   *  \param an integer to set the value of number
   */
  inline void setNumber(const int& newNumber)
  {
    this->number = newNumber;
  }

  /** \fn const string getId(void) const
   *  \brief allows to get the id of the DynamicalSystem
   *  \return the value of ths id
   */
  inline const string getId(void) const
  {
    return this->id;
  }

  /** \fn void setId(const string&)
   *  \brief allows to set the value of id
   *  \param a string to set the value of id
   */
  inline void setId(const string& newId)
  {
    this->id = newId;
  }

  /** \fn const int getN(void) const;
   *  \brief allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline const int getN(void) const
  {
    return this->n;
  }

  /** \fn void setN(const int&)
   *  \brief allows to set the value of n
   *  \param an integer to set the value of n
   */
  inline void setN(const int& newN)
  {
    this->n = newN;
  }

  // --- X0 ---

  /** \fn  const SimpleVector getX0(void) const
   *  \brief get the value of x0, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getX0() const
  {
    return *(this->x0);
  }

  /** \fn SiconosVector* getX0Ptr(void) const
   *  \brief get x0, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getX0Ptr(void) const
  {
    return this->x0;
  }

  /** \fn void setX0(const SiconosVector& newValue)
   *  \brief set the value of x0 to newValue
   *  \param SiconosVector newValue
   */
  inline void setX0(const SiconosVector& newValue)
  {
    *(this->x0) = newValue;
  }

  /** \fn void setX0Ptr(SiconosVector* newPtr)
   *  \brief set x0 to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  inline void setX0Ptr(SiconosVector* newPtr)
  {
    delete x0;
    x0 = 0;
    this->x0 = newPtr;
  }

  // --- X ---

  /** \fn const SimpleVector getX(void) const
   *  \brief get the value of x, the state of the DynamicalSystem
   *  \return SimpleVector
    * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getX(void) const
  {
    return *(this->x);
  }

  /** \fn SiconosVector* getXPtr(void) const
   *  \brief get x, the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXPtr(void) const
  {
    return this->x;
  }

  /** \fn void setX (const SiconosVector& newValue)
   *  \brief set the value of x to newValue
   *  \param SiconosVector newValue
   */
  inline void setX(const SiconosVector& newValue)
  {
    *(this->x) = newValue;
  }

  /** \fn void setXPtr(SiconosVector* newPtr)
   *  \brief set x to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  inline void setXPtr(SiconosVector *newPtr)
  {
    delete x ;
    x = 0 ;
    this->x = newPtr;
  }

  // ---  XDot ---

  /** \fn  const SimpleVector getXDot(void) const
   *  \brief get the value of xDot derivative of the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getXDot(void) const
  {
    return *(this->xDot);
  }

  /** \fn SiconosVector* getXDotPtr(void) const
   *  \brief get xDot, the derivative of the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXDotPtr(void) const
  {
    return this->xDot;
  }

  /** \fn void setXDot (const SiconosVector& newValue)
   *  \brief set the value of xDot to newValue
   *  \param SiconosVector newValue
   */
  inline void setXDot(const SiconosVector& newValue)
  {
    *(this->xDot) = newValue;
  }

  /** \fn void setXDotPtr(SiconosVector* newPtr)
   *  \brief set xDot to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  inline void setXDotPtr(SiconosVector *newPtr)
  {
    delete xDot;
    xDot = 0;
    this->xDot = newPtr;
  }

  // --- XFree ---

  /** \fn  const SimpleVector getXFree(void) const
   *  \brief get the value of xFree
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getXFree(void) const
  {
    return *(this->xFree);
  }

  /** \fn SiconosVector* getXFreePtr(void) const
   *  \brief get xFree
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXFreePtr(void) const
  {
    return this->xFree;
  }

  /** \fn void setXFree (const SiconosVector& newValue)
   *  \brief set the value of xFree to newValue
   *  \param SiconosVector newValue
   */
  inline void setXFree(const SiconosVector& newValue)
  {
    *(this->xFree) = newValue;
  }

  /** \fn void setXFreePtr(SiconosVector* newPtr)
   *  \brief set xFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  inline void setXFreePtr(SiconosVector *newPtr)
  {
    delete xFree;
    xFree = 0;
    this->xFree = newPtr;
  }

  // --- R ---

  /** \fn  const SimpleVector getR(void) const
   *  \brief get the value of r
   *  \return SimpleVector
   */
  inline const SimpleVector getR(void) const
  {
    return *(this->r);
  }

  /** \fn SimpleVector* getRPtr(void) const
   *  \brief get r
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getRPtr(void) const
  {
    return this->r;
  }

  /** \fn void setR (const SimpleVector& newValue)
   *  \brief set the value of r to newValue
   *  \param SimpleVector newValue
   */
  inline void setR(const SimpleVector& newValue)
  {
    *(this->r) = newValue;
  }

  /** \fn void setRPtr(SimpleVector* newPtr)
   *  \brief set R to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setRPtr(SimpleVector *newPtr)
  {
    delete r;
    r = 0;
    this->r = newPtr;
  }

  // --- Memory ---

  /** \fn const int getStepsInMemory(void) const
   *  \brief allows to get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline const int getStepsInMemory(void) const
  {
    return this->stepsInMemory;
  }

  /** \fn void setStepsInMemory(const int&)
   *  \brief allows to set the value of stepsInMemory
   *  \param int steps : the value to set stepsInMemory
   */
  inline void setStepsInMemory(const int& steps)
  {
    this->stepsInMemory = steps;
  }

  // X memory

  /** \fn  const SiconosMemory getXMemory(void) const
   *  \brief get the value of xMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getXMemory(void) const
  {
    return this->xMemory;
  }

  /** \fn SiconosMemory getXMemoryPtr(void) const
   *  \brief allows to get all the values of the state vector x stored in memory
   *  \return the memory object which stores previous values of x
   */
  inline SiconosMemory* getXMemoryPtr(void)
  {
    return &(this->xMemory);
  }

  /** \fn void setXMemory(const SiconosMemory &)
   *  \brief set the value of xMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setXMemory(const SiconosMemory& xMem)
  {
    this->xMemory = xMem;
  }

  // xDot memory

  /** \fn  const SiconosMemory getXDotMemory(void) const
   *  \brief get the value of xDotMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getXDotMemory(void) const
  {
    return this->xDotMemory;
  }

  /** \fn SiconosMemory* getXDotMemoryPtr(void)
   *  \brief get all the values of old xDot vectors
   *  \return a pointer on the memory object which contains previous values of xDot
   */
  inline SiconosMemory* getXDotMemoryPtr(void)
  {
    return &(this->xDotMemory);
  }

  /** \fn void setXDotMemory(const SiconosMemory &)
   *  \brief set the value of xDotMemory
   *  \param SiconosMemory &xDotMem : a memory object.
   */
  inline void setXDotMemory(const SiconosMemory& xDotMem)
  {
    this->xDotMemory = xDotMem;
  }

  // r memory
  /** \fn  const SiconosMemory getRMemory(void) const
   *  \brief get the value of RMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getRMemory(void) const
  {
    return this->rMemory;
  }

  /** \fn SiconosMemory* getRMemoryPtr(void)
   *  \brief get all the values of old r vectors
   *  \return the memory object which contains previous values of r
   */
  inline SiconosMemory* getRMemoryPtr(void)
  {
    return &(this->rMemory);
  }

  /** \fn void setRMemory(const SiconosMemory &)
   *  \brief set the value of rMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setRMemory(const SiconosMemory& rMem)
  {
    this->rMemory = rMem;
  }

  // --- JacobianX ---

  /** \fn  const SiconosMatrix getJacobianX(void) const
   *  \brief get the value of JacobianX
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianX(void) const
  {
    return *(this->jacobianX);
  }

  /** \fn SiconosMatrix* getJacobianXPtr(void) const
   *  \brief get JacobianX
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianXPtr(void) const
  {
    return this->jacobianX;
  }

  /** \fn void setJacobianX (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianX to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setJacobianX(const SiconosMatrix& newValue)
  {
    *(this->jacobianX) = newValue;
  }

  /** \fn void setJacobianXPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianX to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setJacobianXPtr(SiconosMatrix *newPtr)
  {
    delete jacobianX;
    jacobianX = 0;
    jacobianX = newPtr;
  }

  // --- Boundary Conditions ---
  //\todo : getter and setter to be reviewed when implement BC correctly
  /** \fn  const BoundaryCondition getBoundaryCondition(void) const
   *  \brief get the value of BoundaryCondition
   *  \return an object BoundaryCondition
   */
  //  inline BoundaryCondition getBoundaryCondition(void) { return &(this->BC); }

  /** \fn BoundaryCondition getBoundaryConditionPtr(void) const
   *  \brief get the BoundaryCondition
   *  \return a pointer on the BoundaryCondition object
   */
  inline BoundaryCondition* getBoundaryConditionPtr(void) const
  {
    return this->BC;
  }

  /** \fn void setBoundaryCondition(const BoundaryCondition&)
   *  \brief set the Boundary Conditions
   *  \param ref on an object BoundaryCondition
   */
  //inline void setBoundaryCondition(const BoundaryCondition& newBC) {*(this->BC) = newBC; }

  /** \fn void setBoundaryConditionPtr(BoundaryCondition*)
   *  \brief set the BoundaryCondition pointer
   *  \param BoundaryCondition *bc : the BoundaryCondition to set BC
   */
  inline void setBoundaryConditionPtr(BoundaryCondition *newBC)
  {
    delete BC;
    BC = 0;
    BC = newBC;
  }

  // --- dsxml ---
  /** \fn inline const DSXML* getDynamicalSystemXMLPtr() const
   *  \brief get the object DSXML of the DynamicalSystem
   *  \return a pointer on the DSXML of the DynamicalSystem
   */
  inline const DSXML* getDynamicalSystemXMLPtr() const
  {
    return this->dsxml;
  }

  /** \fn inline void setDynamicalSystemXMLPtr(DSXML *dsxml)
   *  \brief set the DSXML of the DynamicalSystem
   *  \param DSXML* dsxml : the address of theDSXML to set
   */
  inline void setDynamicalSystemXMLPtr(DSXML *newDsxml)
  {
    dsxml = newDsxml;
  }


  // --- Vector field ---
  /** \fn void setVectorFieldFunction(const string&, const string&)
   *  \brief allow to set a specified function to compute vector field
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setVectorFieldFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeJacobianXFunction(const string&, const string&)
   *  \brief allow to set a specified function to compute jacobianX
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianXFunction(const string& pluginPath, const string& functionName);

  // --- type of DS ---
  /** \fn inline string getType()
   *  \brief allows to get the type of a DynamicalSystem
   *  \return string : the type of the DynamicalSystem
   */
  inline const string getType() const
  {
    return this->DSType;
  }

  // --- DS input-output ---
  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the DynamicalSystem
   *  \return the vector of DSInputOutput
   */
  inline vector<DSInputOutput*> getDSInputOutputs(void)
  {
    return dsioVector;
  }

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
  inline void setDSInputOutputs(vector<DSInputOutput*> newDsioVect)
  {
    this->dsioVector = newDsioVect;
  }

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the DynamicalSystem
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  void addDSInputOutput(DSInputOutput* dsio)
  {
    this->dsioVector.push_back(dsio);
  }

  // --- ---

  /** \fn void initMemory(const int& steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  virtual void initMemory(const int& steps) ;

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, xDot and r in the stored previous values
   *  xMemory, xDotMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  virtual void swapInMemory(void);

  /** \fn void vectorField (const double& time)
   * \brief Default function for computing the vector field \f$ f: (x,t) \in R^{n} \times R  \mapsto  R^{n}\f$
   * \param double time : the time for the computation
   *  \exception RuntimeException
   */
  virtual void vectorField(const double& time);

  /** \fn static void computeJacobianX (const double& time)
   *  \brief Default function for computing the gradient of the vector field with the respect
   *  to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : the time for the computation
   *  \exception RuntimeException
   */
  virtual void computeJacobianX(const double& time);

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** \fn void saveDSDataToXML()
   *  \brief copy the data common to each system in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSDataToXML();

  /** \fn void saveBCToXML()
   *  \brief copy the Boundary Conditions data in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveBCToXML();

  /** \fn void saveDSIOToXML()
   *  \brief copy the DS Input-Output data in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSIOToXML();

  /** \fn void display()
   *  \brief print the data of the dynamical system on the standard output
   */
  void display() const;

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
  BoundaryCondition* createLinearBC(SiconosVector* omega = 0,
                                    SiconosMatrix* omega0 = 0, SiconosMatrix* omegaT = 0);

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

  /** \fn double dsConvergenceIndicator()
   *  \brief Default function for computing an indicator of convergence
   *   \brief return a double
   */
  virtual double dsConvergenceIndicator() const ;

protected:

  // Default constructor
  DynamicalSystem();

  /** \fn void fillBoundaryConditionsFromXml()
   *  \brief uses the DSXML of the DynamicalSystem to fill BoundaryCondition fields
   *  \exception RuntimeException
   */
  virtual void fillBoundaryConditionsFromXml();

  /** \fn void fillDsioFromXml()
   *  \brief uses the DSXML of the DynamicalSystem to fill DSIO vector
   *  \exception RuntimeException
   */
  virtual void fillDsioFromXml();

  /** Dynamical System type: General Dynamical System (NLDS) LagrangianDS (LNLDS),
      LagrangianLinearTIDS (LTIDS), LinearSystemDS (LDS)*/
  string DSType;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** this number defines in a single way the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc.)*/
  string id;

  /** the dimension of the system (i.e. size of the state vector x, or the vector r, ...)*/
  int n;

  /** initial state of the system */
  SiconosVector *x0;

  /** state of the system, \f$  x \in R^{n}\f$ */
  SiconosVector *x;

  /** the  previous state vectors stored in memory*/
  SiconosMemory xMemory;

  /** the time derivative of the state x (the velocity) */
  SiconosVector *xDot;

  /** the  previous xDot vectors */
  SiconosMemory xDotMemory;

  /** the  free state vector (state vector for r=0) */
  SiconosVector *xFree;

  /** the  input vector due to the non-smooth law \f$  r \in R^{n}\f$ (multiplier, force, ...)*/
  SimpleVector *r;

  /**  the previous r vectors */
  SiconosMemory rMemory;

  /** number of previous states stored in memory */
  int stepsInMemory;

  /* the name of the plugin used to compute the vectorField */
  string vectorFieldFunctionName;

  /* the name of the plugin used to compute the JacobianX */
  string computeJacobianXFunctionName;

  /** Gradient of the vectorfield \f$ f(x,t) \f$ with respect to \f$ x\f$*/
  SiconosMatrix *jacobianX;

  /** boundary conditions defined if the DynamicalSystem has some */
  BoundaryCondition *BC;

  /** the XML object linked to the DynamicalSystem  */
  DSXML *dsxml;

  /** class for plugin managing (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** adress of the plugin function that computes vectorfield */
  vfPtr vectorFieldPtr;

  /** vector of the DS Inputs-Outputs of the Dynamical System */
  vector<DSInputOutput*> dsioVector;

  /** \fn void (*computeJacobianXPtr) (int* sizeOfX, double* time, double* xPtr, double* jacobianXPtr)
   *  \brief  Pointer on function to compute the gradient of the vector field with the respect to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param int* sizeOfX : size of vector X
   *  \param double* time : time for computation
   *  \param double* xPtr : pointer to the first element of X
   *  \param double* jacobianXPtr : pointer to the first element of jacobianX matrix (in-out parameter)
   */
  void (*computeJacobianXPtr)(int* sizeOfX, double* time, double* xPtr, double* jacobianXPtr);

private :
};

#endif // DYNAMICALSYSTEM_H



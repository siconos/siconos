#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosConst.h"
#include "RuntimeException.h"
#include "check.h"

#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "SiconosMemory.h"
#include "SiconosSharedLibrary.h"

#include "NonSmoothDynamicalSystem.h"
#include "DSInputOutput.h"
#include "BoundaryCondition.h"
#include "DSXML.h"

#include <string>
#include <vector>

class NonSmoothDynamicalSystem;
class BoundaryCondition;
class DSInputOutput;
class DSXML;
class SiconosVector;
class SiconosMatrix;
class SiconosMemory;
class SiconosSharedLibrary;

/** \class DynamicalSystem
 *  \brief  General first order non linear dynamical systems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) April 29, 2004
 *
 *
 * This class defines and computes a generic n-dimensional
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
 * and the initial conditions are given by
 *  * \f[
 *  x(t_0)=x_0
 * \f]
 * To define an boundary Value Problem, the pointer on  a BoundaryCondition must be set.
 *
 * \todo One word on the bilateral constraint
 *
 *
 * Particular cases such as linear system (LinearDS) or
 * Lagrangian Non Linear System (LagrangianDS)  are specialization of this class.
 *
 * \todo Add a pointer to an object Constraint .
 */
class DynamicalSystem
{
public:

  // --- Constructors ---

  /** \fn DynamicalSystem(DSXML * nsdsXML, NonSmoothDynamicalSystem* =NULL)
   *  \brief xml constructor
   *  \param DSXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  DynamicalSystem(DSXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** \fn DynamicalSystem(DSXML * nsdsXML, const unsigned int& number, const unsigned int& n,
      const SiconosVector& x0, const string& vectorFieldPlugin)
      *  \brief constructor from a set of data
      *  \param int : reference number for this DynamicalSystem
      *  \param int : dimension of this DynamicalSystem
      *  \param SiconosVector : initial state of this DynamicalSystem
      *  \param string : plugin name for vectorField of this DynamicalSystem
      *  \exception RuntimeException
      */
  DynamicalSystem(const int&, const unsigned int&,
                  const SiconosVector&, const std::string& = "BasicPlugin:vectorField");

  /** \fn DynamicalSystem(const DynamicalSystem &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  DynamicalSystem(const DynamicalSystem &);

  // --- Destructor ---
  virtual ~DynamicalSystem();

  // ---  Getters and setters ---

  /** \fn NonSmoothDynamicalSystem* getNSDSPtr(void) const;
   *  \brief get the NonSmoothDynamicalSystem containing this DynamicalSystem
   *  \return NonSmoothDynamicalSystem*
   */
  inline NonSmoothDynamicalSystem* getNSDSPtr() const
  {
    return nsds;
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
  inline const int getNumber() const
  {
    return number;
  }

  /** \fn void setNumber(const int&)
   *  \brief allows to set the value of number
   *  \param an integer to set the value of number
   */
  inline void setNumber(const int& newNumber)
  {
    number = newNumber;
  }

  /** \fn const string getId(void) const
   *  \brief allows to get the id of the DynamicalSystem
   *  \return the value of ths id
   */
  inline const std::string  getId() const
  {
    return id;
  }

  /** \fn void setId(const string&)
   *  \brief allows to set the value of id
   *  \param a string to set the value of id
   */
  inline void setId(const std::string & newId)
  {
    id = newId;
  }

  /** \fn const int getN(void) const;
   *  \brief allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline const unsigned int getN(void) const
  {
    return n;
  }

  /** \fn void setN(const int&)
   *  \brief allows to set the value of n
   *  \param an integer to set the value of n
   */
  inline void setN(const int& newN)
  {
    n = newN;
  }

  // --- X0 ---

  /** \fn  const SimpleVector getX0(void) const
   *  \brief get the value of x0, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getX0() const
  {
    return *x0;
  }

  /** \fn SiconosVector* getX0Ptr(void) const
   *  \brief get x0, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getX0Ptr() const
  {
    return x0;
  }

  /** \fn void setX0(const SiconosVector& newValue)
   *  \brief set the value of x0 to newValue
   *  \param SiconosVector newValue
   */
  inline void setX0(const SiconosVector& newValue)
  {
    *x0 = newValue;
  }

  /** \fn void setX0Ptr(SiconosVector* newPtr)
   *  \brief set x0 to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setX0Ptr(SiconosVector*);

  // --- X ---

  /** \fn const SimpleVector getX(void) const
   *  \brief get the value of x, the state of the DynamicalSystem
   *  \return SimpleVector
    * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getX() const
  {
    return *x;
  }

  /** \fn SiconosVector* getXPtr(void) const
   *  \brief get x, the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXPtr() const
  {
    return x;
  }

  /** \fn void setX (const SiconosVector& newValue)
   *  \brief set the value of x to newValue
   *  \param SiconosVector newValue
   */
  inline void setX(const SiconosVector& newValue)
  {
    *x = newValue;
  }

  /** \fn void setXPtr(SiconosVector* newPtr)
   *  \brief set x to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXPtr(SiconosVector *);

  // X memory

  /** \fn  const SiconosMemory getXMemory(void) const
   *  \brief get the value of xMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getXMemory() const
  {
    return *xMemory;
  }

  /** \fn SiconosMemory getXMemoryPtr(void) const
   *  \brief get all the values of the state vector x stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getXMemoryPtr() const
  {
    return xMemory;
  }

  /** \fn void setXMemory(const SiconosMemory &)
   *  \brief set the value of xMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setXMemory(const SiconosMemory& newValue)
  {
    *xMemory = newValue;
  }

  /** \fn void setXMemory(SiconosMemory * newPtr)
   *  \brief set xMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setXMemoryPtr(SiconosMemory *);

  // ---  XDot ---

  /** \fn  const SimpleVector getXDot(void) const
   *  \brief get the value of xDot derivative of the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getXDot() const
  {
    return *xDot;
  }

  /** \fn SiconosVector* getXDotPtr(void) const
   *  \brief get xDot, the derivative of the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXDotPtr() const
  {
    return xDot;
  }

  /** \fn void setXDot (const SiconosVector& newValue)
   *  \brief set the value of xDot to newValue
   *  \param SiconosVector newValue
   */
  inline void setXDot(const SiconosVector& newValue)
  {
    *xDot = newValue;
  }

  /** \fn void setXDotPtr(SiconosVector* newPtr)
   *  \brief set xDot to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXDotPtr(SiconosVector *);

  // XDot memory

  /** \fn  const SiconosMemory getXDotMemory(void) const
   *  \brief get the value of xDotMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getXDotMemory() const
  {
    return *xDotMemory;
  }

  /** \fn SiconosMemory getXDotMemoryPtr(void) const
   *  \brief get all the values of the state vector xDot stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getXDotMemoryPtr() const
  {
    return xDotMemory;
  }

  /** \fn void setXDotMemory(const SiconosMemory &)
   *  \brief set the value of xDotMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setXDotMemory(const SiconosMemory& newValue)
  {
    *xDotMemory = newValue;
  }

  /** \fn void setXDotMemory(SiconosMemory * newPtr)
   *  \brief set xDotMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setXDotMemoryPtr(SiconosMemory *);

  // --- XFree ---

  /** \fn  const SimpleVector getXFree(void) const
   *  \brief get the value of xFree
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getXFree() const
  {
    return *xFree;
  }

  /** \fn SiconosVector* getXFreePtr(void) const
   *  \brief get xFree
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXFreePtr() const
  {
    return xFree;
  }

  /** \fn void setXFree (const SiconosVector& newValue)
   *  \brief set the value of xFree to newValue
   *  \param SiconosVector newValue
   */
  inline void setXFree(const SiconosVector& newValue)
  {
    *xFree = newValue;
  }

  /** \fn void setXFreePtr(SiconosVector* newPtr)
   *  \brief set xFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXFreePtr(SiconosVector *);

  /** \fn  const SimpleVector* getUPtr(void) const
   *  \brief get the value of u -> since u is not yet implemented in DS, this is only an interface for LinearDS
   *  \return pointer on SimpleVector
   */
  virtual inline SimpleVector* getUPtr() const
  {
    return NULL;
  }

  // --- R ---

  /** \fn  const SimpleVector getR(void) const
   *  \brief get the value of r
   *  \return SimpleVector
   */
  inline const SimpleVector getR() const
  {
    return *r;
  }

  /** \fn SimpleVector* getRPtr(void) const
   *  \brief get r
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getRPtr() const
  {
    return r;
  }

  /** \fn void setR (const SimpleVector& newValue)
   *  \brief set the value of r to newValue
   *  \param SimpleVector newValue
   */
  inline void setR(const SimpleVector& newValue)
  {
    *r = newValue;
  }

  /** \fn void setRPtr(SimpleVector* newPtr)
   *  \brief set R to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setRPtr(SimpleVector *);

  // r memory

  /** \fn  const SiconosMemory getRMemory(void) const
   *  \brief get the value of rMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getRMemory() const
  {
    return *rMemory;
  }

  /** \fn SiconosMemory getRMemoryPtr(void) const
   *  \brief get all the values of the state vector r stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getRMemoryPtr() const
  {
    return rMemory;
  }

  /** \fn void setRMemory(const SiconosMemory &)
   *  \brief set the value of rMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setRMemory(const SiconosMemory& newValue)
  {
    *rMemory = newValue;
  }

  /** \fn void setRMemory(SiconosMemory * newPtr)
   *  \brief set rMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setRMemoryPtr(SiconosMemory *);

  // --- Steps in memory ---

  /** \fn const int getStepsInMemory(void) const
   *  \brief get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline const int getStepsInMemory() const
  {
    return stepsInMemory;
  }

  /** \fn void setStepsInMemory(const int&)
   *  \brief set the value of stepsInMemory
   *  \param int steps : the value to set stepsInMemory
   */
  inline void setStepsInMemory(const int& steps)
  {
    stepsInMemory = steps;
  }

  // --- JacobianX ---

  /** \fn  const SiconosMatrix getJacobianX(void) const
   *  \brief get the value of JacobianX
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianX() const
  {
    return *jacobianX;
  }

  /** \fn SiconosMatrix* getJacobianXPtr(void) const
   *  \brief get JacobianX
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianXPtr() const
  {
    return jacobianX;
  }

  /** \fn void setJacobianX (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianX to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setJacobianX(const SiconosMatrix& newValue)
  {
    *jacobianX = newValue;
  }

  /** \fn void setJacobianXPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianX to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianXPtr(SiconosMatrix *newPtr);

  /** \fn  std::string getVectorFieldFunctionName() const
   *  \brief get name of function that computes vectorField (if vectorField from plugin)
   *  \return a string
   */
  inline const std::string getVectorFieldFunctionName() const
  {
    return vectorFieldFunctionName;
  }

  /** \fn  std::string getComputeJacobianXFunctionName() const
   *  \brief get name of function that computes computeJacobianX (if computeJacobianX from plugin)
   *  \return a string
   */
  inline const std::string getComputeJacobianXFunctionName() const
  {
    return computeJacobianXFunctionName;
  }

  // --- Boundary Conditions ---

  /*\todo: to be finished when BC class will be allright */
  /** \fn  const BoundaryCondition getBoundaryCondition(void) const
   *  \brief get the value of BoundaryCondition
   *  \return an object BoundaryCondition
   */
  //inline BoundaryCondition getBoundaryCondition() const { return *BC; }

  /** \fn BoundaryCondition getBoundaryConditionPtr(void) const
   *  \brief get the BoundaryCondition
   *  \return a pointer on the BoundaryCondition object
   */
  inline BoundaryCondition* getBoundaryConditionPtr() const
  {
    return BC;
  }

  /** \fn void setBoundaryCondition(const BoundaryCondition&)
   *  \brief set the Boundary Conditions
   *  \param ref on an object BoundaryCondition
   */
  //inline void setBoundaryCondition(const BoundaryCondition& newBC) {*BC = newBC; }

  /** \fn void setBoundaryConditionPtr(BoundaryCondition*)
   *  \brief set the BoundaryCondition pointer
   *  \param BoundaryCondition *bc : the BoundaryCondition to set BC
   */
  void setBoundaryConditionPtr(BoundaryCondition *newBC);

  // --- dsxml ---
  /** \fn inline const DSXML* getDynamicalSystemXMLPtr() const
   *  \brief get the object DSXML of the DynamicalSystem
   *  \return a pointer on the DSXML of the DynamicalSystem
   */
  inline const DSXML* getDynamicalSystemXMLPtr() const
  {
    return dsxml;
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
  void setVectorFieldFunction(const std::string & pluginPath, const std::string& functionName);

  /** \fn void setComputeJacobianXFunction(const string&, const string&)
   *  \brief allow to set a specified function to compute jacobianX
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianXFunction(const std::string & pluginPath, const std::string & functionName);

  // --- type of DS ---
  /** \fn inline string getType()
   *  \brief get the type of a DynamicalSystem
   *  \return string : the type of the DynamicalSystem
   */
  inline const std::string  getType() const
  {
    return DSType;
  }

  /** \fn inline string setType()
  *  \brief set the type of a DynamicalSystem
  *  \param string : the type of the DynamicalSystem
  */
  inline void getType(const std::string newType)
  {
    DSType = newType;
  }

  // --- DS input-output ---
  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the DynamicalSystem
   *  \return the vector of DSInputOutput
   */
  inline std::vector<DSInputOutput*> getDSInputOutputs(void)
  {
    return dsioVector;
  }

  /** \fn DSInputOutput* getDSInputOutput(int)
   *  \brief allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
   *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the DynamicalSystem
   *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
   */
  DSInputOutput* getDSInputOutput(const unsigned int&);

  /** \fn void setDSInputOutputs(vector<DSInputOutput*>)
   *  \brief allows to set all the DSInputOutputs of the DynamicalSystem
   *  \param vector<DSInputOutput*> : the vector to set
   */
  inline void setDSInputOutputs(std::vector<DSInputOutput*> newDsioVect)
  {
    dsioVector = newDsioVect;
  }

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the DynamicalSystem
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  void addDSInputOutput(DSInputOutput* dsio)
  {
    dsioVector.push_back(dsio);
  }

  /** \fn  std::map<std::string , SimpleVector*> getTmpWorkVector()
   *  \brief get the vector of temporary saved vector
   *  \return a std vector
   */
  inline std::map<const std::string , SimpleVector*> getTmpWorkVector()
  {
    return tmpWorkVector;
  }

  /** \fn  SimpleVector getTmpWorkVector(const std::string& id)
   *  \brief get a temporary saved vector, ref by id
   *  \return a std vector
   */
  inline SimpleVector* getTmpWorkVector(const std::string & id)
  {
    return tmpWorkVector[id];
  }

  /** \fn void set(map<std::string , SimpleVector*>)
   *  \brief set TmpWorkVector
   *  \param a map<std::string , SimpleVector*>
   */
  inline void setTmpWorkVector(std::map<const std::string , SimpleVector*> newVect)
  {
    tmpWorkVector = newVect;
  }

  /** \fn void addTmpWorkVector(SimpleVector*, const string&)
  *  \brief to add a temporary vector
  *  \param a SimpleVector*
  *  \param a string id
  */
  void addTmpWorkVector(SimpleVector* newVal, const std::string& id)
  {
    *tmpWorkVector[id] = *newVal;
  }

  /** \fn void allocateTmpWorkVector(const std::string&, const int&)
   *  \brief to allocate memory for a new vector in tmp map
   *  \param the id of the SimpleVector
   *  \param an int to set the size
   */
  void allocateTmpWorkVector(const std::string& id, const int& size)
  {
    tmpWorkVector[id] = new SimpleVector(size);
  }

  /** \fn freeTmpWorkVector(const std::string& )
   *  \brief to free memory in the map
   *  \param the id of the SimpleVector to free
   */
  void freeTmpWorkVector(const std::string& id)
  {
    delete tmpWorkVector[id];
  }

  // --- ---

  /** \fn void initMemory(const int& steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory
   */
  virtual void initMemory(const unsigned int&) ;

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, xDot and r in the stored previous values
   *  xMemory, xDotMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  virtual void swapInMemory();

  /** \fn void vectorField (const double& time)
   * \brief Default function for computing the vector field \f$ f: (x,t) \in R^{n} \times R  \mapsto  R^{n}\f$
   * \param double time : the time for the computation
   *  \exception RuntimeException
   */
  virtual void computeVectorField(const double&);

  /** \fn static void computeJacobianX (const double& time)
   *  \brief Default function for computing the gradient of the vector field with the respect
   *  to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : the time for the computation
   *  \exception RuntimeException
   */
  virtual void computeJacobianX(const double&);

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
  virtual void display() const;

  /** \var typedef void (*vfPtr) (int* sizeOfX, double* time, double* xPtr, double* xdotPtr);
   *  \brief signature of plugin function computing the vectorfield
   *  \param int* sizeOfX : the size of the vector X
   *  \param double* time : the time for the computation
   *  \param double* xPtr : the pointer to the first element of the vector X
   *  \param double* jacobianXPtr : the pointer to the first element of the matrix jacobianX (in-out parameter)
   */
  typedef void (*vfPtr)(unsigned int* sizeOfX, const double* time, double* xPtr, double* xdotPtr);

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
  virtual double dsConvergenceIndicator()  ;

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
      LagrangianLinearTIDS (LTIDS), LinearDS (LDS)*/
  std::string  DSType;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** this number defines in a single way the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc.)*/
  std::string  id;

  /** the dimension of the system (i.e. size of the state vector x, or the vector r, ...)*/
  unsigned int n;

  /** initial state of the system */
  SiconosVector *x0;

  /** state of the system, \f$  x \in R^{n}\f$ */
  SiconosVector *x;

  /** the  previous state vectors stored in memory*/
  SiconosMemory *xMemory;

  /** the time derivative of the state x (the velocity) */
  SiconosVector *xDot;

  /** the  previous xDot vectors */
  SiconosMemory *xDotMemory;

  /** the  free state vector (state vector for r=0) */
  SiconosVector *xFree;

  /** the  input vector due to the non-smooth law \f$  r \in R^{n}\f$ (multiplier, force, ...)*/
  SimpleVector *r;

  /**  the previous r vectors */
  SiconosMemory *rMemory;

  /** number of previous states stored in memory */
  unsigned int stepsInMemory;

  /** A container of vectors to save temporary values (for Newton convergence computation for example)*/
  std::map<const std::string, SimpleVector*> tmpWorkVector;

  /** Gradient of the vectorfield \f$ f(x,t) \f$ with respect to \f$ x\f$*/
  SiconosMatrix *jacobianX;

  /* the name of the plugin used to compute the vectorField */
  std::string  vectorFieldFunctionName;

  /* the name of the plugin used to compute the JacobianX */
  std::string  computeJacobianXFunctionName;

  /** boundary conditions defined if the DynamicalSystem has some */
  BoundaryCondition *BC;

  /** the XML object linked to the DynamicalSystem  */
  DSXML *dsxml;

  /** class for plugin managing (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** \fn void (*vectorFieldPtr) (int* sizeOfX, double* time, double* xPtr, double* xdotPtr)
   *  \brief pointer on function to compute vectorfield
   *  \param int* sizeOfX : the size of the vector x
   *  \param double* time : current time
   *  \param double* xPtr : the pointer to the first element of the vector x
   *  \param double* xdotPtr : the pointer to the first element of the vector xDot
   */
  void (*vectorFieldPtr)(unsigned int* sizeOfX, const double* time, double* xPtr, double* xdotPtr);

  /** vector of the DS Inputs-Outputs of the Dynamical System */
  std::vector<DSInputOutput*> dsioVector;

  /** \fn void (*computeJacobianXPtr) (int* sizeOfX, double* time, double* xPtr, double* jacobianXPtr)
   *  \brief  Pointer on function to compute the gradient of the vector field with the respect to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param int* sizeOfX : size of vector X
   *  \param double* time : time for computation
   *  \param double* xPtr : pointer to the first element of X
   *  \param double* jacobianXPtr : pointer to the first element of jacobianX matrix (in-out parameter)
   */
  void (*computeJacobianXPtr)(unsigned int* sizeOfX, const double* time, double* xPtr, double* jacobianXPtr);

  /** Flags to know if pointers have been allocated inside constructors or not */

  bool isX0AllocatedIn;
  bool isXAllocatedIn;
  bool isXMemoryAllocatedIn;
  bool isXDotAllocatedIn;
  bool isXDotMemoryAllocatedIn;
  bool isXFreeAllocatedIn;
  bool isRAllocatedIn;
  bool isRMemoryAllocatedIn;
  bool isJacobianXAllocatedIn;
  bool isBCAllocatedIn;
  std::vector<bool> isDsioAllocatedIn;
};

#endif // DYNAMICALSYSTEM_H



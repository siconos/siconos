#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosConst.h"
#include "RuntimeException.h"
#include "check.h"

#include "SiconosMatrix.h"
#include "SiconosVector.h"
#include "SiconosMemory.h"
#include "SiconosSharedLibrary.h"

#include "NonSmoothDynamicalSystem.h"
#include "DSInputOutput.h"
#include "BoundaryCondition.h"
#include "DynamicalSystemXML.h"

#include <string>
#include <vector>
#include <iostream>

class NonSmoothDynamicalSystem;
class BoundaryCondition;
class DSInputOutput;
class DynamicalSystemXML;
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
 * \dot x = f(x,t) + T(x) u(x, \dot x, t) + r,
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$ u \in R^{uSize}\f$ a "control" term
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

  // ===== CONSTRUCTORS =====

  /** \fn DynamicalSystem(DynamicalSystemXML * nsdsXML, NonSmoothDynamicalSystem* =NULL)
   *  \brief xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** \fn DynamicalSystem(DynamicalSystemXML * nsdsXML, const unsigned int& number, const unsigned int& n,
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

  // ===== DESTRUCTOR =====

  virtual ~DynamicalSystem();

  // ===== GETTERS/SETTERS =====

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
  inline void setType(const std::string newType)
  {
    DSType = newType;
  }

  // --- NonSmoothDynamicalSystem ---

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

  // --- Number ---

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

  // --- Id ---

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

  // --- n ---

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
  void setX0(const SiconosVector&);

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
  void setX(const SiconosVector&);

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
  void setXMemory(const SiconosMemory&);

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
  void setXDot(const SiconosVector&);

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
  void setXDotMemory(const SiconosMemory&);

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
  void setXFree(const SiconosVector&);

  /** \fn void setXFreePtr(SiconosVector* newPtr)
   *  \brief set xFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXFreePtr(SiconosVector *);

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
  void setR(const SimpleVector&);

  /** \fn void setRPtr(SimpleVector* newPtr)
   *  \brief set R to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setRPtr(SimpleVector *);

  // rMemory

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
  void setRMemory(const SiconosMemory&);

  /** \fn void setRMemory(SiconosMemory * newPtr)
   *  \brief set rMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setRMemoryPtr(SiconosMemory *);

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
  void setJacobianX(const SiconosMatrix&);

  /** \fn void setJacobianXPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianX to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianXPtr(SiconosMatrix *newPtr);

  // uSize

  /** \fn const int getUSize(void) const;
   *  \brief to get uSize, size of u
   *  \return the value of uSize
   */
  inline const unsigned int getUSize(void) const
  {
    return uSize;
  }

  /** \fn void setUSize(const unsigned int&)
   *  \brief to set the value of uSize
   *  \param an integer to set the value of uSize
   */
  void setUSize(const unsigned int&);

  // ---  U ---

  /** \fn  const SimpleVector getU(void) const
   *  \brief get the value of u, control term
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getU() const
  {
    return *u;
  }

  /** \fn SiconosVector* getUPtr(void) const
   *  \brief get u, the "control" term
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getUPtr() const
  {
    return u;
  }

  /** \fn void setU (const SiconosVector& newValue)
   *  \brief set the value of u to newValue
   *  \param SiconosVector newValue
   */
  void setU(const SiconosVector&);

  /** \fn void setUPtr(SiconosVector* newPtr)
   *  \brief set u to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setUPtr(SiconosVector *);

  // --- T ---

  /** \fn  const SiconosMatrix getT(void) const
   *  \brief get the value of T
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getT() const
  {
    return *T;
  }

  /** \fn SiconosMatrix* getTPtr(void) const
   *  \brief get T
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getTPtr() const
  {
    return T;
  }

  /** \fn void setT (const SiconosMatrix& newValue)
   *  \brief set the value of T to newValue
   *  \param SiconosMatrix newValue
   */
  void setT(const SiconosMatrix&);

  /** \fn void setTPtr(SiconosMatrix* newPtr)
   *  \brief set T to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setTPtr(SiconosMatrix *newPtr);

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

  /** \fn inline const DynamicalSystemXML* getDynamicalSystemXMLPtr() const
   *  \brief get the object DynamicalSystemXML of the DynamicalSystem
   *  \return a pointer on the DynamicalSystemXML of the DynamicalSystem
   */
  inline const DynamicalSystemXML* getDynamicalSystemXMLPtr() const
  {
    return dsxml;
  }

  /** \fn inline void setDynamicalSystemXMLPtr(DynamicalSystemXML *dsxml)
   *  \brief set the DynamicalSystemXML of the DynamicalSystem
   *  \param DynamicalSystemXML* dsxml : the address of theDynamicalSystemXML to set
   */
  inline void setDynamicalSystemXMLPtr(DynamicalSystemXML *newDsxml)
  {
    dsxml = newDsxml;
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

  // ===== TMP WORK VECTOR =====

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

  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** \fn void initMemory(const int& steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory
   */
  virtual void initMemory(const unsigned int&) ;

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, xDot and r in the stored previous values
   *  xMemory, xDotMemory, rMemory,
   */
  virtual void swapInMemory();

  // ===== COMPUTE PLUGINS FUNCTIONS =====

  // --- getters for plugin functions names ---

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

  /** \fn  std::string getComputeUFunctionName() const
   *  \brief get name of function that computes u (if u from plugin)
   *  \return a string
   */
  inline const std::string getComputeUFunctionName() const
  {
    return computeUFunctionName;
  }

  /** \fn  std::string getComputeTFunctionName() const
   *  \brief get name of function that computes T (if T from plugin)
   *  \return a string
   */
  inline const std::string getComputeTFunctionName() const
  {
    return computeTFunctionName;
  }

  // --- setters for functions to compute plugins ---

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

  /** \fn void setComputeUFunction(const string&, const string&)
   *  \brief allow to set a specified function to compute u
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeTFunction(const string&, const string&)
   *  \brief allow to set a specified function to compute T
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeTFunction(const std::string & pluginPath, const std::string & functionName);

  // --- compute plugin functions ---

  /** \fn void vectorField (const double& time)
   * \brief Default function to compute the vector field \f$ f: (x,t) \in R^{n} \times R  \mapsto  R^{n}\f$
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeVectorField(const double&);

  /** \fn static void computeJacobianX (const double& time)
   *  \brief Default function to compute the gradient of the vector field with the respect
   *  to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianX(const double&);

  /** \fn static void computeU (const double&)
   *  \brief Default function to compute u
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeU(const double&);

  /** \fn static void computeT ()
   *  \brief Default function to compute T
   *  \exception RuntimeException
   */
  virtual void computeT();

  // --- isPlugin ---

  /** \fn const std::vector<bool> getIsPlugin() const;
   *  \brief get flag that checks if u and T are loaded from plugin or not
   *  \return a vector of bool
   */
  inline const std::vector<bool> getIsPlugin() const
  {
    return isPlugin;
  }

  // ===== XML MANAGEMENT FUNCTIONS =====

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

  // ===== MISCELLANEOUS ====

  /** \fn void display()
   *  \brief print the data of the dynamical system on the standard output
   */
  virtual void display() const;

  /** \var typedef void (*vfPtr) (int* sizeOfX, double* time, double* xPtr, double* xdotPtr);
   *  \brief signature of plugin function computing the vectorfield
   *  \param int* sizeOfX : the size of the vector X
   *  \param double* time : current time
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

  /** \fn void dsConvergenceIndicator()
   *  \brief Default function for computing an indicator of convergence
   *  \return a double when DS is a Lagrangian
   */
  virtual double dsConvergenceIndicator();

protected:

  // Default constructor
  DynamicalSystem();

  /** \fn void fillBoundaryConditionsFromXml()
   *  \brief uses the DynamicalSystemXML of the DynamicalSystem to fill BoundaryCondition fields
   *  \exception RuntimeException
   */
  virtual void fillBoundaryConditionsFromXml();

  /** \fn void fillDsioFromXml()
   *  \brief uses the DynamicalSystemXML of the DynamicalSystem to fill DSIO vector
   *  \exception RuntimeException
   */
  virtual void fillDsioFromXml();

  // ===== PROTECTED MEMBERS =====

  /** Dynamical System type: General Dynamical System (NLDS) LagrangianDS (LNLDS),
      LagrangianLinearTIDS (LTIDS), LinearDS (LDS)*/
  std::string  DSType;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** this number defines in a single way the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc.)*/
  std::string  id;

  /** the dimension of the system (i.e. size of the state vector x)*/
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

  /** Gradient of the vectorfield \f$ f(x,t) \f$ with respect to \f$ x\f$*/
  SiconosMatrix *jacobianX;

  /** size of vector u */
  unsigned int uSize;

  /** "control" term */
  SiconosVector *u;

  /** Matrix coefficient of u */
  SiconosMatrix *T;

  /** number of previous states stored in memory */
  unsigned int stepsInMemory;

  /** A container of vectors to save temporary values (for Newton convergence computation for example)*/
  std::map<const std::string, SimpleVector*> tmpWorkVector;

  /** boundary conditions defined if the DynamicalSystem has some */
  BoundaryCondition *BC;

  /** the XML object linked to the DynamicalSystem  */
  DynamicalSystemXML *dsxml;

  /** vector of the DS Inputs-Outputs of the Dynamical System */
  std::vector<DSInputOutput*> dsioVector;

  // --- plugins ---

  /** class for plugin managing (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /* the name of the plugin used to compute the vectorField */
  std::string  vectorFieldFunctionName;

  /* the name of the plugin used to compute the JacobianX */
  std::string  computeJacobianXFunctionName;

  /* the name of the plugin used to compute u */
  std::string  computeUFunctionName;

  /* the name of the plugin used to compute T */
  std::string  computeTFunctionName;

  /** Flag to check if u and T are plugins or not */
  std::vector<bool> isPlugin;

  /** \fn void (*vectorFieldPtr) (int* sizeOfX, double* time, double* xPtr, double* xdotPtr)
   *  \brief pointer on function to compute vectorfield
   *  \param int* sizeOfX : the size of the vector x
   *  \param double* time : current time
   *  \param double* xPtr : the pointer to the first element of the vector x
   *  \param double* xdotPtr : the pointer to the first element of the vector xDot
   */
  void (*vectorFieldPtr)(unsigned int* sizeOfX, const double* time, double* xPtr, double* xdotPtr);

  /** \fn void (*computeJacobianXPtr) (int* sizeOfX, double* time, double* xPtr, double* jacobianXPtr)
   *  \brief  Pointer on function to compute the gradient of the vector field with the respect to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param int* sizeOfX : size of vector x
   *  \param double* time : current time
   *  \param double* xPtr : pointer to the first element of x
   *  \param double* jacobianXPtr : pointer to the first element of jacobianX matrix (in-out parameter)
   */
  void (*computeJacobianXPtr)(unsigned int* sizeOfX, const double* time, double* xPtr, double* jacobianXPtr);

  /** \fn void (*computeUPtr) (int* sizeOfX, double* time, double* xPtr, double* xDotPtr, double* TPtr)
   *  \brief  Pointer on function to compute u
   *  \param double* time : current time
   *  \param int* sizeOfU : size of vector u
   *  \param int* sizeOfX : size of vector x
   *  \param double* xPtr : pointer to the first element of x
   *  \param double* xDotPtr : pointer to the first element of xDot
   *  \param double* TPtr : pointer to the first element of u vector (in-out parameter)
   */
  void (*computeUPtr)(unsigned int* sizeOfU, unsigned int* sizeOfX, const double* time, double* xPtr, double* xDotPtr, double* UPtr);

  /** \fn void (*computeTPtr) (int* sizeOfX, double* xPtr, double* TPtr)
   *  \brief  Pointer on function to compute T
   *  \param int* sizeOfX : size of vector X
   *  \param double* xPtr : pointer to the first element of X
   *  \param double* TPtr : pointer to the first element of T matrix (in-out parameter)
   */
  void (*computeTPtr)(unsigned int* sizeOfU, unsigned int* sizeOfX, double* xPtr, double* TPtr);

  /** Flags to know if pointers have been allocated inside constructors or not */

  /** index corresponds to the following order:  x0, x, xMemory, xDot, xDotMemory, xFree, jacobianX */
  std::vector<bool> isXAllocatedIn;
  /** r and rMemory */
  std::vector<bool> isRAllocatedIn;
  /** u and T */
  std::vector<bool> isControlAllocatedIn;
  /** Boundary conditions */
  bool isBCAllocatedIn;
  /** dsio */
  std::vector<bool> isDsioAllocatedIn;
};

#endif // DYNAMICALSYSTEM_H



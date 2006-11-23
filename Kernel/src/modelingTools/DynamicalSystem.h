/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
 * Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 * Siconos is a free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * Siconos is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siconos; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
 */

/*! \file DynamicalSystem.h

*/

#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosConst.h"
#include "RuntimeException.h"
#include "check.h"

#include "SimpleMatrix.h"
#include "BlockMatrix.h"
#include "SimpleVector.h"
#include "SiconosMemory.h"
#include "SiconosSharedLibrary.h"

#include "NonSmoothDynamicalSystem.h"
#include "DynamicalSystemXML.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>

const std::string LNLDS = "LagrangianDS";
const std::string LTIDS = "LagrangianLinearTIDS";
const std::string LDS = "LinearDS";
const std::string LITIDS = "LinearTIDS";
const std::string NLDS = "NonLinearDS";

class NonSmoothDynamicalSystem;
class DynamicalSystemXML;
class SiconosVector;
class SimpleMatrix;
class BlockMatrix;
class SiconosMemory;
class SiconosSharedLibrary;

//!  General First Order Non Linear Dynamical Systems
/**
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
 *  \date (Creation) April 29, 2004
 *
 * This class defines and computes a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * \dot x = f(x,t) + T(x) u(x,t) + r,
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state.
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$ u \in R^{uSize}\f$ a "control" term.
 *
 *  with \f$ f : R^{n} \times R  \mapsto  R^{n}   \f$ .
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  \f[
 *  x(t_0)=x_0
 * \f]
 * To define a boundary Value Problem, the pointer on  a BoundaryCondition must be set.
 *
 * \f$ f(x,t) \f$ is a plug-in function, and can be computed using computeF(t).
 * Its Jacobian according to x is denoted jacobianXF, and computed thanks to computeJacobianXF(t).
 * f and jacobianXF can be plugged to external functions thanks to setComputeFFunction/setComputeJacobianXFFunction.
 *
 * Right-hand side of the equation is denoted rhs and is computed thanks to computeRhs(t).
 *
 * \f[
 *    rhs(x,t) =  f(x,t) + T(x) u(x,t) + r
 * \f]
 *
 * Its Jacobian according to x is jacobianXRhs:
 *
 *  \f[
 *   jacobianXRhs = \nabla_xrhs(x,t) = \nabla_xf(x,t) + \nabla_xT(x)U(x,t) + \nabla_x r
 *  \f]
 *
 * At the time, \f$\nabla_xT(x)\f$ is never taken into account.
 * Moreover, in order to avoid useless memory allocation,
 * we try to re-use the same memory location for f and rhs (and for jacobianXRhs and jacobianXF).
 *
 * \todo Add a pointer to an object Constraint .
 */
class DynamicalSystem
{
protected:

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

  /** the right hand side of the equation */
  SiconosVector *rhs;

  /** jacobian according to X of this rhs */
  SiconosMatrix * jacobianXRhs;

  /** the  free state vector (state vector for r=0) */
  SiconosVector *xFree;

  /** the  input vector due to the non-smooth law \f$  r \in R^{n}\f$ (multiplier, force, ...)*/
  SiconosVector *r;

  /**  the previous r vectors */
  SiconosMemory *rMemory;

  /** f(x,t) */
  SiconosVector *f;

  /** Gradient of \f$ f(x,t) \f$ with respect to \f$ x\f$*/
  SiconosMatrix *jacobianXF;

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

  /** the XML object linked to the DynamicalSystem  */
  DynamicalSystemXML *dsxml;

  /** Parameters list, last argument of plug-in functions. What are those parameters depends on user´s choice.
   *  This a list of pointer to SimpleVector. Each one is identified thanks to a key which is the plug-in name.
   * A flag is also added in the isAllocatedIn map to check inside-class memory allocation for this object.*/
  std::map<std::string, SimpleVector*> parametersList;

  // --- plugins ---

  /** class for plugin managing (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /* the name of the plugin used to compute f(x,t) */
  std::string  computeFFunctionName;

  /* the name of the plugin used to compute the jacobian of f according to x */
  std::string  computeJacobianXFFunctionName;

  /* the name of the plugin used to compute u */
  std::string  computeUFunctionName;

  /* the name of the plugin used to compute T */
  std::string  computeTFunctionName;

  /** Flag to check if operators are plugged or not (and thus constant)
   * For example isPlugin["jacobianXF"] = false, means that jacobianXF is constant,
   * then computeJacobianXF does not change its value, and not plugged.*/
  std::map<std::string, bool> isPlugin;

  /** DynamicalSystem plug-in to compute f(x,t) - id="f".
   *  @param  : the size of the vector x
   *  @param  : current time
   *  @param  : the pointer to the first element of the vector x
   *  @param  : the pointer to the first element of the vector f(x,t)
   *  @param  : a vector of user-defined parameters
   */
  void (*computeFPtr)(unsigned int, double, const double*, double*, double*);

  /** DynamicalSystem plug-in to compute the gradient of f(x,t) with respect to the state: \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   * @param sizeOfX : size of vector x
   * @param time : current time
   * @param x : pointer to the first element of x
   * @param[in,out] jacob : pointer to the first element of jacobianXF matrix
   * @param[in,out] param   : a vector of user-defined parameters
   */
  void (*computeJacobianXFPtr)(unsigned int, double, const double*, double*, double*);

  /** DynamicalSystem plug-in to compute u(x,t)
   * @param unsigned int sizeOfU : size of vector u
   * @param unsigned int sizeOfX : size of vector x
   * @param double time : current time
   * @param double* x : pointer to the first element of x
   * @param[in,out] double* u : pointer to the first element of u vector (in-out parameter)
   * @param[in,out] double* param   : a vector of user-defined parameters
   */
  void (*computeUPtr)(unsigned int, unsigned int, double, const double*, double*, double*);

  /** DynamicalSystem plug-in to compute T(x)
   * @param unsigned int sizeOfU : size of vector u
   * @param unsigned int sizeOfX : size of vector X
   * @param double* x : pointer to the first element of X
   * @param[in,out] double* T: pointer to the first element of T matrix
   * @param[in,out] double* param   : a vector of user-defined parameters
   */
  void (*computeTPtr)(unsigned int, unsigned int, const double*, double*, double*);

  /** Flags to know if pointers have been allocated inside constructors or not */

  /** Flags to know if pointers have been allocated inside constructors or not */
  std::map<std::string, bool> isAllocatedIn;

  /** set all allocation flags (isAllocated map)
   *  \param bool: = if true (default) set default configuration, else set all to false
   */
  virtual void initAllocationFlags(const bool  = true);

  /** set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  virtual void initPluginFlags(const bool);

  /** init parameter vector corresponding to id to a SimpleVector* of size 1
   *  \param a string, id of the plug-in
   */
  void initParameter(const std::string);

  /** default constructor
   */
  DynamicalSystem();

public:

  // ===== CONSTRUCTORS =====

  /** xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** constructor from a set of data
   *  \param int : reference number for this DynamicalSystem
   *  \param int : dimension of this DynamicalSystem
   *  \param SiconosVector : initial state of this DynamicalSystem
   *  \param string : plugin name for f of this DynamicalSystem (optional)
   *  \param string : plugin name for jacobianXF of this DynamicalSystem (optional)
   *  \exception RuntimeException
   */
  DynamicalSystem(const int, const unsigned int, const SiconosVector&,
                  const std::string = "DefaultPlugin:computeF", const std::string = "DefaultPlugin:computeJacobianXF");

  /** copy constructor
   *  \param a Dynamical system to copy
   */
  DynamicalSystem(const DynamicalSystem &);

  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~DynamicalSystem();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

  /** initialization of xFree
   *  \param a string: the simulation type. For TimeStepping: memory allocation. For EventDriven: links (pointers) to q and velocity.
   */
  virtual void initFreeVectors(const std::string);

  // ===== GETTERS/SETTERS =====

  // --- type of DS ---

  /** get the type of a DynamicalSystem
   *  \return string : the type of the DynamicalSystem
   */
  inline const std::string  getType() const
  {
    return DSType;
  }

  /** set the type of a DynamicalSystem
   *  \param string : the type of the DynamicalSystem
   */
  inline void setType(const std::string newType)
  {
    DSType = newType;
  }

  // --- NonSmoothDynamicalSystem ---

  /** get the NonSmoothDynamicalSystem containing this DynamicalSystem
   *  \return NonSmoothDynamicalSystem*
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return nsds;
  }

  /** set the NonSmoothDynamicalSystem containing the DynamicalSystem
   *  \param NonSmoothDynamicalSystem*
   */
  inline void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newNsds)
  {
    nsds = newNsds;
  }

  // --- Number ---

  /** to get the number of the DynamicalSystem
   *  \return the value of number
   */
  inline const int getNumber() const
  {
    return number;
  }

  /** same as getNumber, but return an unsigned long int, used for set<DynamicalSystem*> in OSI, NSDS ...
   *   as sorting criterion.
   *  \return the value of number
   */
  inline const unsigned long int getNumberForSorting() const
  {
    return number;
  }

  /** allows to set the value of number
   *  \param an integer to set the value of number
   */
  inline void setNumber(const int newNumber)
  {
    number = newNumber;
  }

  // --- Id ---

  /** allows to get the id of the DynamicalSystem
   *  \return the value of ths id
   */
  inline const std::string  getId() const
  {
    return id;
  }

  /** allows to set the value of id
   *  \param a string to set the value of id
   */
  inline void setId(const std::string  newId)
  {
    id = newId;
  }

  // --- n ---

  /** allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline const unsigned int getN(void) const
  {
    return n;
  }

  /** allows to set the value of n
   *  \param an integer to set the value of n
   */
  inline void setN(const unsigned int newN)
  {
    n = newN;
  }

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline const unsigned int getDim(void) const
  {
    return n;
  }

  // --- X0 ---

  /** get the value of x0, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getX0() const
  {
    return *x0;
  }

  /** get x0, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getX0Ptr() const
  {
    return x0;
  }

  /** set the value of x0 to newValue
   *  \param SiconosVector newValue
   */
  void setX0(const SiconosVector&);

  /** set x0 to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setX0Ptr(SiconosVector*);

  // --- X ---

  /** get the value of x, the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getX() const
  {
    return *x;
  }

  /** get x, the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXPtr() const
  {
    return x;
  }

  /** set the value of x to newValue
   *  \param SiconosVector newValue
   */
  void setX(const SiconosVector&);

  /** set x to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXPtr(SiconosVector *);

  // X memory

  /** get the value of xMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getXMemory() const
  {
    return *xMemory;
  }

  /** get all the values of the state vector x stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getXMemoryPtr() const
  {
    return xMemory;
  }

  /** set the value of xMemory
   *  \param a ref on a SiconosMemory
   */
  void setXMemory(const SiconosMemory&);

  /** set xMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setXMemoryPtr(SiconosMemory *);

  // ---  Rhs ---

  /** get the value of rhs derivative of the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getRhs() const
  {
    return *rhs;
  }

  /** get rhs, the derivative of the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getRhsPtr() const
  {
    return rhs;
  }

  /** set the value of rhs to newValue
   *  \param SiconosVector newValue
   */
  void setRhs(const SiconosVector&);

  /** set rhs to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setRhsPtr(SiconosVector *);

  // --- JacobianXRhs ---

  /** get the value of JacobianXRhs
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXRhs() const
  {
    return *jacobianXRhs;
  }

  /** get JacobianXRhs
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianXRhsPtr() const
  {
    return jacobianXRhs;
  }

  /** set the value of JacobianXRhs to newValue
   *  \param SiconosMatrix newValue
   */
  virtual void setJacobianXRhs(const SiconosMatrix&);

  /** set JacobianXRhs to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  virtual void setJacobianXRhsPtr(SiconosMatrix *newPtr);

  // --- XFree ---

  /** get the value of xFree
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getXFree() const
  {
    return *xFree;
  }

  /** get xFree
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXFreePtr() const
  {
    return xFree;
  }

  /** set the value of xFree to newValue
   *  \param SiconosVector newValue
   */
  void setXFree(const SiconosVector&);

  /** set xFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXFreePtr(SiconosVector *);

  // --- R ---

  /** get the value of r
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   *  \return a SiconosVector
   */
  inline const SimpleVector getR() const
  {
    return *r;
  }

  /** get r
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getRPtr() const
  {
    return r;
  }

  /** set the value of r to newValue
   *  \param SiconosVector newValue
   */
  void setR(const SiconosVector&);

  /** set R to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setRPtr(SiconosVector *);

  // rMemory

  /** get the value of rMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getRMemory() const
  {
    return *rMemory;
  }

  /** get all the values of the state vector r stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getRMemoryPtr() const
  {
    return rMemory;
  }

  /** set the value of rMemory
   *  \param a ref on a SiconosMemory
   */
  void setRMemory(const SiconosMemory&);

  /** set rMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setRMemoryPtr(SiconosMemory *);

  // ---  F ---

  /** get the value of f derivative of the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getF() const
  {
    return *f;
  }

  /** get f, the derivative of the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getFPtr() const
  {
    return f;
  }

  /** set the value of f to newValue
   *  \param SiconosVector newValue
   */
  void setF(const SiconosVector&);

  /** set f to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setFPtr(SiconosVector *);

  // --- JacobianXF ---

  /** get the value of JacobianXF
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXF() const
  {
    return *jacobianXF;
  }

  /** get JacobianXF
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianXFPtr() const
  {
    return jacobianXF;
  }

  /** set the value of JacobianXF to newValue
   *  \param SiconosMatrix newValue
   */
  virtual void setJacobianXF(const SiconosMatrix&);

  /** set JacobianXF to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  virtual void setJacobianXFPtr(SiconosMatrix *newPtr);

  // uSize

  /** to get uSize, size of u
   *  \return the value of uSize
   */
  inline const unsigned int getUSize(void) const
  {
    return uSize;
  }

  /** to set the value of uSize
   *  \param an integer to set the value of uSize
   */
  void setUSize(const unsigned int);

  // ---  U ---

  /** get the value of u, control term
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getU() const
  {
    return *u;
  }

  /** get u, the "control" term
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getUPtr() const
  {
    return u;
  }

  /** set the value of u to newValue
   *  \param SiconosVector newValue
   */
  void setU(const SiconosVector&);

  /** set u to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setUPtr(SiconosVector *);

  // --- T ---

  /** get the value of T
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getT() const
  {
    return *T;
  }

  /** get T
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getTPtr() const
  {
    return T;
  }

  /** set the value of T to newValue
   *  \param SiconosMatrix newValue
   */
  void setT(const SiconosMatrix&);

  /** set T to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setTPtr(SiconosMatrix *newPtr);

  // --- Steps in memory ---

  /** get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline const int getStepsInMemory() const
  {
    return stepsInMemory;
  }

  /** set the value of stepsInMemory
   *  \param int steps : the value to set stepsInMemory
   */
  inline void setStepsInMemory(const int steps)
  {
    stepsInMemory = steps;
  }

  // --- dsxml ---

  /** get the object DynamicalSystemXML of the DynamicalSystem
   *  \return a pointer on the DynamicalSystemXML of the DynamicalSystem
   */
  inline const DynamicalSystemXML* getDynamicalSystemXMLPtr() const
  {
    return dsxml;
  }

  /** set the DynamicalSystemXML of the DynamicalSystem
   *  \param DynamicalSystemXML* dsxml : the address of theDynamicalSystemXML to set
   */
  inline void setDynamicalSystemXMLPtr(DynamicalSystemXML *newDsxml)
  {
    dsxml = newDsxml;
  }

  // ===== TMP WORK VECTOR =====

  /** get the vector of temporary saved vector
   *  \return a std vector
   */
  inline std::map<const std::string , SimpleVector*> getTmpWorkVector()
  {
    return tmpWorkVector;
  }

  /** get a temporary saved vector, ref by id
   *  \return a std vector
   */
  inline SimpleVector* getTmpWorkVector(const std::string  id)
  {
    return tmpWorkVector[id];
  }

  /** set TmpWorkVector
   *  \param a map<std::string , SimpleVector*>
   */
  inline void setTmpWorkVector(std::map<const std::string , SimpleVector*> newVect)
  {
    tmpWorkVector = newVect;
  }

  /** to add a temporary vector
   *  \param a SimpleVector*
   *  \param a string id
   */
  inline void addTmpWorkVector(SimpleVector* newVal, const std::string id)
  {
    *tmpWorkVector[id] = *newVal;
  }

  /** to allocate memory for a new vector in tmp map
   *  \param the id of the SimpleVector
   *  \param an int to set the size
   */
  inline void allocateTmpWorkVector(const std::string id, const int size)
  {
    tmpWorkVector[id] = new SimpleVector(size);
  }

  /** to free memory in the map
   *  \param the id of the SimpleVector to free
   */
  inline void freeTmpWorkVector(const std::string id)
  {
    delete tmpWorkVector[id];
  }

  /** dynamical system initialization function: mainly set memory and compute value for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const std::string&, double = 0, unsigned int = 1) ;

  /** dynamical system update: mainly call compute for all time or state depending functions
   *  \param current time
   */
  virtual void update(const double);

  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory
   */
  virtual void initMemory(const unsigned int) ;

  /** push the current values of x and r in the stored previous values
   *  xMemory and rMemory,
   */
  virtual void swapInMemory();

  // ===== COMPUTE PLUGINS FUNCTIONS =====

  // --- getters for plugin functions names ---

  /** get name of function that computes f (if f from plugin)
   *  \return a string
   */
  inline const std::string getComputeFFunctionName() const
  {
    return computeFFunctionName;
  }

  /** get name of function that computes computeJacobianXF (if computeJacobianXF from plugin)
   *  \return a string
   */
  inline const std::string getComputeJacobianXFFunctionName() const
  {
    return computeJacobianXFFunctionName;
  }

  /** get name of function that computes u (if u from plugin)
   *  \return a string
   */
  inline const std::string getComputeUFunctionName() const
  {
    return computeUFunctionName;
  }

  /** get name of function that computes T (if T from plugin)
   *  \return a string
   */
  inline const std::string getComputeTFunctionName() const
  {
    return computeTFunctionName;
  }

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute f(x,t)
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeFFunction(const std::string  pluginPath, const std::string functionName);

  /** to set a specified function to compute jacobianXF
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeJacobianXFFunction(const std::string  pluginPath, const std::string  functionName);

  /** to set a specified function to compute u
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string  pluginPath, const std::string  functionName);

  /** to set a specified function to compute T
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeTFunction(const std::string  pluginPath, const std::string  functionName);

  // -- parametersList --

  /** get the full map of parameters
   *  \return a map<string,SimpleVector*>
   */
  inline std::map<std::string, SimpleVector*> getParameters() const
  {
    return parametersList;
  };

  /** get the vector of parameters corresponding to plug-in function named id
   *  \return a SimpleVector
   */
  inline const SimpleVector getParameter(const std::string id)
  {
    return *(parametersList[id]);
  };

  /** get the pointer to the vector of parameters corresponding to plug-in function named id
   *  \return a pointer on a SimpleVector
   */
  inline SimpleVector* getParameterPtr(const std::string id)
  {
    return parametersList[id];
  };

  /** set the map for parameters
   *  \param a map<string, SimpleVector*>
   */
  void setParameters(const std::map<std::string, SimpleVector*>&);

  /** set vector corresponding to plug-in function named id to newValue
   *  \param a SimpleVector
   *  \param a string
   */
  void setParameter(const SimpleVector&, const std::string);

  /** set vector corresponding to plug-in function named id to newPtr (!! pointer link !!)
   *  \param a pointer to SimpleVector
   *  \param a string
   */
  void setParameterPtr(SimpleVector *, const std::string);

  // --- compute plugin functions ---

  /** Default function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeF(const double);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXF(const double, const bool  = false);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeRhs(const double, const bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXRhs(const double, const bool  = false);

  /** Default function to compute u
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeU(const double);

  /** function to compute u when x is not those of the current object.
   *  \param double time : current time
   *  \param SiconosVector* : pointer to a x value
   *  \exception RuntimeException
   */
  virtual void computeU(const double,  SiconosVector* xx);

  /** Default function to compute T
   *  \exception RuntimeException
   */
  virtual void computeT();

  // --- isPlugin ---

  /** get isPlugin, map of flags to check if operators are plugged or not
   *  \return a map of bool
   */
  inline const std::map<std::string, bool> getIsPlugin() const
  {
    return isPlugin;
  }

  /** return true if "name" is plugged, else false (ie name is constant)
   *  \return a map of bool
   */
  inline const bool isPlugged(const std::string name)
  {
    return isPlugin[name];
  }

  // ===== XML MANAGEMENT FUNCTIONS =====

  /** copy the data of the DS in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** copy the data common to each system in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSDataToXML();

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const;

  /** Default function for computing an indicator of convergence
   *  \return a double when DS is a Lagrangian
   */
  virtual double dsConvergenceIndicator();

  /** set R to zero
   */
  virtual void resetNonSmoothPart();


};

#endif // DYNAMICALSYSTEM_H



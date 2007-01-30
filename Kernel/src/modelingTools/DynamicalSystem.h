/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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

Abstract class - General interface for all Dynamical Systems.
*/

#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosConst.h"
#include "RuntimeException.h"
#include "SimpleMatrix.h"
#include "SimpleVector.h"
#include "SiconosMemory.h"
#include "SiconosSharedLibrary.h"
#include "NonSmoothDynamicalSystem.h"
#include "DynamicalSystemXML.h"
#include <string>
#include <vector>
#include <iostream>
#include <map>

class NonSmoothDynamicalSystem;
class DynamicalSystemXML;
class SiconosVector;
class SimpleMatrix;
class SiconosMemory;
class SiconosSharedLibrary;

/** Names used to identify Dynamical Systems */
/** First Order Non Linear DS */
const std::string FONLDS = "FirstOrderNonLinearDS";
/** First Order Linear DS */
const std::string FOLDS = "FirstOrderLinearDS";
/** First Order Linear and Time-Invariant Coefficients DS */
const std::string FOLTIDS = "FirstOrderLinearDS";
/** Lagrangian, Second Order,  Non Linear DS */
const std::string LNLDS = "LagrangianDS";
/** Lagrangian, Second Order,  Linear and Time-Invariant Coefficients DS */
const std::string LLTIDS = "LagrangianLinearTIDS";

/** A map to save temporary working vectors */
typedef std::map<const std::string , SiconosVector*> WorkMap;

/** A map to save temporary working matrices */
typedef std::map<std::string , SiconosMatrix*> WorkMap2;

/** A map to link string to bool (for plug-in flags)  */
typedef std::map<std::string, bool> BoolMap;

/**  General First Order Non Linear Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) April 29, 2004
 *
 * This class defines and computes a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * g(\dot x, x, t, z) = 0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state.
 *    - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *  For example, z may be used to set some perturbation parameters, or to control the system (z will be set by some actuators) or anything else.
 *
 *  with \f$ g : R^{n} \times R  \mapsto  R^{n}   \f$ .
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  \f[
 *  x(t_0)=x_0
 * \f]
 * To define a boundary Value Problem, the pointer on a BoundaryCondition must be set (not yet implemented).
 *
 * If \f[ \nabla_{\dot x} g \f] is invertible, the system can be written as:
 *
 * \f[
 * \dot x = rhs(x, t, z)
 * \f]
 *
 * Right-hand side (\f$ \dot x \f$) of the equation is computed thanks to computeRhs(t).
 *
 * And its Jacobian according to x, named jacobianXRhs, with computeJacobianXRhs(t).
 *
 * Those two functions are pure virtual and must then be implemented in all the derived classes.
 *
 * Dynamical System types (followed by derived classes names):
 *  - First Order Non Linear Dynamical Systems (FirstOrderNonLinearDS)
 *  - First Order Linear DS (FirstOrderLinearDS)
 *  - First Order Linear and Time Invariant Coefficient DS (FirstOrderLinearTIDS)
 *  - Lagrangian DS (LagrangianDS)
 *  - Lagrangian Linear and Time Invariant coefficients DS (LagrangianLinearTIDS)
 *
 * Remarks:
 *  - the copy constructor is declared as a private function => copy and pass-by-value of DS is forbidden.
 *
 */
class DynamicalSystem
{
protected:

  /** Dynamical System type - See possible types in description of this file.*/
  std::string  DSType;

  /** this number defines in a single way the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc.)*/
  std::string  id;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** the dimension of the system (i.e. size of the state vector x) */
  unsigned int n;

  /** initial state of the system */
  SiconosVector *x0;

  /** state of the system, \f$  x \in R^{n}\f$ - Container, with \f$ x[0]=\f$ x \f$ , x[1]= \f$ \dot x \f$ . */
  VectorOfVectors x;

  /** jacobian according to x of the right-hand side (\f$ \dot x = f(x,t) + Tu +r \f$) */
  SiconosMatrix *jacobianXRhs;

  /** Arbitrary algebraic values vector, z. Discret state of the system. */
  SiconosVector * z;

  /** the  previous state vectors stored in memory*/
  SiconosMemory *xMemory;

  /** number of previous states stored in memory */
  unsigned int stepsInMemory;

  /** A container of vectors to save temporary values (for Newton convergence computation for example)*/
  WorkMap workVector;

  /** A container of matrices to save temporary values (zero-matrix, Id-matrix, inverse of Mass or any tmp work matrix ...)
   * No get-set functions at the time. Only used as a protected member.*/
  WorkMap2 workMatrix;

  /** the XML object linked to the DynamicalSystem  */
  DynamicalSystemXML *dsxml;

  // --- plugins utilities---

  /** class for plugin managing (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /** Flag to check if operators are plugged or not (and thus constant)
   * For example isPlugin["jacobianXF"] = false, means that jacobianXF is constant,
   * then computeJacobianXF does not change its value, and not plugged.*/
  BoolMap isPlugin;

  /** Flags to know if pointers have been allocated inside constructors or not */
  AllocationFlagsMap isAllocatedIn;

  /** set all allocation flags (isAllocated map)
   *  \param bool: = if true (default) set default configuration, else set all to false
   */
  virtual void initAllocationFlags(bool  = true);

  // ===== CONSTRUCTORS =====

  /** default constructor
   * \param string: the type of the system, default=FONLDS, non-linear first order system.
   */
  DynamicalSystem(const std::string& = FONLDS);

  /** copy constructor
   *  \param a Dynamical system to copy
   */
  DynamicalSystem(const DynamicalSystem &);

public:

  /** xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** constructor from a set of data
   *  \param string : type of the system
   *  \param int : reference number for this DynamicalSystem
   *  \param int : size of the system (n)
   */
  DynamicalSystem(const std::string&, int, unsigned int newN);

  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~DynamicalSystem();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

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
  inline void setType(const std::string& newType)
  {
    DSType = newType;
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
  inline void setNumber(int newNumber)
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
  inline void setId(const std::string&  newId)
  {
    id = newId;
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
  inline void setN(unsigned int newN)
  {
    n = newN;
  }

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline const unsigned int getDim() const
  {
    return n;
  };

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

  /** get the value of \f$ x \f$, the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getX() const
  {
    return *(x[0]);
  }

  /** get \f$ x \f$ (pointer), the state of the DynamicalSystem.
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXPtr() const
  {
    return x[0];
  }

  /** set the value of \f$ x \f$ (ie (*x)[0]) to newValue
   *  \param SiconosVector newValue
   */
  void setX(const SiconosVector&);

  /** set \f$ x \f$ (ie (*x)[0]) to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXPtr(SiconosVector *);

  // ---  Rhs ---

  /** get the value of the right-hand side, \f$ \dot x \f$, derivative of the state of the DynamicalSystem.
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getRhs() const
  {
    return *(x[1]);
  }

  /** get the right-hand side, \f$ \dot x \f$, the derivative of the state of the DynamicalSystem.
   *  \return a pointer on a SiconosVector
   */
  inline SiconosVector* getRhsPtr() const
  {
    return x[1];
  }

  /** set the value of the right-hand side, \f$ \dot x \f$, to newValue
   *  \param SiconosVector newValue
   */
  void setRhs(const SiconosVector&);

  /** set right-hand side, \f$ \dot x \f$, to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setRhsPtr(SiconosVector *);

  // --- JacobianXRhs ---

  /** get the value of the gradient according to \f$ x \f$ of the right-hand side
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXRhs() const
  {
    return *jacobianXRhs;
  }

  /** get gradient according to \f$ x \f$ of the right-hand side (pointer)
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

  // -- z --

  /** get the value of \f$ z \f$, the vector of algebraic parameters.
   * \return a SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getZ() const
  {
    return *z;
  }

  /** get \f$ z \f$ (pointer), the vector of algebraic parameters.
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getZPtr() const
  {
    return z;
  }

  /** set the value of \f$ z \f$ to newValue
   *  \param SiconosVector newValue
   */
  void setZ(const SiconosVector&);

  /** set \f$ z \f$ to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setZPtr(SiconosVector *);

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
  inline void setStepsInMemory(int steps)
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

  // ===== WORK VECTOR =====

  /** get the vector of temporary saved vector
   *  \return a WorkMap (map that links string to vectors)
   */
  inline WorkMap getWorkVector()
  {
    return workVector;
  }

  /** get a temporary saved vector, ref by id
   *  \return a std vector
   */
  inline SiconosVector* getWorkVector(const std::string&  id)
  {
    return workVector[id];
  }

  /** set WorkVector
   *  \param a WorkMap
   */
  inline void setWorkVector(const WorkMap& newVect)
  {
    workVector = newVect;
  }

  /** to add a temporary vector
   *  \param a SiconosVector*
   *  \param a string id
   */
  inline void addWorkVector(SiconosVector* newVal, const std::string& id)
  {
    *workVector[id] = *newVal;
  }

  /** to allocate memory for a new vector in tmp map
   *  \param the id of the SimpleVector
   *  \param an int to set the size
   */
  inline void allocateWorkVector(const std::string& id, int size)
  {
    workVector[id] = new SimpleVector(size);
  }

  /** to free memory in the map
   *  \param the id of the SimpleVector to free
   */
  inline void freeWorkVector(const std::string& id)
  {
    delete workVector[id];
  }

  /** Initialization function for the rhs and its jacobian (including memory allocation).
   *  \param time of initialization
   */
  virtual void initRhs(double) = 0 ;

  /** dynamical system initialization function: mainly set memory and compute value for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const std::string&, double = 0, unsigned int = 1) = 0;

  /** dynamical system update: mainly call compute for all time or state depending functions
   *  \param current time
   */
  void update(double);

  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** initialize the SiconosMemory objects: reserve memory for i vectors in memory and reset all to zero.
   *  \param the size of the SiconosMemory (i)
   */
  virtual void initMemory(unsigned int) ;

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  virtual void swapInMemory() = 0;

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeRhs(double, bool  = false) = 0;

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXRhs(double, bool  = false) = 0;

  // --- isPlugin ---

  /** get isPlugin, map of flags to check if operators are plugged or not
   *  \return a map of bool
   */
  inline const BoolMap getIsPlugin() const
  {
    return isPlugin;
  }

  /** return true if "name" is plugged, else false (ie name is constant)
   *  \return a map of bool
   */
  inline const bool isPlugged(const std::string& name)
  {
    return isPlugin[name];
  }

  // ===== XML MANAGEMENT FUNCTIONS =====

  /** copy the data of the DS into the XML tree
   */
  virtual void saveDSToXML();

  /** copy the data specific to each system into the XML tree
   */
  virtual void saveSpecificDataToXML() = 0;

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const = 0;

  /** Default function for computing an indicator of convergence
   *  \return a double when DS is a Lagrangian
   */
  virtual double dsConvergenceIndicator();

  /** set R to zero
   */
  virtual void resetNonSmoothPart() = 0;


};

#endif // DYNAMICALSYSTEM_H



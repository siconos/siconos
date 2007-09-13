/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
  \brief Abstract class - General interface for all Dynamical Systems.
*/

#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosSharedLibrary.h"
#include "RuntimeException.h"
#include "Tools.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"
#include "SiconosMemory.h"

class NonSmoothDynamicalSystem;
class DynamicalSystemXML;
class SiconosVector;
class SiconosMatrix;
class SimpleMatrix;
class SimpleVector;
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

/** Pointer to function for plug-in. */
typedef void (*FPtr6)(double, unsigned int, const double*, const double*, double*, unsigned int, double*);

/**  Abstract class to handle Dynamical Systems => interface for derived classes (First Order or Lagrangian systems)
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) January 15, 2007
 *
 * This class is used to describe dynamical systems of the form :
 * \f[
 * g(\dot x, x, t, z) = 0
 * \f]
 * where
 *    - \f$ x \in R^{n} \f$ is the state.
 *    - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *  For example, z may be used to set some perturbation parameters, or to control the system (z will be set by some actuators) or anything else.
 *
 *  with \f$ g : R^{n} \times R  \mapsto  R^{n}   \f$ .
 *
 * Operators and the functions used to compute them:
 *
 * - g: computeG(t)
 * - jacobianXG[0] = \f$ \nabla_x g(t,\dot x,x,z) \f$: computeJacobianG(0,...)
 * - jacobianXG[1] = \f$ \nabla_{\dot x} g(t,\dot x,x,z) \f$: computeJacobianG(1,...)
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  \f[
 *  x(t_0)=x_0
 * \f]
 * To define a boundary Value Problem, a pointer on a BoundaryCondition must be set (not yet implemented).
 *
 * Under some specific conditions, the system can be written as:
 *
 * \f[
 * \dot x = rhs(x, t, z)
 * \f]
 *
 * In that case, \f$ \nabla_{\dot x} g \f$ must be invertible.
 * Right-hand side (\f$ \dot x \f$) of the equation is computed thanks to computeRhs(t).
 *
 * And its Jacobian according to x, named jacobianXRhs, with computeJacobianXRhs(t).
 *
 * <b> Those two functions (computeRhs and computeJacobianXRhs) are pure virtual and must be implemented in all the derived classes. </b>
 *
 * Dynamical System types (followed by derived classes names):
 *  - First Order Non Linear Dynamical Systems (FirstOrderNonLinearDS)
 *  - First Order Linear DS (FirstOrderLinearDS)
 *  - First Order Linear and Time Invariant Coefficient DS (FirstOrderLinearTIDS)
 *  - Lagrangian DS (LagrangianDS)
 *  - Lagrangian Linear and Time Invariant coefficients DS (LagrangianLinearTIDS)
 *
 * About members:
 *  - A DynamicalSystem is identified thanks to a number and an id.
 *  - It handles a pointer to the owner NonSmoothDynamicalSystem and a pointer to a DynamicalSystemXML object for data i/o.
 *  - A VectorOfVectors, x,  is used to saved the state: x[0]=\f$ x \f$ and x[1]=\f$ \dot x \f$ = right-hand side.
 *
 * Remarks:
 *  - the copy constructor is declared as a private function => copy and pass-by-value of DS is forbidden.
 *
 * Warning:
 *  - At the time, nothing is implemented in simulation to proceed with systems written as \f$ g(...) = 0 \f$. Then use only the form
 *    \f$ \dot x = rhs(...) \f$.
 *
 */
class DynamicalSystem
{
protected:

  /** Dynamical System type - See possible types in description of this file.*/
  std::string  DSType;

  /** An id number for the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc)*/
  std::string  id;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** the dimension of the system (\e ie size of the state vector x) */
  unsigned int n;

  /** initial state of the system */
  SiconosVector *x0;

  /** state of the system, \f$  x \in R^{n}\f$ - With \f$ x[0]=\f$ x \f$ , x[1]= \f$ \dot x \f$ . */
  VectorOfVectors x;

  /** jacobian according to x of the right-hand side (\f$ \dot x = f(x,t) + r \f$) */
  SiconosMatrix *jacobianXRhs;

  /** Arbitrary algebraic values vector, z, discret state of the system. */
  SiconosVector * z;

  /** \f$ g(t,\dot x,x,z) \f$ */
  SiconosVector *g;

  /** jacobianXG[0] = \f$ \nabla_x g(t,\dot x,x,z) \f$, jacobianXG[1] = \f$ \nabla_{\dot x} g(t,\dot x,x,z) \f$  */
  VectorOfMatrices jacobianG;

  /** contains the names of the various plug-in. Example: pluginNames["mass"] is the function used to compute the mass.
     Id are the names of the member (mass, fInt ...) */
  NamesList pluginNames;

  /** DynamicalSystem plug-in to compute \f$ g(t,\dot x,x,z) \f$
   *  @param   current time
   *  @param   the size of the vector x
   *  @param   the pointer to the first element of the vector x[0]=\f$ x \f$
   *  @param   the pointer to the first element of the vector x[1]=\f$ \dot x \f$
   *  @param   the pointer to the first element of the vector g(t, ...)
   *  @param   the size of the vector z
   *  @param   a vector of parameters, z
   */
  FPtr6 computeGPtr;

  /** Plug-in to compute jacobianG (computeJacobianGPtr[i] for jacobianG[i]).
   *  @param   current time
   *  @param   the size of the vector x
   *  @param   the pointer to the first element of the vector x[0]=\f$ x \f$
   *  @param   the pointer to the first element of the vector x[1]=\f$ \dot x \f$
   *  @param   the pointer to the first element of the vector g(t, ...)
   *  @param   the size of the vector z
   *  @param   a vector of parameters, z
   */
  std::vector<FPtr6> computeJacobianGPtr;

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

  /** Flag to check if operators are plugged or not (and thus constant).
   * For example isPlugin["jacobianXF"] = false, means that jacobianXF is constant, not plugged ,
   * then computeJacobianXF does not change its value. */
  BoolMap isPlugin;

  /** Flags to know if pointers have been allocated inside constructors or not */
  BoolMap isAllocatedIn;

  /** set all allocation flags (isAllocated map) to true or false
   *  \param a bool
   */
  virtual void initAllocationFlags(bool);

  /** set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  virtual void initPluginFlags(bool);

  // ===== CONSTRUCTORS =====

  /** default constructor
   * \param string: the type of the system, default=FONLDS, non-linear first order system.
   */
  DynamicalSystem(const std::string& = FONLDS);

  /** copy constructor => private, no copy nor pass-by-value.
   */
  DynamicalSystem(const DynamicalSystem &);

public:

  /** xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
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

  /** function used to sort DynamicalSystem in SiconosSet<DynamicalSystem*>
   *  \return an int
   */
  inline const int getSort() const
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
  void setJacobianXRhs(const SiconosMatrix&);

  /** set JacobianXRhs to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianXRhsPtr(SiconosMatrix *newPtr);

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

  // --- g ---

  /** get the value of \f$ g \f$
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getG() const
  {
    return *g;
  }

  /** get \f$ g \f$ (pointer)
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getGPtr() const
  {
    return g;
  }

  /** set the value of \f$ g \f$ to newValue
   *  \param SiconosVector newValue
   */
  void setG(const SiconosVector&);

  /** set \f$ g \f$ to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setGPtr(SiconosVector *);

  // -- Jacobian g --

  /** get the value of jacobianG
   *  \param index (0: \f$ \nabla_x \f$, 1: \f$ \nabla_{\dot x} \f$ )
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianG(unsigned int i) const
  {
    return *jacobianG[i];
  }

  /** get JacobianG
   *  \param index (0: \f$ \nabla_x \f$, 1: \f$ \nabla_{\dot x} \f$ )
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianGPtr(unsigned int i) const
  {
    return jacobianG[i];
  }

  /** set the value of JacobianG to newValue
   *  \param index (0: \f$ \nabla_x \f$, 1: \f$ \nabla_{\dot x} \f$ )
   *  \param SiconosMatrix newValue
   */
  void setJacobianG(unsigned int, const SiconosMatrix&);

  /** set JacobianG to pointer newPtr
   *  \param index (0: \f$ \nabla_x \f$, 1: \f$ \nabla_{\dot x} \f$ )
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianGPtr(unsigned int, SiconosMatrix *newPtr);

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

  /** get name of function that computes "name"
   *  \return a string
   */
  inline const std::string getFunctionName(const std::string& name) const
  {
    return (pluginNames.find(name))->second;
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

  // --- get plugin functions names ---

  /** get the list of plug-in names
   *  \return a NamesList
   */
  inline NamesList getComputeFunctionNames() const
  {
    return pluginNames;
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

  /** to set a specified function to compute g
   *  \param index (0: \f$ \nabla_x \f$, 1: \f$ \nabla_{\dot x} \f$ )
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   */
  void setComputeGFunction(const std::string&  pluginPath, const std::string& functionName);

  /** to set a specified function to compute jacobianG
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   */
  void setComputeJacobianGFunction(unsigned int, const std::string&  pluginPath, const std::string&  functionName);

  /** Default function to compute g
   *  \param double, the current time
   */
  void computeG(double);

  /** default function to compute the gradient of g
   *  \param double time : the current time
   *  \param index (0: \f$ \nabla_x \f$, 1: \f$ \nabla_{\dot x} \f$ )
   */
  void computeJacobianG(unsigned int, double);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeRhs(double, bool  = false) = 0;

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
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



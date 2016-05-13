/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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

/*! \file DynamicalSystem.hpp
  \brief Abstract class - General interface for all Dynamical Systems.
*/

#ifndef DynamicalSystem_H
#define DynamicalSystem_H

#include "SiconosPointers.hpp"
#include "SiconosFwd.hpp"

#include "SSLH.hpp"
#include "RuntimeException.hpp"

#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosMemory.hpp"
#include "DynamicalSystemTypes.hpp"
#include "PluggedObject.hpp"
#include "PluginTypes.hpp"
#include "SiconosVisitor.hpp"

/**  Abstract class to handle Dynamical Systems => interface for
   derived classes (First Order or Lagrangian systems)

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) January 15, 2007

  This class is used to describe dynamical systems of the form :
  \f[
  g(\dot x, x, t, z) = 0
  \f]
  where

     - \f$ x \in R^{n} \f$ is the state.

     - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic
   variables, some sort of discret state.  For example, z may be used
   to set some perturbation parameters, or to control the system (z
   will be set by some actuators) or anything else.

   with \f$ g : R^{n} \times R  \to  R^{n}   \f$ .

  Operators and the functions used to compute them:

  - g: computeg(t)

  - jacobianXG[0] = \f$ \nabla_x g(t,\dot x,x,z) \f$:
    computeJacobiang(0,...)

  - jacobianXG[1] = \f$ \nabla_{\dot x} g(t,\dot x,x,z) \f$:
    computeJacobiang(1,...)

  By default, the DynamicalSystem is considered to be an Initial Value
  Problem (IVP) and the initial conditions are given by

   \f[
   x(t_0)=x_0
  \f]

  Under some specific conditions, the system can be written as:

  \f[
  \dot x = rhs(x, t, z)
  \f]

  In that case, \f$ \nabla_{\dot x} g \f$ must be invertible.
  Right-hand side (\f$ \dot x \f$) of the equation is computed thanks
  to computeRhs(t)

  and its Jacobian according to x, named jacobianRhsx, with
  computeJacobianRhsx(t).

  <b> Those two functions (computeRhs and computeJacobianRhsx) are
  pure virtual and must be implemented in all the derived
  classes. </b>

  Dynamical System types (followed by derived classes names):

   - First Order Non Linear Dynamical Systems (FirstOrderNonLinearDS)

   - First Order Linear DS (FirstOrderLinearDS)

   - First Order Linear and Time Invariant Coefficient DS
     (FirstOrderLinearTIDS)

   - Lagrangian DS (LagrangianDS)

   - Lagrangian Linear and Time Invariant coefficients DS
     (LagrangianLinearTIDS)

  About members:

  - A DynamicalSystem is identified thanks to a number.

   - A VectorOfVectors, x, is used to saved the state: x[0]=\f$ x \f$
     and x[1]=\f$ \dot x \f$ = right-hand side.

   - number is set automatically using count static variable except
   Warning:

   - At the time, nothing is implemented in simulation to proceed with
     systems written as \f$ g(...) = 0 \f$. Then use only the form \f$
     \dot x = rhs(...) \f$.

 */

class DynamicalSystem
{

public:
  /** List of indices used to save tmp work vectors
   * The last value is the size of the present list, so you HAVE to leave it at the end position.
   */
  enum WorkNames {local_buffer, freeresidu, free, qtmp, acce_memory, acce_like, free_tdg, sizeWorkV};

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(DynamicalSystem);

  /** used to set ds number */
  static unsigned int count;

  /** copy constructor => private, no copy nor pass-by-value.
   */
  //  DynamicalSystem(const DynamicalSystem & );

protected:

  /** An id number for the DynamicalSystem */
  unsigned int _number;

  /** the dimension of the system (\e ie size of the state vector x) */
  unsigned int _n;

  /** initial state of the system */
  SP::SiconosVector _x0;

  /** the input vector due to the non-smooth law \f$ r \in R^{n}\f$
   * (multiplier, force, ...)
   * \remark V.A. 17/09/2011 :
   * This should be a VectorOfVectors as for _x when higher relative degree
   * systems will be simulated
   */
  SP::SiconosVector _r;

  /** used by the relative convergence criteron*/
  double _normRef;

  /** state of the system,
   *  \f$  x \in R^{n}\f$ - With x[0]=\f$ x \f$ , x[1]= \f$ \dot{x} \f$ . */
  VectorOfVectors _x;

  /** jacobian according to x of the right-hand side (\f$ \dot x =
      f(x,t) + r \f$) */
  SP::SiconosMatrix _jacxRhs;

  //  SP::SiconosMatrix _jacgx;
  //  SP::SiconosMatrix _jacxDotG;
  //  SP::SiconosMatrix jacobianZG;

  /** Arbitrary algebraic values vector, z, discret state of the
      system. */
  SP::SiconosVector _z;


  /** DynamicalSystem plug-in to compute \f$ g(t,\dot x,x,z) \f$
   *  @param   current time
   *  @param   the size of the vector x
   *  @param   the pointer to the first element of the vector x[0]=\f$ x \f$
   *  @param   the pointer to the first element of the vector x[1]=\f$ \dot x \f$
   *  @param   the pointer to the first element of the vector g(t, ...)
   *  @param   the size of the vector z
   *  @param   a vector of parameters, z
   */
  SP::PluggedObject _pluging;

  /** Plug-in to compute jacobianG (computeJacobiangPtr[i] for jacobianG[i]).
   *  @param   current time
   *  @param   the size of the vector x
   *  @param   the pointer to the first element of the vector x[0]=\f$ x \f$
   *  @param   the pointer to the first element of the vector x[1]=\f$ \dot x \f$
   *  @param   the pointer to the first element of the vector g(t, ...)
   *  @param   the size of the vector z
   *  @param   a vector of parameters, z
   */
  SP::PluggedObject _pluginJacgx;
  SP::PluggedObject _pluginJacxDotG;


  /** the  previous state vectors stored in memory*/
  SP::SiconosMemory _xMemory;

  /** number of previous states stored in memory */
  unsigned int _stepsInMemory;

  /** A container of vectors to save temporary values (for Newton convergence computation for example)*/
  VectorOfVectors _workspace;

  /** A container of matrices to save temporary values (zeroMatrix, idMatrix, inverse of Mass or any tmp work matrix ...)
   * No get-set functions at the time. Only used as a protected member.*/
  VectorOfSMatrices _workMatrix;

  // ===== CONSTRUCTORS =====

  /** default constructors/destructor
   */
  DynamicalSystem();
  virtual void zeroPlugin();

protected:
  /** a vector reserved to compute the freeState.*/
//  SP::SiconosVector _workFree;

public:

  /** constructor from a set of data
   *  \param newN int : size of the system (n)
   */
  DynamicalSystem(unsigned int newN);

  /** Copy constructor
   * \param ds the DynamicalSystem to copy
   */
  DynamicalSystem(const DynamicalSystem & ds);
  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~DynamicalSystem() {}

  //@}

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem() = 0;

  /*! @name Members access */
  //@{

  // --- Number ---

  /** to get the number of the DynamicalSystem
   *  \return the value of number
   */
  inline int number() const
  {
    return _number;
  }

  /** function used to sort DynamicalSystem in SiconosSet<SP::DynamicalSystem>
   *  \return an int (warning: must be const, despite intel compilers warning, because of SiconosSet Cmp function arguments)
   */
  inline int getSort() const
  {
    return _number;
  }

  // --- n ---

  /** allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline unsigned int getN() const
  {
    return _n;
  }

  /** allows to set the value of n
   *  \param newN an integer to set the value of n
   */
  inline void setN(unsigned int newN)
  {
    _n = newN;
  }

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline unsigned int getDim() const
  {
    return _n;
  };

  // --- X0 ---

  /** get the value of x0, the initial state of the DynamicalSystem
   *  \return SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */
  inline const SiconosVector getX0() const
  {
    return *_x0;
  }

  /** get x0, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector x0() const
  {
    return _x0;
  };

  /** get _normRef
   * \return a reference to _normRef
   */
  inline double normRef() const
  {
    return _normRef;
  };

  // --- R ---

  /** get the value of r
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   *  \return a vector
   */
  inline const SiconosVector getR() const
  {
    return *_r;
  }

  /** get r
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector r() const
  {
    return _r;
  }

  /** set the value of r to newValue
   *  \param newValue SiconosVector 
   */
  void setR(const SiconosVector& newValue );

  /** set R to pointer newPtr
   *  \param newPtr SP::SiconosVector newPtr
   */
  void setRPtr(SP::SiconosVector newPtr);

  /** set the value of x0 to newValue
   *  \param newValue SiconosVector newValue
   */
  void setX0(const SiconosVector& newValue);

  /** set x0 to pointer newPtr
   *  \param newPtr SP::SiconosVector newPtr
   */
  void setX0Ptr(SP::SiconosVector newPtr);

  // --- X ---

  /** get the value of \f$ x \f$, the state of the DynamicalSystem
   * \return SiconosVector
   */
  inline const SiconosVector& getx() const
  {
    return *(_x[0]);
  }

  /** get \f$ x \f$ (pointer), the state of the DynamicalSystem.
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector x() const
  {
    return _x[0];
  }

  /** set the value of \f$ x \f$ (ie (*x)[0]) to newValue
   *  \param newValue SiconosVector 
   */
  void setX(const SiconosVector& newValue);

  /** set \f$ x \f$ (ie (*x)[0]) to pointer newPtr
   *  \param newPtr SP::SiconosVector 
   */
  void setXPtr(SP::SiconosVector newPtr);

  // ---  Rhs ---

  /** get the value of the right-hand side, \f$ \dot x \f$, derivative of the state of the DynamicalSystem.
   *  \return SiconosVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */
  inline SiconosVector& getRhs() const
  {
    return *(_x[1]);
  }

  /** get the right-hand side, \f$ \dot x \f$, the derivative of the state of the DynamicalSystem.
   *  \return a pointer on a SiconosVector
   */
  inline SP::SiconosVector rhs() const
  {
    return _x[1];
  }

  /** set the value of the right-hand side, \f$ \dot x \f$, to newValue
   *  \param newValue SiconosVector newValue
   */
  void setRhs(const SiconosVector& newValue);

  /** set right-hand side, \f$ \dot x \f$, to pointer newPtr
   *  \param newPtr SP::SiconosVector newPtr
   */
  void setRhsPtr(SP::SiconosVector newPtr);

  // --- JacobianRhsx ---

  /** get the value of the gradient according to \f$ x \f$ of the right-hand side
   *  \return SimpleMatrix
   */
  inline SiconosMatrix& getJacobianRhsx() const
  {
    return *_jacxRhs;
  }

  /** get gradient according to \f$ x \f$ of the right-hand side (pointer)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianRhsx() const
  {
    return _jacxRhs;
  }

  /** set the value of JacobianRhsx to newValue
   *  \param newValue SiconosMatrix newValue
   */
  void setJacobianRhsx(const SiconosMatrix& newValue);

  /** set JacobianRhsx to pointer newPtr
   *  \param newPtr SP::SiconosMatrix  
   */
  void setJacobianRhsxPtr(SP::SiconosMatrix newPtr);

  // -- z --

  /** get the value of \f$ z \f$, the vector of algebraic parameters.
   * \return a SiconosVector
   */
  inline const SiconosVector& getz() const
  {
    return *_z;
  }

  /** get \f$ z \f$ (pointer), the vector of algebraic parameters.
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector z() const
  {
    return _z;
  }

  /** set the value of \f$ z \f$ to newValue
   *  \param newValue SiconosVector 
   */
  void setz(const SiconosVector& newValue) ;

  /** set \f$ z \f$ to pointer newPtr
   *  \param newPtr SP::SiconosVector 
   */
  void setzPtr(SP::SiconosVector newPtr);

  // X memory

  /** get all the values of the state vector x stored in a SiconosMemory object
   *  \return a pointer to the SiconosMemory object
   */
  inline SP::SiconosMemory xMemory() const
  {
    return _xMemory;
  }

  // --- Steps in memory ---

  /** get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline int getStepsInMemory() const
  {
    return _stepsInMemory;
  }

  /** set the value of stepsInMemory
   *  \param steps  the value to set stepsInMemory
   */
  inline void setStepsInMemory(int steps)
  {
    _stepsInMemory = steps;
  }

  // ===== WORK VECTOR =====

  /** get the vector of temporary saved vector
   *  \return a VectorOfVectors (map that links std::string to vectors)
   */
  inline VectorOfVectors workspace() const
  {
    return _workspace;
  }

  /** get a temporary saved vector, ref by id
   * \param id  WorkNames
   * \return a SP::SiconosVector
   */
  inline SP::SiconosVector workspace(const WorkNames& id) const
  {
    return _workspace[id];
  }

  /** set WorkVector
   *  \param newVect a VectorOfVectors
   */
  inline void setWorkVector(const VectorOfVectors& newVect)
  {
    _workspace = newVect;
  }

  /** to add a temporary vector
   *  \param newVal a SP::SiconosVector
   *  \param id a std::string id
   */
  inline void addWorkVector(SP::SiconosVector newVal, const WorkNames& id)
  {
    *_workspace[id] = *newVal;
  }
  /** sub a vector to a temporary one
   *  \param newVal a SP::SiconosVector
   *  \param id a std::string id
   */
  inline void subWorkVector(SP::SiconosVector newVal, const WorkNames& id)
  {
    *_workspace[id] -= *newVal;
  }

  /** to allocate memory for a new vector in tmp map
   *  \param id the id of the SiconosVector
   *  \param size an int to set the size
   */
  inline void allocateWorkVector(const WorkNames& id, int size)
  {
    _workspace[id].reset(new SiconosVector(size));
  }

  //@}

  /** Determine whether this is a linear DS
   * \return true if the Dynamical system is linear.
   */
  virtual bool isLinear()
  {
    return false;
  }

  /** Initialization function for the rhs and its jacobian (including
   *  memory allocation). 
   * \param time of initialization
   */
  virtual void initRhs(double time) = 0 ;

  /** dynamical system initialization function except for _r :
   *  mainly set memory and compute value for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param size the size of the memory, default size = 1.
   */
  virtual void initialize(double time = 0, unsigned int size = 1) = 0;

  /** dynamical system initialization function for NonSmoothInput _r
   *  \param level of _r.
   */
  virtual void initializeNonSmoothInput(unsigned int level) = 0;

  /** dynamical system update: mainly call compute for all time or state depending functions
   *  \param time  current time
   */
  void update(double time);

  /*! @name Memory vectors management  */
  //@{
  /** initialize the SiconosMemory objects: reserve memory for i vectors in memory and reset all to zero.
   *  \param steps the size of the SiconosMemory (i)
   */
  virtual void initMemory(unsigned int steps);

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  virtual void swapInMemory() = 0;

  //@}
  /*! @name Plugins management  */
  //@{
  /** to set a specified function to compute g
   *  \param  pluginPath std::string pluginPath : the complete path to the plugin
   *  \param functionName std::string functionName : the function name to use in this library
   */
  void setComputegFunction(const std::string&  pluginPath, const std::string& functionName);

  /** set a specified function to compute g
   *  \param fct a pointer on the plugin function
   */
  void setComputegFunction(FPtr6 fct);

  /** to set a specified function to compute jacobianG
   *  \param pluginPath std::string  : the complete path to the plugin
   *  \param functionName the std::string functionName : function name to use in this library
   */
  void setComputeJacobianXGFunction(const std::string&  pluginPath, const std::string&  functionName);

  
 /** to set a specified function to compute jacobianDotXG
   *  \param pluginPath std::string  : the complete path to the plugin
   *  \param functionName the std::string functionName : function name to use in this library
   */
  void setComputeJacobianDotXGFunction(const std::string&  pluginPath, const std::string&  functionName);
  //  void setComputeJacobianZGFunction( const std::string&  pluginPath, const std::string&  functionName);
  
  /** set a specified function to compute jacobianG
   *  \param newPtr a pointer on the plugin function
   */
  virtual void setComputeJacobianXGFunction(FPtr6 newPtr) {};
  
  /** set a specified function to compute jacobianG
   *  \param newPtr a pointer on the plugin function
   */
  virtual void setComputeJacobianDotXGFunction(FPtr6 newPtr) {};
  
  /** Default function to compute g
   *  \param time double, the current time
   */
//  void computeg(double time);

  /**default function to update the plugins functions using a new time:
   * \param time  the current time
   */
  virtual void updatePlugins(double time)
  {
    ;
  }

  //@}


  /*! @name Right-hand side computation */
  //@{
  /** Default function to the right-hand side term
   *  \param time  (double)  current time
   *  \param isDSup (bool)  flag to avoid recomputation of operators
   */
  virtual void computeRhs(double time, bool isDSup = false) = 0;

  /** Default function to jacobian of the right-hand side term according to x
   *  \param time  (double)  current time
   *  \param isDSup (bool)  flag to avoid recomputation of operators
   */
  virtual void computeJacobianRhsx(double time , bool isDSup = false) = 0;

  //@}

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const = 0;

  /** Default function for computing an indicator of convergence
   *  \return a double when DS is a Lagrangian
   */
  virtual double dsConvergenceIndicator() ;

  /** set R to zero
   */
  virtual void resetAllNonSmoothPart() = 0;

  /** set R to zero for a given level
   * \param level
   */
  virtual void resetNonSmoothPart(unsigned int level) = 0;

  virtual void endStep() {};

  /** Get _pluging
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginG() const
  {
    return _pluging;
  };

  /** Get _pluginJacgx
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginJacGX() const
  {
    return _pluginJacgx;
  };

  /** Get _pluginJacxDotG
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginJacXDotG() const
  {
    return _pluginJacxDotG;
  };

  /** Initialize the workspace elements
   * \param workVector the vectors needed for the integration
   * \param workMatrices the matrices needed for the integration
   */
  virtual void initWorkSpace(VectorOfVectors& workVector, VectorOfMatrices& workMatrices) {};
//  virtual void initWorkSpace(VectorOfVectors& workVector, VectorOfMatrices& workMatrices) = 0;

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(DynamicalSystem);

};


#endif // DynamicalSystem_H

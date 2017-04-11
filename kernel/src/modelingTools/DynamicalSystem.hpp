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

#include <iostream>
/** Abstract interface to Dynamical Systems

    \author SICONOS Development Team - copyright INRIA
    \date (Creation) January 15, 2007


    This class is used to describe dynamical systems of the form :
    \f[
    g(\dot x, x, t, z) = 0
    \f]
    where

     - \f$ x \in R^{n} \f$ is the state.

     - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic
     variables, some sort of discret state.  For example, z may be used
     to set some perturbation parameters, to control the system (z
     set by actuators) and so on.
     - \f$ g : R^{n} \times R  \to  R^{n}   \f$ .

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

*/

class DynamicalSystem
{

public:
  /** List of indices used to save tmp work vectors
   * The last value is the size of the present list, so you HAVE to leave it at the end position.
   */
  enum DSWorkVectorId {local_buffer, freeresidu, free, acce_memory, acce_like, sizeWorkV};

private:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(DynamicalSystem);

  /** used to set ds number */
  static unsigned int __count;

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

  /** state of the system,
   *  \f$  x \in R^{n}\f$ - With _x[0]=\f$ x \f$ , _x[1]= \f$ \dot{x} \f$ . */
  VectorOfVectors _x;

  /** jacobian according to x of the right-hand side (\f$ rhs = \dot x =
      f(x,t) + r \f$) */
  SP::SiconosMatrix _jacxRhs;

  /** Arbitrary algebraic values vector, z, discrete state of the
      system. */
  SP::SiconosVector _z;
  
  /** the  previous state vectors stored in memory 
   */
  SP::SiconosMemory _xMemory;

  /** number of previous states stored in memory */
  unsigned int _stepsInMemory;

  // ===== CONSTRUCTORS =====

  /** default constructor */
  DynamicalSystem();

  /** minimal constructor, from state dimension
      result in \f$ \dot x = r \f$
   *  \param dimension size of the system (n)
   */
  DynamicalSystem(unsigned int dimension);

  /** Copy constructor
   * \param ds the DynamicalSystem to copy
   */
  DynamicalSystem(const DynamicalSystem & ds);

  /** Initialize all PluggedObject whether they are used or not.
   */
  virtual void _zeroPlugin() = 0;

  /** Common code for constructors
      should be replaced in C++11 by delegating constructors
      \param intial_state vector of initial values for state
   */
  void _init();

public:

  /** destructor */
  virtual ~DynamicalSystem() {};

  /*! @name Right-hand side computation */
  //@{

  /** allocate (if needed)  and compute rhs and its jacobian.
   * \param time of initialization
   */
  virtual void initRhs(double time) = 0 ;

  /** set nonsmooth input to zero
   *  \param int input-level to be initialized.
   */
  virtual void initializeNonSmoothInput(unsigned int level) = 0;

  /** compute all component of the dynamical system, for the current state.
   *  \param double time of interest
   */
  void update(double time);

  /** update right-hand side for the current state
   *  \param double time of interest
   *  \param bool isDSup flag to avoid recomputation of operators
   */
  virtual void computeRhs(double time, bool isDSup = false) = 0;

  /** update \f$\nabla_x rhs\f$ for the current state
   *  \param double time of interest
   *  \param bool isDSup flag to avoid recomputation of operators
   */
  virtual void computeJacobianRhsx(double time , bool isDSup = false) = 0;

  /** reset nonsmooth part of the rhs, for all 'levels' */
  virtual void resetAllNonSmoothParts() = 0;

  /** set nonsmooth part of the rhs to zero for a given level
   * \param level
   */
  virtual void resetNonSmoothPart(unsigned int level) = 0;

  ///@}

  /*! @name Attributes access

    For each 'Member' : \n
    - Member() returns a pointer to the object
    - getMember() a copy of the object
    - setMember() set the content of Member with a copy
    - setMemberPtr() set a pointer link to Member

    @{ */

  /** returns the id of the dynamical system */
  inline int number() const
  {
    return _number;
  }

  /** set the id of the DynamicalSystem
   *  \return the previous value of number
   */
  inline int setNumber(int new_number)
  {
    int old_n = _number;
    _number = new_number;
    return old_n;
  }

  /** returns the size of the vector state x */
  inline unsigned int n() const
  {
    return _n;
  }

  /** returns the dimension of the system (n for first order, ndof for Lagrangian).
   * Useful to avoid if(typeOfDS) when size is required.
   */
  virtual inline unsigned int dimension() const
  {
    return _n;
  };

  /** returns a pointer to the initial state vector */
  inline SP::SiconosVector x0() const
  {
    return _x0;
  };

  /** get a copy of the initial state vector */
  inline const SiconosVector getX0() const
  {
    return *_x0;
  }

  /** set initial state (copy)
   *  \param a SiconosVector
   */
  void setX0(const SiconosVector& newValue);

  /** set initial state (pointer link)
   *  \param newPtr SP::SiconosVector
   */
  void setX0Ptr(SP::SiconosVector newPtr);

  /** returns a pointer to the state vector \f$ x \f$
   *  \return SP::SiconosVector
   */
  inline SP::SiconosVector x() const
  {
    return _x[0];
  }
 
  /** get a copy of the current state vector \f$ x \f$
   * \return SiconosVector
   */
  inline const SiconosVector& getx() const
  {
    return *(_x[0]);
  }

  /** set content of current state vector \f$ x \f$
   *  \param newValue SiconosVector 
   */
  void setX(const SiconosVector& newValue);

  /** set state vector \f$ x \f$ (pointer link)
   *  \param newPtr SP::SiconosVector 
   */
  void setXPtr(SP::SiconosVector newPtr);

  /** returns a pointer to r vector (input due to nonsmooth behavior)
   *  \return SP::SiconosVector
   */
  inline SP::SiconosVector r() const
  {
    return _r;
  }

  /** get a copy of r vector (input due to nonsmooth behavior)
   *  \return a SiconosVector
   */
  inline const SiconosVector getR() const
  {
    return *_r;
  }

  /** set r vector (input due to nonsmooth behavior) content (copy)
   *  \param newValue SiconosVector 
   */
  void setR(const SiconosVector& newValue );

  /** set r vector (input due to nonsmooth behavior) (pointer link)
   *  \param newPtr SP::SiconosVector newPtr
   */
  void setRPtr(SP::SiconosVector newPtr);

  /** returns a pointer to the right-hand side vector, (i.e. \f$ \dot x \f$)
   *  \return SP::SiconosVector
   */
  inline SP::SiconosVector rhs() const
  {
    return _x[1];
  }

  /** get a copy of the right-hand side vector, (i.e. \f$ \dot x \f$)
   *  \return SiconosVector
   */
  inline SiconosVector& getRhs() const
  {
    return *(_x[1]);
  }

  /** set the value of the right-hand side, \f$ \dot x \f$
   *  \param newValue SiconosVector
   */
  virtual void setRhs(const SiconosVector& newValue);

  /** set right-hand side, \f$ \dot x \f$ (pointer link)
   *  \param newPtr SP::SiconosVector
   */
  virtual void setRhsPtr(SP::SiconosVector newPtr);

  /** returns a pointer to $\nabla_x rhs()$
   *  \return SP::SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianRhsx() const
  {
    return _jacxRhs;
  }

  /** set the value of \f$\nabla_x rhs()\f$$
   *  \param newValue SiconosMatrix
   */
  void setJacobianRhsx(const SiconosMatrix& newValue);

  /** set \f$\nabla_x rhs()\f$, pointer link
   *  \param newPtr SP::SiconosMatrix  
   */
  void setJacobianRhsxPtr(SP::SiconosMatrix newPtr);

  /** returns a pointer to \f$ z \f$, the vector of algebraic parameters.
   *  \return SP::SiconosVector
   */
  inline SP::SiconosVector z() const
  {
    return _z;
  }

  /** get a copy of \f$ z \f$, the vector of algebraic parameters.
   * \return a SiconosVector
   */
  inline const SiconosVector& getz() const
  {
    return *_z;
  }

  /** set the value of \f$ z \f$ (copy)
   *  \param newValue SiconosVector 
   */
  void setz(const SiconosVector& newValue) ;

  /** set \f$ z \f$ (pointer link)
   *  \param newPtr SP::SiconosVector 
   */
  void setzPtr(SP::SiconosVector newPtr);

  /** @} end of members access group. */

  /*! @name Memory vectors management  */
  //@{

  
  /** returns saved values of state vector, if any
   *  \return SP::SiconosMemory
   */
  inline SP::SiconosMemory xMemory() const
  {
    return _xMemory;
  }

  /** returns the number of step saved in memory for state vector
   *  \return int
   */
  inline int stepsInMemory() const
  {
    return _stepsInMemory;
  }
  
  /** set number of steps to be saved
   *  \param int steps
   */
  inline void setStepsInMemory(int steps)
  {
    _stepsInMemory = steps;
  }
  
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
  
  /** call all plugged functions for the current state
   * \param time  the current time
   */
  virtual void updatePlugins(double time) = 0;
  
  ///@}

  /*! @name Miscellaneous public methods */
  //@{

  /** reset the global DynamicSystem counter (for ids)
   *  \return the previous value of count
   */
  static inline int resetCount(int new_count=0)
  {
    int old_count = __count;
    __count = new_count;
    return old_count;
  };

  /** reset the state x() to the initial state x0 */
  virtual void resetToInitialState();
  
  /** True if the system is linear.
   * \return a boolean
   */
  virtual bool isLinear()
  {
    return false;
  };

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const = 0;

  ///@}

  //visitors hook
  VIRTUAL_ACCEPT_VISITORS(DynamicalSystem);
  
};


#endif // DynamicalSystem_H

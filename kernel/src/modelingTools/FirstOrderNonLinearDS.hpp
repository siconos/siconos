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

/*! \file FirstOrderNonLinearDS.hpp
  \brief First Order Non Linear Dynamical Systems
*/

#ifndef FIRSTORDERNONLINEARDS_H
#define FIRSTORDERNONLINEARDS_H

#include "DynamicalSystem.hpp"


/**  General First Order Non Linear Dynamical Systems - \f$ M(t) \dot{x} = f(x,t,z) + r, \quad x(t_0) = x_0 \f$
     
     \author SICONOS Development Team - copyright INRIA
     \date (Creation) April 29, 2004
     
     This class defines and computes a generic n-dimensional
     dynamical system of the form :
     
     \f[
     M \dot x = f(x,t,z) + r, \quad x(t_0) = x_0
     \f]

     where

     - \f$ x \in R^{n} \f$ is the state.
     - \f$ M \in R^{n\times n}\f$ a "mass matrix"
     - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
     - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic
     variables, some sort of discret state.  For example, z may be used
     to set some perturbation parameters, to control the system (z
     set by actuators) and so on.

     - \f$ f : R^{n} \times R  \mapsto  R^{n}\f$ the vector field.
 
  By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
  and the initial conditions are given by
   \f[
   x(t_0)=x_0
  \f]
  To define a Boundary Value Problem, a pointer on a BoundaryCondition must be set.
   
  The right-hand side and its jacobian (from base classe) are defined as

  \f[
  rhs = \dot x =  M^{-1}(f(x,t,z)+ r) \\
  jacobianRhsx = \nabla_x rhs(x,t,z) = M^{-1}\nabla_x f(x,t,z)
  \f]


  The following operators can be plugged, in the usual way (see User Guide)
  
  - \f$f(x,t,z)\f$
  - \f$\nabla_x f(x,t,z)\f$
  - \f$M(t)\f$
  
 
 */
class FirstOrderNonLinearDS : public DynamicalSystem
{

private:

  /** plugin signature */
  typedef void (*FNLDSPtrfct)(double, unsigned int, const double*, double*, unsigned int, double*);
 
protected:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(FirstOrderNonLinearDS);

  /** Common code for constructors
      should be replaced in C++11 by delegating constructors
      \param intial_state vector of initial values for state
  */
  void _init(SP::SiconosVector initial_state);

  /** Matrix coefficient of \f$ \dot x \f$ */
  SP::SiconosMatrix _M;

  /** value of f(x,t,z) */
  SP::SiconosVector _f;

  // Note FP: isn't it strange to define b in this class here rather than in Linear derived class?
  /** strength vector */
  SP::SiconosVector _b;

  /** to store f(x_k,t_k,z_k)*/
  SP::SiconosVector _fold;

  /** Gradient of \f$ f(x,t,z) \f$ with respect to \f$ x\f$*/
  SP::SiconosMatrix _jacobianfx;

  /** DynamicalSystem plug-in to compute f(x,t,z)
   *  \param current time
   *  \param size of the vector _x
   *  \param[in,out] pointer to the first element of the vector _x
   *  \param[in,out] the pointer to the first element of the vector _f
   *  \param the size of the vector _z
   *  \param a vector of parameters _z
   */
  SP::PluggedObject _pluginf;

  /** DynamicalSystem plug-in to compute the gradient of f(x,t,z) with respect to the state: \f$ \nabla_x f: (x,t,z) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   * \param time current time
   * \param sizeOfX size of vector x
   * \param x pointer to the first element of x
   * \param[in,out] jacob pointer to the first element of jacobianfx matrix
   * \param  the size of the vector z
   * \param[in,out]  a vector of parameters, z
   */
  SP::PluggedObject _pluginJacxf;

  SP::PluggedObject _pluginM;

  /**  the previous r vectors */
  SP::SiconosMemory _rMemory;

  /** Copy of M Matrix, LU-factorized, used to solve systems like Mx = b with LU-factorization.
      (Warning: may not exist, used if we need to avoid factorization in place of M) */
  SP::SiconosMatrix _invM;

  /** default constructor */
  FirstOrderNonLinearDS(): DynamicalSystem() {};

  /** Reset the PluggedObjects */
  virtual void _zeroPlugin();


public:

  // ===== CONSTRUCTORS =====

  /** constructor from initial state, leads to \f$ \dot x = r\f$
      \param newX0 initial state
      \warning you need to set explicitely the plugin for f and its jacobian if needed (e.g. if used with
      an EventDriven scheme)
  */
  FirstOrderNonLinearDS(SP::SiconosVector newX0);

  /** constructor from initial state and f (plugins), \f$\dot x = f(x, t, z) + r\f$
   *  \param newX0 initial state
   *  \param fPlugin name of the plugin function to be used for f(x,t,z)
   *  \param jacobianfxPlugin name of the plugin to be used for the jacobian of f(x,t,z)
   */
  FirstOrderNonLinearDS(SP::SiconosVector newX0, const std::string& fPlugin, const std::string& jacobianfxPlugin);

  /** Copy constructor
   * \param FONLDS the FirstOrderNonLinearDS to copy
   */
  FirstOrderNonLinearDS(const FirstOrderNonLinearDS & FONLDS);

  /** destructor */
  virtual ~FirstOrderNonLinearDS() {};

  /*! @name Right-hand side computation */
  //@{

  /** allocate (if needed)  and compute rhs and its jacobian.
   * \param time of initialization
   */
  void initRhs(double time);

  /** set nonsmooth input to zero
   *  \param int input-level to be initialized.
   */
  void initializeNonSmoothInput(unsigned int level) ;

  /** update right-hand side for the current state
   *  \param double time of interest
   *  \param bool isDSup flag to avoid recomputation of operators
   */
  void computeRhs(double time, bool isDSUp = false);

  /** update \f$\nabla_x rhs\f$ for the current state
   *  \param double time of interest
   *  \param bool isDSup flag to avoid recomputation of operators
   */
  void computeJacobianRhsx(double time, bool isDSUp = false);

  /** reset non-smooth part of the rhs (i.e. r), for all 'levels' */
  virtual void resetAllNonSmoothParts();

  /** set nonsmooth part of the rhs (i.e. r) to zero for a given level
   * \param level
   */
  virtual void resetNonSmoothPart(unsigned int level);

  ///@}

  /*! @name Attributes access 
    @{ */

  /** returns a pointer to M, matrix coeff. on left-hand side
  */
  inline SP::SiconosMatrix M() const
  {
    return _M;
  }

  /** set M, matrix coeff of left-hand side (pointer link)
   *  \param newM the new M matrix
   */
  inline void setMPtr(SP::SiconosMatrix newM)
  {
    _M = newM;
  }

  // --- invM ---
  /** get a copy of the LU factorisation of M operator
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getInvM() const
  {
    return *_invM;
  }

  /** get the inverse of LU fact. of M operator (pointer link)
   *  \return pointer to a SiconosMatrix
   */
  inline SP::SiconosMatrix invM() const
  {
    return _invM;
  }

  /** returns f(x,t,z) (pointer link)
   */
  inline SP::SiconosVector f() const
  {
    return _f;
  }

  /** set f(x,t,z) (pointer link)
   *  \param newPtr a SP::SiconosVector
   */
  inline void setFPtr(SP::SiconosVector newPtr)
  {
    _f = newPtr;
  }

  /** get jacobian of f(x,t,z) with respect to x (pointer link)
   *  \return SP::SiconosMatrix
   */
  virtual SP::SiconosMatrix jacobianfx() const
  {
    return _jacobianfx;
  }

  /** set jacobian of f(x,t,z) with respect to x (pointer link)
   *  \param newPtr the new value
   */
  inline void setJacobianfxPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianfx = newPtr;
  }

  /** get b vector (pointer link)
   *  \return a SP::SiconosVector
   */
  inline SP::SiconosVector b() const
  {
    return _b;
  }

  /** set b vector (pointer link)
   *  \param b a SiconosVector
   */
  inline void setbPtr(SP::SiconosVector b)
  {
    _b = b;
  }

  /** set b vector (copy)
   *  \param b a SiconosVector
   */
  void setb(const SiconosVector& b);

  /** @} end of members access group. */

  /*! @name Memory vectors management  */
  //@{

  /** get all the values of the state vector r stored in memory
   *  \return a memory vector
   */
  inline SP::SiconosMemory rMemory() const
  {
    return _rMemory;
  }

  /** returns previous value of rhs -->OSI Related!!*/
  inline SP::SiconosVector fold() const
  {
    return _fold;
  }

  /** initialize the SiconosMemory objects: reserve memory for i
      vectors in memory and reset all to zero.
      *  \param steps the size of the SiconosMemory (i)
      */
  void initMemory(unsigned int steps);

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  void swapInMemory();

  //@}

  /*! @name Plugins management  */
  
  //@{
  
  /** Call all plugged-function to initialize plugged-object values
      \param time value
  */
  virtual void updatePlugins(double time);

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute M
   *  \param pluginPath the complete path to the plugin
   *  \param functionName function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeMFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute M
   *  \param fct a pointer on the plugin function
   */
  void setComputeMFunction(FPtr1 fct);

  /** to set a specified function to compute f(x,t)
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute the vector f
   *  \param fct a pointer on the plugin function
   */
  void setComputeFFunction(FPtr1 fct);

  /** to set a specified function to compute jacobianfx
   *  \param pluginPath the complete path to the plugin
   *  \param functionName function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianfxFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobianfx
   *  \param fct a pointer on the plugin function
   */
  void setComputeJacobianfxFunction(FPtr1 fct);

  // --- compute plugin functions ---

  /** Default function to compute \f$ M: (x,t)\f$
   * \param time time instant used in the computations
   */
  void computeM(double time);

  /** Default function to compute \f$ f: (x,t)\f$
   * \param time time instant used in the computations
   */
  //virtual void computef(double time);

  /** function to compute \f$ f: (x,t)\f$
   * \param time time instant used in the computations
   * \param state x value
   */
  virtual void computef(double time, SP::SiconosVector state);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n}
   *   \times R \mapsto R^{n \times n} \f$ with x different from
   *   current saved state.
   *  \param time instant used in the computations
   *  \param x2 a SiconosVector to store the resuting value
   */
  virtual void computeJacobianfx(double time, SP::SiconosVector state);

  /** Get _pluginf
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginF() const
  {
    return _pluginf;
  };

  /** Get _pluginJacxf
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginJacxf() const
  {
    return _pluginJacxf;
  };

  /** Get _pluginM
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginM() const
  {
    return _pluginM;
  };

  ///@}

  /*! @name Miscellaneous public methods */
  //@{

  /** print the data of the dynamical system on the standard output
   */
  void display() const;
  ///@}


  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderNonLinearDS)

#endif



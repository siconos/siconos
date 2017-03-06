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

/*! \file FirstOrderLinearDS.hpp
 */
#ifndef FirstOrderLinearDS_H
#define FirstOrderLinearDS_H

#include "FirstOrderNonLinearDS.hpp"


/** First Order Linear Systems - \f$M(t) \dot x = A(t)x(t)+ b(t) + r, \quad x(t_0)=x_0\f$.
 
    \author SICONOS Development Team - copyright INRIA
    \date (Creation) Apr 29, 2004
 
    This class represents first order linear systems of the form:

    \f[
    M(t) \dot x = A(t)x(t)+ b(t) + r,
    x(t_0)=x_0
    \f]

    where

    - \f$x \in R^{n} \f$ is the state,
    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
    - \f$M \in R^{n\times n} \f$ is an invertible matrix
    - \f$A \in R^{n\times n}\f$
    - \f$b \in R^{n} \f$
 

    The following operators can be plugged, in the usual way (see User Guide)
  
    - \f$A(t)\f$
    - \f$b(t)\f$
    - \f$M(t)\f$

*/
class FirstOrderLinearDS : public FirstOrderNonLinearDS
{
protected:

  /* serialization hooks */
  ACCEPT_SERIALIZATION(FirstOrderLinearDS);


  /** matrix specific to the FirstOrderLinearDS \f$ A \in R^{n \times n}  \f$*/
  SP::SiconosMatrix _A;
  
  /** FirstOrderLinearDS plug-in to compute A(t,z), id = "A"
   * @param time : current time
   * @param sizeOfA : size of square-matrix A
   * @param[in,out] A : pointer to the first element of A
   * @param size of vector z
   * @param[in,out] z a vector of user-defined parameters
   */
  SP::PluggedObject _pluginA;

  /** FirstOrderLinearDS plug-in to compute b(t,z), id = "b"
   * @param time : current time
   * @param sizeOfB : size of vector b
   * @param[in,out] b : pointer to the first element of b
   * @param size of vector z
   * @param[in,out] param  : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginb;

  /** default constructor
   */
  FirstOrderLinearDS(): FirstOrderNonLinearDS() {};

  /** Reset all the plugins */
  virtual void _zeroPlugin();

public:

  /** plugin signature */
  typedef   void (*LDSPtrFunction)(double, unsigned int, double*, unsigned int, double*);
  
  /** constructor from initial state and plugins
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param APlugin plugin for A
   *  \param bPlugin plugin for b
   */
  FirstOrderLinearDS(SP::SiconosVector newX0, const std::string& APlugin, const std::string& bPlugin);

  /** constructor from initial state and plugin for A
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param newA matrix A
   */
  FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA);

  /** constructor from initial state
   *  \param newX0 the initial state of this DynamicalSystem
   */
  FirstOrderLinearDS(SP::SiconosVector newX0);

  /** constructor from a initial state and constant matrices
   *  \param newX0 the initial state of this DynamicalSystem
   *  \param newA matrix A
   *  \param newB b
   */
  FirstOrderLinearDS(SP::SiconosVector newX0, SP::SiconosMatrix newA, SP::SiconosVector newB);

  /** Copy constructor
   * \param FOLDS the original FirstOrderLinearDS we want to copy
   */
  FirstOrderLinearDS(const FirstOrderLinearDS & FOLDS);

  /** destructor */
  virtual ~FirstOrderLinearDS() {};

  /*! @name Right-hand side computation */
  
  /** Initialization function for the rhs and its jacobian.
   *  \param time time of initialization.
   */
  void initRhs(double time) ;

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

  ///@}

  /*! @name Attributes access 
    @{ */
  
  /** get the matrix \f$A\f$
   *  \return pointer (SP) on a matrix
   */
  inline SP::SiconosMatrix A() const
  {
    return _A;
  }

  /** get jacobian of f(x,t,z) with respect to x (pointer link)
   *  \return SP::SiconosMatrix
   */
  virtual SP::SiconosMatrix jacobianfx() const
  {
    return _A;
  };

  /** set A to pointer newPtr
   *  \param newA the new A matrix
   */
  inline void setAPtr(SP::SiconosMatrix newA)
  {
    _A = newA;
  }

  /** set A to a new matrix
   * \param newA the new A matrix
   **/
  void setA(const SiconosMatrix& newA);

  // --- plugins related functions
  /*! @name Plugins management  */
  
  //@{

  /** Call all plugged-function to initialize plugged-object values
      \param time value
  */
  virtual void updatePlugins(double time);

  /** set a specified function to compute the matrix A => same action as setComputeJacobianfxFunction
   *  \param pluginPath the complete path to the plugin
   *  \param functionName the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute the matrix A
   *  \param fct a pointer on a function
   */
  void setComputeAFunction(LDSPtrFunction fct);

  /** set a specified function to compute the vector b
   *  \param pluginPath the complete path to the plugin file
   *  \param functionName the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputebFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute the vector b
   *  \param fct a pointer on a function
   */
  void setComputebFunction(LDSPtrFunction fct);

  /** default function to compute matrix A => same action as
      computeJacobianfx
      \param time time instant used to compute A
  */
  void computeA(double time);

  /** default function to compute vector b
   * \param time time instant used to compute b
   */
  virtual void computeb(double time);
  
  /** Get _pluginA
   * \return the plugin for A
   */
  inline SP::PluggedObject getPluginA() const
  {
    return _pluginA;
  };

  /** Get _pluginb
   * \return the plugin for b
   */
  inline SP::PluggedObject getPluginb() const
  {
    return _pluginb;
  };

  /** Set _pluginA
   * \param newPluginA the new plugin
   */
  inline void setPluginA(SP::PluggedObject newPluginA)
  {
    _pluginA = newPluginA;
  };

  /** Set _pluginb
   * \param newPluginB the new plugin
   */
  inline void setPluginB(SP::PluggedObject newPluginB)
  {
    _pluginb = newPluginB;
  };
  
  ///@}

  /*! @name Miscellaneous public methods */

  /** data display on screen
   */
  void display() const;

  /** True if the system is linear.
   * \return a boolean
   */
  virtual bool isLinear()
  {
    return true;
  }

  ///@}
  
  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderLinearDS)

#endif // FirstOrderLinearDS_H

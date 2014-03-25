/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file FirstOrderNonLinearDS.hpp
  First Order Non Linear Dynamical Systems
*/

#ifndef FIRSTORDERNONLINEARDS_H
#define FIRSTORDERNONLINEARDS_H

#include "DynamicalSystem.hpp"


typedef void (*FNLDSPtrfct)(double, unsigned int, const double*, double*, unsigned int, double*);

/**  General First Order Non Linear Dynamical Systems - \f$ M \dot{x} = f(x,t,z) + r, \quad x(t_0) = x_0\f$
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) April 29, 2004
 *
 * This class defines and computes a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * M \dot x = f(x,t,z) + r, \quad
 * x(t_0) = x_0
 * \f]
 * where
 *    - \f[ x \in R^{n} \f] is the state.
 *    - \f$ M \in R^{n\times n} a "mass matrix"
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discrete state.
 *      For example, z may be used to set some perturbation parameters, or anything else.
 *
 *  with \f$ f : R^{n} \times R  \mapsto  R^{n}   \f$ the vector field.
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  \f[
 *  x(t_0)=x_0
 * \f]
 * To define a Boundary Value Problem, the pointer on  a BoundaryCondition must be set.
 *
 * \f$ f(x,t) \f$ is a plug-in function, and can be computed using computef(t).
 * Its Jacobian according to x is denoted jacobianfx, and computed thanks to computeJacobianfx(t).
 * f and jacobianfx can be plugged to external functions thanks to setComputeFFunction/setComputeJacobianfxFunction.
 *
 * Right-hand side of the equation is computed thanks to computeRhs(t).
 *
 * \f[
 *    \dot x =  M^{-1}(f(x,t,z)+ r)
 * \f]
 *
 * Its Jacobian according to x is jacobianRhsx:
 *
 *  \f[
 *   jacobianRhsx = \nabla_x rhs(x,t,z) = M^{-1}\nabla_x f(x,t,z)
 *  \f]
 *
 * At the time:
 *  - M is considered to be constant. (ie no plug-in, no jacobian ...)
 *  - M is not allocated by default. The only way to use M is setM or setMPtr.
 *
 */
class FirstOrderNonLinearDS : public DynamicalSystem
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(FirstOrderNonLinearDS);


  /** Matrix coefficient of \f$ \dot x \f$ */
  SP::SiconosMatrix _M;

  /** value of f(x,t,z) */
  SP::SiconosVector _f;

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

  /** Residu r*/
  SP::SiconosVector _residur;

  /** g_alpha*/
  SP::SiconosVector _g_alpha;
  //  SP::SiconosVector mXfree;

  SP::SiconosVector _xp;
  SP::SiconosVector _xq;




  /** Copy of M Matrix, used to solve systems like Mx = b with LU-factorization.
      (Warning: may not exist, used if we need to avoid factorization in place of M) */
  SP::SiconosMatrix _invM;

  /** default constructor
   */
  FirstOrderNonLinearDS(): DynamicalSystem() {};


public:

  // ===== CONSTRUCTORS =====

  /** constructor from a set of data
      \param newX0 initial state of this DynamicalSystem
      \warning you need to set yoursel the plugin for f and also for the
      jacobian if you use a EventDriven scheme
  */
  FirstOrderNonLinearDS(SP::SiconosVector newX0);

  /** constructor from a set of data
   *  \param newX0 initial state of this DynamicalSystem
   *  \param fPlugin plugin name for f of this DynamicalSystem
   *  \param jacobianfxPlugin plugin name for jacobianfx of this DynamicalSystem
   */
  FirstOrderNonLinearDS(SP::SiconosVector newX0, const std::string& fPlugin, const std::string& jacobianfxPlugin);

  /** Copy consctructor
   * \param FONLDS the FirstOrderNonLinearDS to copy
   */
  FirstOrderNonLinearDS(const FirstOrderNonLinearDS & FONLDS);

  /** destructor
   */
  virtual ~FirstOrderNonLinearDS() {};

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  inline SP::SiconosVector gAlpha() const
  {
    return _g_alpha;
  }


  /* Get _residur
   * \return _residur
   */
  inline SP::SiconosVector residur() const
  {
    return _residur;
  }


  // rMemory

  /** get all the values of the state vector r stored in memory
   *  \return a memory vector
   */
  inline SP::SiconosMemory rMemory() const
  {
    return _rMemory;
  }

  /** get M
   *  \return pointer on a plugged-matrix
   */
  inline SP::SiconosMatrix M() const
  {
    return _M;
  }

  /** set M to a new value
   *  \param newM the new M matrix
   */
  inline void setMPtr(SP::SiconosMatrix newM)
  {
    _M = newM;
  }

  // --- invM ---
  /** get the value of invM
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getInvMSimple() const
  {
    return *_invM;
  }

  /** get invM
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix invM() const
  {
    return _invM;
  }

  /** set the value of invM to newValue
   *  \param newValue the new value
   */
  void setInvM(const SiconosMatrix& newValue);

  /** link invM with a new pointer
   *  \param newPtr the new value
   */
  void setInvMPtr(SP::SiconosMatrix newPtr);

  // --- f ---

  /** get f
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector f() const
  {
    return _f;
  }
  inline SP::SiconosVector fold() const
  {
    return _fold;
  }

  /** set f to pointer newPtr
   *  \param newPtr a SP::SiconosVector
   */
  inline void setFPtr(SP::SiconosVector newPtr)
  {
    _f = newPtr;
  }

  /** get jacobianfx
   *  \return SP::SiconosMatrix
   */
  virtual SP::SiconosMatrix jacobianfx() const
  {
    return _jacobianfx;
  }

  /** set jacobianfx to pointer newPtr
   *  \param newPtr the new value
   */
  inline void setJacobianfxPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianfx = newPtr;
  }

  /** Initialization function for the rhs and its jacobian.
   *  \param time the time of initialization
   */
  void initRhs(double time);

  /** Call all plugged-function to initialize plugged-object values
      \param time the time used in the computations
   */
  virtual void updatePlugins(double time);

  /** dynamical system initialization function except for _r :
   *  mainly set memory and compute value for initial state values.
   *  \param time time of initialisation, default value = 0
   *  \param sizeOfMemory the size of the memory, default size = 1.
   */
  void initialize(double time = 0, unsigned int sizeOfMemory = 1);

  /** dynamical system initialization function for NonSmoothInput _r
   *  \param level for _r
   */
  void initializeNonSmoothInput(unsigned int level) ;


  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** initialize the SiconosMemory objects: reserve memory for i
      vectors in memory and reset all to zero.
   *  \param steps the size of the SiconosMemory (i)
   */
  void initMemory(unsigned int steps);

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  void swapInMemory();

  // ===== COMPUTE PLUGINS FUNCTIONS =====

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute M
   *  \param pluginPath the complete path to the plugin
   *  \param functionName function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeMFunction(const std::string& pluginPath, const std::string& functionName);

  /** set a specified function to compute M
   *  \param FPtr1 a pointer on the plugin function
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

  /** function to compute \f$ M: (x,t)\f$ with x different from current saved state.
   * \param time time instant used in the computations
   * \param x2 a SiconosVector to store the resuting value
   */
  void computeM(double time, SP::SiconosVector x2);

  /** Default function to compute \f$ f: (x,t)\f$
   * \param time time instant used in the computations
   */
  virtual void computef(double time);

  /** function to compute \f$ f: (x,t)\f$ with x different from current saved state.
   * \param time time instant used in the computations
   * \param x2 a SiconosVector to store the resuting value
   */
  virtual void computef(double time, SP::SiconosVector x2);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param time time instant used in the computations
   *  \param isDSup flag to avoid recomputation of operators
   *
   */
  virtual void computeJacobianfx(double time, bool isDSup = false);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n}
   *   \times R \mapsto R^{n \times n} \f$ with x different from
   *   current saved state.
   *  \param time instant used in the computations
   *  \param x2 a SiconosVector to store the resuting value
   */
  virtual void computeJacobianfx(double time, SP::SiconosVector x2);

  /** Default function to the right-hand side term
   *  \param time time instant used in the computations
   *  \param isDSup flag to avoid recomputation of operators
   *
   */
  void computeRhs(double time, bool isDSUp = false);

  /** Default function to jacobian of the right-hand side term according to x.
   *  Required when using an EventDriven Simulation.
   *  \param time instant used in the computations
   *  \param isDSup flag to avoid recomputation of operators
   *
   */
  void computeJacobianRhsx(double time, bool isDSUp = false);

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  void display() const;

  /** set R to zero
   */
  virtual void resetAllNonSmoothPart();

  /** set R to zero fo a given level
   * \param level the level to reset
   */
  virtual void resetNonSmoothPart(unsigned int level);

  /*
   * get the Xp work vector.
   */
  inline SP::SiconosVector xp() const
  {
    return _xp;
  };
  /*
   * get the Xq work vector.
   */
  inline SP::SiconosVector xq() const
  {
    return _xq;
  };

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

  /** Reset the PluggedObjects */
  virtual void zeroPlugin();

  /** To compute \f$\frac{|x_{i+1} - xi|}{|x_i|}\f$ where
      \f$x_{i+1}\f$ represents the present state and \f$x_i\f$ the
      previous one
   * \return a double
   */
  /*  double dsConvergenceIndicator(); */

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderNonLinearDS)

#endif



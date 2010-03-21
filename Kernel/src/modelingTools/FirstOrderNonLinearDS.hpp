/* Siconos-Kernel, Copyright INRIA 2005-2010.
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

/*! \file FirstOrderNonLinearDS.h
  First Order Non Linear Dynamical Systems
*/

#ifndef FIRSTORDERNONLINEARDS_H
#define FIRSTORDERNONLINEARDS_H

#include "DynamicalSystem.hpp"
#include "BlockMatrix.hpp"

class DynamicalSystem;
class BlockMatrix;

typedef void (*FNLDSPtrfct)(double, unsigned int, const double*, double*, unsigned int, double*);

/**  General First Order Non Linear Dynamical Systems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) April 29, 2004
 *
 * This class defines and computes a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * M \dot x = f(x,t,z) + r,
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state.
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *  For example, z may be used to set some perturbation parameters, or to control the system (z will be set by some actuators) or anything else.
 *
 *  with \f$ f : R^{n} \times R  \mapsto  R^{n}   \f$ .
 *  and M a nXn matrix.
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  \f[
 *  x(t_0)=x_0
 * \f]
 * To define a boundary Value Problem, the pointer on  a BoundaryCondition must be set.
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
 *   jacobianRhsx = \nabla_xrhs(x,t,z) = M^{-1}\nabla_xf(x,t,z)
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

  /** Matrix coefficient of \f$ \dot x \f$ */
  SP::SiconosMatrix _M;

  /** f(x,t,z) */
  SP::SiconosVector _f;

  /** to store f(x_k,t_k,z_k)*/
  SP::SiconosVector _fold;

  /** Gradient of \f$ f(x,t,z) \f$ with respect to \f$ x\f$*/
  SP::SiconosMatrix _jacobianfx;

  /** DynamicalSystem plug-in to compute f(x,t,z) - id="f".
    *  @param  : current time
    *  @param  : the size of the vector x
    *  @param  : the pointer to the first element of the vector x
    *  @param  : the pointer to the first element of the vector f(x,t)
    *  @param  : the size of the vector z
    *  @param  : a vector of parameters, z
    */
  SP::PluggedObject _pluginf;
  //  FNLDSPtrfct _computeFPtr;
  //  std::string _pluginNameComputeFPtr;


  /** DynamicalSystem plug-in to compute the gradient of f(x,t,z) with respect to the state: \f$ \nabla_x f: (x,t,z) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   * @param time : current time
   * @param sizeOfX : size of vector x
   * @param x : pointer to the first element of x
   * @param[in,out] jacob : pointer to the first element of jacobianfx matrix
   * @param  : the size of the vector z
   * @param[in,out] z: a vector of parameters, z
   */
  SP::PluggedObject _pluginJacxf;
  //   FNLDSPtrfct computeJacobianfxPtr;
  //   std::string pluginNameComputeJacobianfxPtr;

  SP::PluggedObject _pluginM;

  //  FNLDSPtrfct  pluginComputeM;
  //  std::string pluginNamePluginComputeM;

  /**  the previous r vectors */
  SP::SiconosMemory _rMemory;

  /** Residu r*/
  SP::SimpleVector _residur;

  /** g_alpha*/
  SP::SimpleVector _g_alpha;
  //  SP::SiconosVector mXfree;

  SP::SiconosVector _xp;
  SP::SiconosVector _xq;




  /** Copy of M Matrix, used to solve systems like Mx = b with LU-factorization.
      (Warning: may not exist, used if we need to avoid factorization in place of M) */
  SP::SiconosMatrix _invM;

  /** default constructor
   */
  FirstOrderNonLinearDS(): DynamicalSystem() {};

  /** constructor from a set of data
      \param SiconosVector : initial state of this DynamicalSystem
  */
  FirstOrderNonLinearDS(SP::SiconosVector);
  FirstOrderNonLinearDS(const SiconosVector& SiconosVector);
  virtual void zeroPlugin();

public:

  // ===== CONSTRUCTORS =====

  /** xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   */
  FirstOrderNonLinearDS(SP::DynamicalSystemXML dsXML);

  /** constructor from a set of data
   *  \param SiconosVector : initial state of this DynamicalSystem
   *  \param string : plugin name for f of this DynamicalSystem
   *  \param string : plugin name for jacobianfx of this DynamicalSystem
   *  \exception RuntimeException
   */
  FirstOrderNonLinearDS(const SiconosVector&, const std::string&, const std::string&);

  // ===== DESTRUCTOR =====

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



  inline SP::SiconosVector residur() const
  {
    return _residur;
  }


  // rMemory

  /** get all the values of the state vector r stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory rMemory() const
  {
    return _rMemory;
  }

  // --- M ---
  /** get the value of M
   *  \return a plugged-matrix

  inline const PMJF getM() const { return *M; }
  */
  /** get M
   *  \return pointer on a plugged-matrix
   */
  inline SP::SiconosMatrix M() const
  {
    return _M;
  }

  /** set the value of M to newValue
   *  \param plugged-matrix newValue

  void setM(const PMJF&);
  */
  /** set M to pointer newPtr
   *  \param a plugged matrix SP
   */
  inline void setMPtr(SP::SiconosMatrix newPtr)
  {
    _M = newPtr;
  }

  // --- invM ---
  /** get the value of invM
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getInvMSimple() const
  {
    return *_invM;
  }

  /** get the value of invM
   *  \return BlockMatrix
   */
  inline const BlockMatrix getInvMBlock() const
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
   *  \param SimpleMatrix newValue
   */
  void setInvM(const SiconosMatrix&);

  /** link invM with a new pointer
   *  \param a pointer to SiconosMatrix
   */
  void setInvMPtr(SP::SiconosMatrix);

  // --- f ---

  /** get the value of f
   *  \return plugged vector

  inline const PVF getF() const { return *_f; }
  */

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

  /** set the value of f to newValue
   *  \param a plugged vector

  void setF(const PVF&);
  */
  /** set f to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setFPtr(SP::SiconosVector newPtr)
  {
    _f = newPtr;
  }

  // --- jacobianfx ---
  /** get the value of jacobianfx
   *  \return a plugged-matrix

  inline const PMJF getJacobianfx() const { return *jacobianfx; }
  */
  /** get jacobianfx
   *  \return pointer on a plugged-matrix
   */
  virtual SP::SiconosMatrix jacobianfx() const
  {
    return _jacobianfx;
  }

  /** set the value of jacobianfx to newValue
   *  \param plugged-matrix newValue

  void setJacobianfx(const PMJF&);
  */
  /** set jacobianfx to pointer newPtr
   *  \param a plugged matrix SP
   */
  inline void setJacobianfxPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianfx = newPtr;
  }

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization
   */
  void initRhs(double) ;

  /** Call all plugged-function to initialize plugged-object values
      \param time
   */
  virtual void updatePlugins(double);

  /** dynamical system initialization function: mainly set memory and
      compute value for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  void initialize(const std::string&, double = 0, unsigned int = 1) ;

  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** initialize the SiconosMemory objects: reserve memory for i
      vectors in memory and reset all to zero.
   *  \param the size of the SiconosMemory (i)
   */
  void initMemory(unsigned int) ;

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  void swapInMemory();

  // ===== COMPUTE PLUGINS FUNCTIONS =====

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute M
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeMFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute M
   *  \param FPtr1 : a pointer on the plugin function
   */
  void setComputeMFunction(FPtr1 fct);

  /** to set a specified function to compute f(x,t)
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFFunction(const std::string&  pluginPath, const std::string& functionName);

  /** set a specified function to compute the vector f
   *  \param FPtr1 : a pointer on the plugin function
   */
  void setComputeFFunction(FPtr1 fct);

  /** to set a specified function to compute jacobianfx
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianfxFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobianfx
   *  \param FPtr1 : a pointer on the plugin function
   */
  void setComputeJacobianfxFunction(FPtr1 fct);

  // --- compute plugin functions ---

  /** Default function to compute \f$ M: (x,t)\f$
   * \param double time : current time
   */
  void computeM(double);

  /** function to compute \f$ M: (x,t)\f$ with x different from current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  void computeM(double, SP::SiconosVector);

  /** Default function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   */
  virtual void computef(double);

  /** function to compute \f$ f: (x,t)\f$ with x different from current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  virtual void computef(double, SP::SiconosVector);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  void computeJacobianfx(double, bool  = false);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n}
   *   \times R \mapsto R^{n \times n} \f$ with x different from
   *   current saved state.
   *  \param double time : current time
   *  \param SP::SiconosVector
   */
  void computeJacobianfx(double, SP::SiconosVector);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  void computeRhs(double, bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  void computeJacobianRhsx(double, bool  = false);

  // ===== XML MANAGEMENT FUNCTIONS =====

  /** copy the data common to each system in the XML tree
   */
  void saveSpecificDataToXML();

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  void display() const;

  /** set R to zero
   */
  virtual void resetNonSmoothPart();
  /*
   * reset work vectors.
   */
  virtual void preparStep();

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
  /*
   * get the Xfree work vector.
   */
  //  inline SP::SiconosVector xfree() const { return mXfree;};

  /** To compute \f$\frac{|x_{i+1} - xi|}{|x_i|}\f$ where
      \f$x_{i+1}\f$ represents the present state and \f$x_i\f$ the
      previous one
   * \return a double
   */
  /*  double dsConvergenceIndicator(); */

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param SP::DynamicalSystem : the system which must be converted
   * \return a pointer on the system if it is of the right type, NULL otherwise
   */
  static FirstOrderNonLinearDS* convert(DynamicalSystem* ds);


  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(FirstOrderNonLinearDS);

#endif



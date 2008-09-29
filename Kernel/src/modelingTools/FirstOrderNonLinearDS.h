/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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

/*! \file FirstOrderNonLinearDS.h
First Order Non Linear Dynamical Systems
*/

#ifndef FIRSTORDERNONLINEARDS_H
#define FIRSTORDERNONLINEARDS_H

#include "DynamicalSystem.h"
#include "BlockMatrix.h"

class DynamicalSystem;
class BlockMatrix;

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
 * \f$ f(x,t) \f$ is a plug-in function, and can be computed using computeF(t).
 * Its Jacobian according to x is denoted jacobianXF, and computed thanks to computeJacobianXF(t).
 * f and jacobianXF can be plugged to external functions thanks to setComputeFFunction/setComputeJacobianXFFunction.
 *
 * Right-hand side of the equation is computed thanks to computeRhs(t).
 *
 * \f[
 *    \dot x =  M^{-1}(f(x,t,z)+ r)
 * \f]
 *
 * Its Jacobian according to x is jacobianXRhs:
 *
 *  \f[
 *   jacobianXRhs = \nabla_xrhs(x,t,z) = M^{-1}\nabla_xf(x,t,z)
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
  SiconosMatrixSPtr M;

  /** f(x,t,z) */
  SiconosVectorSPtr f;

  /** Gradient of \f$ f(x,t,z) \f$ with respect to \f$ x\f$*/
  SiconosMatrixSPtr jacobianXF;

  /** DynamicalSystem plug-in to compute f(x,t,z) - id="f".
   *  @param  : current time
   *  @param  : the size of the vector x
   *  @param  : the pointer to the first element of the vector x
   *  @param  : the pointer to the first element of the vector f(x,t)
   *  @param  : the size of the vector z
   *  @param  : a vector of parameters, z
   */
  void (*computeFPtr)(double, unsigned int, const double*, double*, unsigned int, double*);

  /** DynamicalSystem plug-in to compute the gradient of f(x,t,z) with respect to the state: \f$ \nabla_x f: (x,t,z) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   * @param time : current time
   * @param sizeOfX : size of vector x
   * @param x : pointer to the first element of x
   * @param[in,out] jacob : pointer to the first element of jacobianXF matrix
   * @param  : the size of the vector z
   * @param[in,out] z: a vector of parameters, z
   */
  void (*computeJacobianXFPtr)(double, unsigned int, const double*, double*, unsigned int, double*);

  /** the  input vector due to the non-smooth law \f$  r \in R^{n}\f$ (multiplier, force, ...)*/
  SiconosVectorSPtr r;

  /**  the previous r vectors */
  SiconosMemorySPtr rMemory;

  /** Copy of M Matrix, used to solve systems like Mx = b with LU-factorization.
      (Warning: may not exist, used if we need to avoid factorization in place of M) */
  SiconosMatrixSPtr invM;


#ifndef WithSmartPtr
  /** set all allocation flags (isAllocated map)
   *  \param bool
   */
  virtual void initAllocationFlags(bool);
#endif

  /** set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  virtual void initPluginFlags(bool);

  /** default constructor
   * \param the type of the system
   */
  FirstOrderNonLinearDS(DSTYPES);

public:

  // ===== CONSTRUCTORS =====

  /** xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothSP::DynamicalSystem (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  FirstOrderNonLinearDS(DynamicalSystemXMLSPtr dsXML,
                        SP::NonSmoothDynamicalSystem = SP::NonSmoothDynamicalSystem());

  /** constructor from a set of data
   *  \param int : reference number for this DynamicalSystem
   *  \param SiconosVector : initial state of this DynamicalSystem
   *  \param string : plugin name for f of this DynamicalSystem (optional)
   *  \param string : plugin name for jacobianXF of this DynamicalSystem (optional)
   *  \exception RuntimeException
   */
  FirstOrderNonLinearDS(int, const SiconosVector&, const std::string& = "DefaultPlugin:computeF", const std::string& = "DefaultPlugin:computeJacobianXF");

  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~FirstOrderNonLinearDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

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
  inline SiconosVectorSPtr getRPtr() const
  {
    return r;
  }

  /** set the value of r to newValue
   *  \param SiconosVector newValue
   */
  void setR(const SiconosVector&);

  /** set R to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setRPtr(SiconosVectorSPtr);

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
  inline SiconosMemorySPtr getRMemoryPtr() const
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
  void setRMemoryPtr(SiconosMemorySPtr);

  // --- M ---
  /** get the value of M
  *  \return SimpleMatrix
  */
  inline const SimpleMatrix getMSimple() const
  {
    return *M;
  }

  /** get the value of M
  *  \return BlockMatrix
  */
  inline const BlockMatrix getMBlock() const
  {
    return *M;
  }

  /** get M
  *  \return pointer on a SiconosMatrix
  */
  inline SiconosMatrixSPtr getMPtr() const
  {
    return M;
  }

  /** set the value of M to newValue
  *  \param SimpleMatrix newValue
  */
  void setM(const SiconosMatrix&);

  /** link M with a new pointer
   *  \param a pointer to SiconosMatrix
   */
  void setMPtr(SiconosMatrixSPtr);

  // --- invM ---
  /** get the value of invM
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getInvMSimple() const
  {
    return *invM;
  }

  /** get the value of invM
  *  \return BlockMatrix
  */
  inline const BlockMatrix getInvMBlock() const
  {
    return *invM;
  }

  /** get invM
  *  \return pointer on a SiconosMatrix
  */
  inline SiconosMatrixSPtr getInvMPtr() const
  {
    return invM;
  }

  /** set the value of invM to newValue
  *  \param SimpleMatrix newValue
  */
  void setInvM(const SiconosMatrix&);

  /** link invM with a new pointer
   *  \param a pointer to SiconosMatrix
   */
  void setInvMPtr(SiconosMatrixSPtr);

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
  inline SiconosVectorSPtr getFPtr() const
  {
    return f;
  }

  /** set the value of f to newValue
   *  \param SiconosVector newValue
   */
  virtual void setF(const SiconosVector&);

  /** set f to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  virtual void setFPtr(SiconosVectorSPtr);

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
  inline SiconosMatrixSPtr getJacobianXFPtr() const
  {
    return jacobianXF;
  }

  /** set the value of JacobianXF to newValue
   *  \param SiconosMatrix newValue
   */
  virtual void setJacobianXF(const SiconosMatrix&);

  /** set JacobianXF to pointer newPtr
   *  \param SP::SiconosMatrix  newPtr
   */
  virtual void setJacobianXFPtr(SiconosMatrixSPtr newPtr);

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization
   */
  virtual void initRhs(double) ;

  /** dynamical system initialization function: mainly set memory and compute value for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  void initialize(const std::string&, double = 0, unsigned int = 1) ;

  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** initialize the SiconosMemory objects: reserve memory for i vectors in memory and reset all to zero.
   *  \param the size of the SiconosMemory (i)
   */
  void initMemory(unsigned int) ;

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  void swapInMemory();

  // ===== COMPUTE PLUGINS FUNCTIONS =====

  // --- setters for functions to compute plugins ---

  /** to set a specified function to compute f(x,t)
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeFFunction(const std::string&  pluginPath, const std::string& functionName);

  /** to set a specified function to compute jacobianXF
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeJacobianXFFunction(const std::string&  pluginPath, const std::string&  functionName);

  // --- compute plugin functions ---

  /** Default function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   */
  virtual void computeF(double);

  /** function to compute \f$ f: (x,t)\f$ with x different from current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  virtual void computeF(double, SiconosVectorSPtr);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXF(double, bool  = false);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$ with x different from current saved state.
   *  \param double time : current time
   *  \param SP::SiconosVector
   */
  virtual void computeJacobianXF(double, SP::SiconosVector);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeRhs(double, bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXRhs(double, bool  = false);

  // ===== XML MANAGEMENT FUNCTIONS =====

  /** copy the data common to each system in the XML tree
   */
  virtual void saveSpecificDataToXML();

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const;

  /** set R to zero
   */
  void resetNonSmoothPart();

  /** To compute \f$\frac{|x_{i+1} - xi|}{|x_i|}\f$ where \f$x_{i+1}\f$ represents the present state and \f$x_i\f$ the previous one
   * \return a double
   */
  virtual double dsConvergenceIndicator();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param SP::DynamicalSystem : the system which must be converted
   * \return a pointer on the system if it is of the right type, NULL otherwise
   */
  static FirstOrderNonLinearDS* convert(DynamicalSystem* ds);

};

#endif



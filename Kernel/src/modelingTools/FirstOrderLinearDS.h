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
/*! \file FirstOrderLinearDS.h

*/
#ifndef FOLINEARDS_H
#define FOLINEARDS_H

#include "FirstOrderNonLinearDS.h"

class FirstOrderNonLinearDS;

/** First order linear systems - Inherits from DynamicalSystems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * M \dot x = A(t)x(t)+ b(t) + r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$M \in R^{n\times n} \f$ is an optional constant invertible matrix
 *  The  right-hand side is described by
 *    - \f$A \in R^{n\times n} \f$
 *    - \f$b \in R^{n} \f$
 *
 * Specific members of this class are A and b.
 *
 * f is not set for such system and thus calls to computeF or other related functions are forbidden.
 *
 *  Thus, the main steps for FirstOrderLinearDS handling consist in:
 *
 *  - Construction: A and b are optional, and can be given as a matrix/vector or a plug-in.
 *  - Initialization: compute values at time=t0 (rhs, jacobianXF, A ...), usually done when calling simulation->initialize.
 *  - Computation at time t, by calling "compute" functions
 *      => computeA
 *      => computeB
 *      => computeRhs, compute \f$ \dot x = M^{-1}(Ax + b + r) \f$
 *
 * Any call to a plug-in requires that it has been set correctly before simulation using one of the following:
 *   => setComputeAFunction
 *   => setComputeBFunction
 *
 **/
typedef   void (*bPtrFunction)(double, unsigned int, double*, unsigned int, double*);


class FirstOrderLinearDS : public FirstOrderNonLinearDS
{
protected:

  /** matrix specific to the FirstOrderLinearDS \f$ A \in R^{n \times n}  \f$*/
  SiconosMatrix *A;

  /** strength vector */
  SimpleVector *b;

  /** FirstOrderLinearDS plug-in to compute A(t,z), id = "A"
   * @param time : current time
   * @param sizeOfA : size of square-matrix A
   * @param[in,out] A : pointer to the first element of A
   * @param size of vector z
   * @param[in,out] z a vector of user-defined parameters
   */
  void (*APtr)(double, unsigned int, double*, unsigned int, double*);

  /** FirstOrderLinearDS plug-in to compute b(t,z), id = "b"
   * @param time : current time
   * @param sizeOfB : size of vector b
   * @param[in,out] b : pointer to the first element of b
   * @param size of vector z
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*bPtr)(double, unsigned int, double*, unsigned int, double*);

  /** set all allocation flags (isAllocated map)
   *  \param bool: = if true (default) set default configuration, else set all to false
   */
  void initAllocationFlags(bool  = true);

  /** set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  void initPluginFlags(bool);

  /** default constructor
   */
  FirstOrderLinearDS();

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** xml constructor
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   */
  FirstOrderLinearDS(DynamicalSystemXML *, NonSmoothDynamicalSystem* = NULL);

  /** constructor from a set of data
   *  \param int : reference number of this DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param string: plugin for A (optional)
   *  \param string: plugin for b (optional)
   */
  FirstOrderLinearDS(int, const SiconosVector&, const std::string& = "DefaultPlugin:computeA",
                     const std::string& = "DefaultPlugin:computeB");

  /** constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   */
  FirstOrderLinearDS(int, const SiconosVector&, const SiconosMatrix&);

  /** constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \param SiconosVector : b
   */
  FirstOrderLinearDS(int, const SiconosVector&, const SiconosMatrix&, const SiconosVector&);

  /** destructor */
  virtual ~FirstOrderLinearDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization.
   */
  virtual void initRhs(double) ;

  // --- getter and setter ---

  // --- A ---
  /** get the value of A
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getA() const
  {
    return *A;
  }

  /** get A
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getAPtr() const
  {
    return A;
  }

  /** set the value of A to newValue
   *  \param SiconosMatrix newValue
   */
  void setA(const SiconosMatrix& newValue);

  /** set A to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setAPtr(SiconosMatrix *);

  // --- b ---

  /** get the value of b
   *  \return SimpleVector
   */
  inline const SimpleVector getB() const
  {
    return *b;
  }

  /** get b
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getBPtr() const
  {
    return b;
  }

  /** set the value of b to newValue
   *  \param SimpleVector newValue
   */
  void setB(const SimpleVector&);

  /** set b to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setBPtr(SimpleVector *);

  // --- plugins related functions

  /** set a specified function to compute the matrix A => same action as setComputeJacobianXFFunction
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string& , const std::string&);

  /** set a specified function to compute the vector b
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeBFunction(const std::string& , const std::string&);
  /** set a specified function to compute the vector b
   *  \param bPtrFunction : a pointer on the plugin function
   */
  void setComputeBFunction(bPtrFunction fct);

  /** default function to compute matrix A => same action as computeJacobianXF
   */
  void computeA(double);

  /** default function to compute vector b
   */
  void computeB(double);

  /** set the value of f to newValue
   *  \param SiconosVector newValue
   */
  inline void setF(const SiconosVector&)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - setF: f is not available for FirstOrderLinearDS.");
  };

  /** set f to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  inline void setFPtr(SiconosVector *)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - setFPtr: f is not available for FirstOrderLinearDS.");
  };

  /** set the value of JacobianXF to newValue: exception for LinearDS since f is not available.
   *  \param SiconosMatrix newValue
   */
  inline void setJacobianXF(const SiconosMatrix&)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - setJacobianXF: f is not available for FirstOrderLinearDS.");
  };

  /** set JacobianXF to pointer newPtr: exception for LinearDS since f is not available.
   *  \param SiconosMatrix * newPtr
   */
  inline void setJacobianXFPtr(SiconosMatrix *newPtr)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - setJacobianXFPtr: f is not available for FirstOrderLinearDS.");
  };

  /** to set a specified function to compute f(x,t): exception for LinearDS since f is not available.
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   */
  inline void setComputeFFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - setComputeFFunction: f is not available for FirstOrderLinearDS.");
  };

  /** to set a specified function to compute jacobianXF: exception for LinearDS since f is not available.
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   */
  inline void setComputeJacobianXFFunction(const std::string&  pluginPath, const std::string&  functionName)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - setComputeJacobianXFFunction: f is not available for FirstOrderLinearDS.");
  };

  // --- compute plugin functions ---

  /** Default function to compute \f$ f: (x,t)\f$: exception for LinearDS since f is not available.
   * \param double time : current time
   */
  inline void computeF(double)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - computeF: f is not available for FirstOrderLinearDS.");
  };

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$: exception for LinearDS since f is not available.
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  inline void computeJacobianXF(double, bool  = false)
  {
    RuntimeException::selfThrow("FirstOrderLinearDS - computeJacobianXF: f is not available for FirstOrderLinearDS.");
  };

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeRhs(double, bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeJacobianXRhs(double, bool  = false);

  // --- xml related functions ---

  /** copy the data specific to each system into the XML tree
   */
  void saveSpecificDataToXML();

  /** data display on screen
   */
  virtual void display() const;

  /** overload LagrangianDS corresponding function
   * \return a double, always zero.
   */
  double dsConvergenceIndicator()
  {
    return 0.0;
  }

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
   */
  static FirstOrderLinearDS* convert(DynamicalSystem* ds);

};

#endif // FOLINEARDS_H

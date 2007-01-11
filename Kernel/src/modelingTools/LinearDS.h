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
/*! \file LinearDS.h

*/
#ifndef LINEARDS_H
#define LINEARDS_H

#include "LinearDSXML.h"
#include "DynamicalSystem.h"

class LinearDSXML;

/** First order linear systems - Inherits from DynamicalSystems
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * M \dot x = A(t)x(t)+T u(t)+b(t)+r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$Mxdot \in R^{n\times n} \f$ is an optional constant invertible matrix
 *  The  right-hand side is described by
 *    - \f$A \in R^{n\times n} \f$
 *    - \f$b \in R^{n} \f$
 *    - \f$u \in R^{uSize} \f$
 *    - \f$T \in R^{n\times uSize} \f$
 *        warning: T and u are members of DynamicalSystem class.
 *
 * The "minimal" form is
 * \f[
 * \dot x = A(t)x(t),
 *  x(t_0)=x_0
 * \f]
 * and so A should always be specified.
 *
 * Links with DynamicalSystem are:
 *
 * \f[
 *   f(x,t) = A(t)x(t) + b(t) \\
 *   jacobianXF = A(t)
 * \f]
 *
 *  Thus, the main steps for LinearDS handling consist in:
 *
 *  - Construction: A is required and must be set as a matrix or a plug-in. b is optional, and can be given as a vector or a plug-in.
 *  - Initialization: compute values at time=t0 (rhs, jacobianXF, A ...), usually done when calling simulation->initialize.
 *  - Computation at time t, by calling "compute" functions
 *      => computeA
 *      => computeB
 *      => computeRhs, compute \f$ \dot x = M^{-1}(Ax + b + Tu + r) \f$
 *      => computeU and computeT (from DynamicalSystem class)
 *
 * Any call to a plug-in requires that it has been set correctly before simulation using one of the following:
 *   => setComputeAFunction
 *   => setComputeBFunction
 *
 **/
class LinearDS : public DynamicalSystem
{
protected:

  /** matrix specific to the LinearDS \f$ A \in R^{n \times n}  \f$*/
  SiconosMatrix *A;

  /** strength vector */
  SimpleVector *b;

  /* the name of the plugin used to compute A */
  std::string  computeAFunctionName;

  /* the name of the plugin used to compute b */
  std::string  computeBFunctionName;

  /** LinearDS plug-in to compute A(t), id = "A"
   * @param sizeOfA : size of square-matrix A
   * @param time : current time
   * @param[in,out] A : pointer to the first element of A
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*APtr)(unsigned int, double, double*, double*);

  /** LinearDS plug-in to compute b(t), id = "b"
   * @param sizeOfB : size of vector b
   * @param time : current time
   * @param[in,out] b : pointer to the first element of b
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*bPtr)(unsigned int, double, double*, double*);

  /** set all allocation flags (isAllocated map)
  *  \param bool: = if true (default) set default configuration, else set all to false
  */
  virtual void initAllocationFlags(const bool  = true);

  /** set all plug-in flags (isPlugin map) to val
  *  \param a bool
  */
  virtual void initPluginFlags(const bool);

  /** default constructor
  */
  LinearDS();

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** xml constructor
  *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
  *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
  *  \exception RuntimeException
  */
  LinearDS(DynamicalSystemXML *, NonSmoothDynamicalSystem* = NULL);

  /** constructor from a set of data
  *  \param int : reference number of this DynamicalSystem
  *  \param int : dimension of this DynamicalSystem
  *  \param SiconosVector : the initial state of this DynamicalSystem
  *  \param string: plugin for A (optional)
  *  \param string: plugin for b (optional)
  *  \exception RuntimeException
  */
  LinearDS(const int, const unsigned int, const SiconosVector&, const std::string = "DefaultPlugin:computeA",
           const std::string = "DefaultPlugin:computeB");

  /** constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \exception RuntimeException
   */
  LinearDS(const int, const SiconosVector&, const SiconosMatrix&);

  /** constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \param SiconosVector : b
   *  \exception RuntimeException
   */
  LinearDS(const int, const SiconosVector&, const SiconosMatrix&, const SiconosVector&);

  /** copy constructor
  *  \param a Dynamical system to copy
  */
  LinearDS(const LinearDS &);

  /** copy constructor
  *  \param a Dynamical system to copy
  */
  LinearDS(const DynamicalSystem &);

  /** destructor */
  virtual ~LinearDS();

  /** check that the system is complete (ie all required data are well set)
  * \return a bool
  */
  virtual bool checkDynamicalSystem();

  /** dynamical system initialization function: mainly set memory and compute value for initial state values.
  *  \param string: simulation type
  *  \param time of initialisation, default value = 0
  *  \param the size of the memory, default size = 1.
  */
  void initialize(const std::string&, double = 0, unsigned int = 1) ;

  /** dynamical system update: mainly call compute for all time or state depending functions
  *  \param current time
  */
  void update(const double);

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

  /** set the value of JacobianXF to newValue
  *  \param SiconosMatrix newValue
  */
  void setJacobianXF(const SiconosMatrix&);

  /** set JacobianXF to pointer newPtr
  *  \param SiconosMatrix * newPtr
  */
  void setJacobianXFPtr(SiconosMatrix *newPtr);

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

  /** get name of function that computes A = jacobianXF
  *  \return a string
  */
  inline const std::string getComputeAFunctionName() const
  {
    return computeAFunctionName;
  }

  /** set a specified function to compute the matrix A => same action as setComputeJacobianXFFunction
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeAFunction(const std::string , const std::string);

  /** get name of function that computes b (if b from plugin)
  *  \return a string
  */
  inline const std::string getComputeBFunctionName() const
  {
    return computeBFunctionName;
  }

  /** set a specified function to compute the vector b
  *  \param string : the complete path to the plugin
  *  \param string : the function name to use in this plugin
  *  \exception SiconosSharedLibraryException
  */
  virtual void setComputeBFunction(const std::string , const std::string);

  /** default function to compute matrix A => same action as computeJacobianXF
  *  \exception RuntimeException
  */
  void computeA(const double);

  /** default function to compute vector b
  *  \exception RuntimeException
  */
  void computeB(const double);

  /** Default function to compute \f$ f: (x,t)\f$
  * \param double time : current time
  *  \exception RuntimeException
  */
  virtual void computeF(const double);

  /** Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
  *  \param double time : current time
  *  \param bool isDSup : flag to avoid recomputation of operators
  *  \exception RuntimeException
  */
  virtual void computeJacobianXF(const double, const bool  = false);

  /** Default function to the right-hand side term
  *  \param double time : current time
  *  \param bool isDSup : flag to avoid recomputation of operators
  *  \exception RuntimeException
  */
  virtual void computeRhs(const double, const bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
  *  \param double time : current time
  *  \param bool isDSup : flag to avoid recomputation of operators
  *  \exception RuntimeException
  */
  virtual void computeJacobianXRhs(const double, const bool  = false);

  // --- xml related functions ---

  /** copy the data of the DS into the XML tree
  *  \exception RuntimeException
  */
  void saveDSToXML();

  /** data display on screen
  */
  virtual void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
  *  \param DynamicalSystem* : the system which must be converted
  * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
  */
  static LinearDS* convert(DynamicalSystem* ds);

};

#endif // LINEARDS_H

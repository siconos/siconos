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

typedef void (*bPtrFunction)(double, unsigned int, double*, unsigned int, double*);

typedef PluggedObject<bPtrFunction, SimpleVector> PVB;
typedef PluggedObject<bPtrFunction, SimpleMatrix> PMA;

TYPEDEF_SPTR(PVB);
TYPEDEF_SPTR(PMA);

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
class FirstOrderLinearDS : public FirstOrderNonLinearDS
{
protected:

  /** matrix specific to the FirstOrderLinearDS \f$ A \in R^{n \times n}  \f$*/
  SP::PMA A;

  /** strength vector */
  SP::PVB b;

  /** default constructor
   */
  FirstOrderLinearDS(): FirstOrderNonLinearDS(DS::FOLDS) {};

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** xml constructor
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   */
  FirstOrderLinearDS(SP::DynamicalSystemXML);

  /** constructor from a set of data
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param string: plugin for A
   *  \param string: plugin for b
   */
  FirstOrderLinearDS(const SiconosVector&, const std::string&, const std::string&);

  /** constructor from a set of data
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   */
  FirstOrderLinearDS(const SiconosVector&, const SiconosMatrix&);

  /** constructor from a set of data
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \param SiconosVector : b
   */
  FirstOrderLinearDS(const SiconosVector&, const SiconosMatrix&, const SiconosVector&);

  /** destructor */
  virtual ~FirstOrderLinearDS() {};

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization.
   */
  void initRhs(double) ;

  // --- getter and setter ---

  // --- A ---
  /** get the value of A
   *  \return a plugged-matrix
   */
  inline const PMA getA() const
  {
    return *A;
  }

  /** get A
   *  \return pointer on a plugged-matrix
   */
  inline SP::PMA getAPtr() const
  {
    return A;
  }

  /** set the value of A to newValue
   *  \param plugged-matrix newValue
   */
  void setA(const PMA&);

  /** set A to pointer newPtr
   *  \param a plugged matrix SP
   */
  inline void setAPtr(SP::PMA newPtr)
  {
    A = newPtr;
  }

  // --- b ---

  /** get the value of b
   *  \return plugged vector
   */
  inline const PVB getB() const
  {
    return *b;
  }

  /** get b
   *  \return pointer on a plugged vector
   */
  inline SP::PVB getBPtr() const
  {
    return b;
  }

  /** set the value of b to newValue
   *  \param a plugged vector
   */
  void setB(const PVB&);

  /** set b to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setBPtr(SP::PVB newPtr)
  {
    b = newPtr;
  }

  // --- plugins related functions

  /** set a specified function to compute the matrix A => same action as setComputeJacobianXFFunction
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string& , const std::string&);

  /** set a specified function to compute the matrix A
   *  \param bPtrFunction : a pointer on the plugin function
   */
  void setComputeAFunction(bPtrFunction fct);

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

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  void computeRhs(double, bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  void computeJacobianXRhs(double, bool  = false);

  // --- xml related functions ---

  /** copy the data specific to each system into the XML tree
   */
  void saveSpecificDataToXML();

  /** data display on screen
   */
  void display() const;

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

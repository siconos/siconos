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

#include "FirstOrderNonLinearDS.hpp"

class FirstOrderNonLinearDS;
typedef   void (*LDSPtrFunction)(double, unsigned int, double*, unsigned int, double*);


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
 * f is not set for such system and thus calls to computeF or other
 * related functions are forbidden.
 *
 *  Thus, the main steps for FirstOrderLinearDS handling consist in:
 *
 * - Construction: A and b are optional, and can be given as a
 *   matrix/vector or a plug-in.
 * - Initialization: compute values at time=t0 (rhs, jacobianfx, A
 *    ...), usually done when calling simulation->initialize.
 * - Computation at time t, by calling "compute" functions
 *      => computeA
 *      => computeb
 *      => computeRhs, compute \f$ \dot x = M^{-1}(Ax + b + r) \f$
 *
 * Any call to a plug-in requires that it has been set correctly
 * before simulation using one of the following:
 *   => setComputeAFunction
 *   => setComputeBFunction
 *
 **/
class FirstOrderLinearDS : public FirstOrderNonLinearDS
{
protected:

  /** matrix specific to the FirstOrderLinearDS \f$ A \in R^{n \times n}  \f$*/
  SP::SiconosMatrix _A;

  /** strength vector */
  SP::SiconosVector _b;

  /** FirstOrderLinearDS plug-in to compute A(t,z), id = "A"
  * @param time : current time
  * @param sizeOfA : size of square-matrix A
  * @param[in,out] A : pointer to the first element of A
  * @param size of vector z
  * @param[in,out] z a vector of user-defined parameters
  */
  //  LDSPtrFunction _APtr;
  //  std::string _pluginNameAPtr;
  SP::PluggedObject _pluginA;

  /** FirstOrderLinearDS plug-in to compute b(t,z), id = "b"
  * @param time : current time
  * @param sizeOfB : size of vector b
  * @param[in,out] b : pointer to the first element of b
  * @param size of vector z
  * @param[in,out] param  : a vector of user-defined parameters
  */
  //  LDSPtrFunction _bPtr;
  //  std::string _pluginNamebPtr;
  SP::PluggedObject _pluginb;

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
  FirstOrderLinearDS(SP::SiconosVector, const std::string&, const std::string&);
  FirstOrderLinearDS(const SiconosVector&, const std::string&, const std::string&);

  /** constructor from a set of data
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   */
  FirstOrderLinearDS(SP::SiconosVector, SP::SiconosMatrix);

  /** constructor from a set of data
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \param SiconosVector : b
   */
  FirstOrderLinearDS(SP::SiconosVector, SP::SiconosMatrix, SP::SiconosVector);

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

  /** Call all plugged-function to initialize plugged-object values
      \param time
   */
  virtual void updatePlugins(double);

  // --- getter and setter ---

  // --- A ---
  /** get the value of A
   *  \return a plugged-matrix

  inline const Plugged_Matrix_FTime getA() const { return *A; }
  */
  /** get A
   *  \return pointer on a plugged-matrix
   */
  inline SP::SiconosMatrix A() const
  {
    return _A;
  }

  virtual SP::SiconosMatrix jacobianfx() const
  {
    return _A;
  };
  /**  function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   */
  virtual void computeF(double);

  /** function to compute \f$ f: (x,t)\f$ with x different from
      current saved state.
   * \param double time : current time
   * \param SP::SiconosVector
   */
  virtual void computeF(double, SP::SiconosVector);


  /** set the value of A to newValue
   *  \param plugged-matrix newValue

  void setA(const Plugged_Matrix_FTime&);
  */

  /** set A to pointer newPtr
   *  \param a plugged matrix SP
   */
  inline void setAPtr(SP::SiconosMatrix newPtr)
  {
    _A = newPtr;
  }
  inline void setA(SP::SiconosMatrix newPtr)
  {
    _A = newPtr;
  }

  // --- b ---

  /** get the value of b
   *  \return plugged vector

  inline const Plugged_Vector_FTime getB() const { return *b; }
  */
  /** get b
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector b() const
  {
    return _b;
  }

  /** set the value of b to newValue
   *  \param a plugged vector

  void setB(const Plugged_Vector_FTime&);
  */

  /** set b to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setb(SP::SiconosVector newPtr)
  {
    _b = newPtr;
  }

  // --- plugins related functions

  /** set a specified function to compute the matrix A => same action as setComputeJacobianfxFunction
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string& , const std::string&);

  /** set a specified function to compute the matrix A
   *  \param VectorFunctionOfTime : a pointer on the plugin function
   */
  void setComputeAFunction(LDSPtrFunction fct);

  /** set a specified function to compute the vector b
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputebFunction(const std::string& , const std::string&);

  /** set a specified function to compute the vector b
   *  \param VectorFunctionOfTime : a pointer on the plugin function
   */
  void setComputebFunction(LDSPtrFunction fct);

  /** default function to compute matrix A => same action as
      computeJacobianfx
   */
  void computeA(double);

  /** default function to compute vector b
   */
  void computeb(double);

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  void computeRhs(double, bool  = false);

  /** Default function to jacobian of the right-hand side term
      according to x

   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  void computeJacobianRhsx(double, bool  = false);

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
    return 1.0;
  }

  /** encapsulates an operation of dynamic casting. Needed by Python
      interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right
   * type, NULL otherwise
   */
  static FirstOrderLinearDS* convert(DynamicalSystem* ds);

};

TYPEDEF_SPTR(FirstOrderLinearDS);

#endif // FOLINEARDS_H

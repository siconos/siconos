/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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
#ifndef LINEARTIDS_H
#define LINEARTIDS_H

#include "LinearDSXML.h"
#include "LinearDS.h"

class LinearDSXML;

/** \class LinearTIDS
 *  \brief First order linear systems - Inherits from DynamicalSystems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * \dot x = Ax(t)+T u(t)+b+r,
 *  x(t_0)=x_0
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *
 *  The  VectorField is specialized by
 *    - \f$A \in R^{n\times n} \f$
 *    - \f$b \in R^{n} \f$
 *    - \f$u \in R^{uSize} \f$
 *    - \f$T \in R^{n\times uSize} \f$
 *        warning: T and u are members of DynamicalSystem class.
 *
 * The "minimal" form is
 * \f[
 * \dot x = Ax(t),
 *  x(t_0)=x_0
 * \f]
 * and so A should always be specified.
 *
 * This class inherits from LinearDS one. The difference is that here A and b do not depends on time.
 *
 * Links with DynamicalSystem are:
 *
 * \f[
 *   f(x,t) = Ax(t) + b \\
 *   jacobianX = A
 * \f]
 *
 * A and b are constant matrix or vector, and thus can not be set using a plug-in.
 * Setting a value to A, automatically set jacobianX.
 *
 * To build and use a linearDS, you first need to call a constructor, with A as a required data.
 * Then, this system has to be initialized -> compute members value at time t0. This is usually done during call to strategy->initialize.
 * Finally, the state of the DS can be obtained by calling "compute" functions. In LinearTIDS case, since A and b are fixed, you can
 * only call computeVectorField(time), to compute xDot = Ax+b+Tu.
 *
 **/

class LinearTIDS : public LinearDS
{
private:

  /** \fn LinearTIDS()
   *  \brief default constructor
   */
  LinearTIDS();

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** \fn LinearTIDS(DynamicalSystemXML * nsdsXML)
   *  \brief xml constructor
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  LinearTIDS(DynamicalSystemXML *, NonSmoothDynamicalSystem* = NULL);

  /** \fn LinearTIDS(const int& ref, const unsigned int& n, const SiconosVector& x0, const SiconosMatrix& A);
   *  \brief constructor from a set of data
   *  \param int : reference number of this DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix: A
   *  \exception RuntimeException
   */
  LinearTIDS(const int&, const SiconosVector&, const SiconosMatrix&);

  /** \fn LinearTIDS(const int& ref, const unsigned int& n, const SiconosVector& x0, const SiconosMatrix& A, const SiconosVector& b);
   *  \brief constructor from a set of data
   *  \param int : reference number of this DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix: A
   *  \param SiconosVector: b
   *  \exception RuntimeException
   */
  LinearTIDS(const int&, const SiconosVector&, const SiconosMatrix&, const SiconosVector&);

  /** \fn LinearTIDS(const LinearTIDS &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  LinearTIDS(const LinearTIDS &);

  /** \fn LinearTIDS(const DynamicalSystem &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  LinearTIDS(const DynamicalSystem &);

  /** \fn ~LinearTIDS()
   *  \brief destructor */
  ~LinearTIDS();

  /** \fn bool checkDynamicalSystem()
   *  \brief check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /** \fn void initialize(const double& = 0, const unsigned int& = 1) ;
   *  \brief dynamical system initialization function: mainly set memory and compute value for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  //  void initialize(const double& = 0, const unsigned int& = 1) ;

  /** \fn void setComputeJacobianXFunction(const string&, const string&)
   *  \brief set a specified function to compute jacobianX=A
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianXFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeAFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the matrix A => same action as setComputeJacobianXFunction
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeAFunction(const std::string &, const std::string &);

  /** \fn void setComputeBFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the vector b
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeBFunction(const std::string &, const std::string &);

  /** \fn void vectorField (const double& time)
   * \brief compute the vector field Ax+b
   * \param double time : current time
   *  \exception RuntimeException
   */
  void computeVectorField(const double&);

  /** \fn void display()
   *  \brief data display on screen
   */
  void display() const;

  /** \fn LinearTIDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
   */
  static LinearTIDS* convert(DynamicalSystem* ds);

};

#endif // LINEARTIDS_H

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
#ifndef LINEARDS_H
#define LINEARDS_H

#include "LinearDSXML.h"
#include "DynamicalSystem.h"

class LinearDSXML;

/** \class LinearDS
 *  \brief First order linear systems - Inherits from DynamicalSystems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Creation) Apr 29, 2004
 *
 *
 *  This class represents first order linear systems of the form:
 * \f[
 * \dot x = A(t)x(t)+T u(t)+b(t)+r,
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
 * \dot x = A(t)x(t),
 *  x(t_0)=x_0
 * \f]
 * and so A should always be specified.
 *
 * Links with DynamicalSystem are:
 *
 * \f[
 *   f(x,t) = A(t)x(t) + b(t) \\
 *   jacobianX = A(t)
 * \f]
 *
 *  Thus, the main steps for LinearDS handling consist in:
 *
 *  - Construction: A is required and must be set as a matrix or a plug-in. b is optional, and can be given as a vector or a plug-in.
 *  - Initialization: compute values at time=t0 (vectorField, jacobianX, A ...), usually done when calling strategy->initialize.
 *  - Computation a time t, by calling "compute" functions
 *      => computeA or computeJacobianX, same result.
 *      => computeB
 *      => computeVectorField, compute xDot = Ax + b + Tu
 *      => computeU and computeT (from DynamicalSystem class)
 *
 * Any call to a plug-in requires that it has been set correctly before simulation using one of the following:
 *   => setComputeJacobianXFunction or setComputeA (same result)
 *   => setComputeBFunction
 *   Warning: setVectorFieldFunction is useless in LinearDS case, since xDot is computed using A, b,u and T.
 *
 **/

class LinearDS : public DynamicalSystem
{
protected:

  /** matrix specific to the LinearDS \f$ A \in R^{n \times n}  \f$*/
  SiconosMatrix *A;
  /** strength vector */
  SimpleVector *b;

  /* the name of the plugin used to compute b */
  std::string  bFunctionName;

  /** \fn void (*bPtr) (const unsigned int & sizeOfB, const double* t, double* b, double* param)
   *  \brief pointer on function to compute b
   *  \param unsigned int sizeOfB : size of vector b
   *  \param double* time : current time
   *  \param double* b : the pointer to the first element of the vector b
   *    \param double* param   : a vector of user-defined parameters
   */
  void (*computeBPtr)(const unsigned int &, const double*, double*, double*);

  /** \fn LinearDS()
   *  \brief default constructor
   */
  LinearDS();

public:

  /** === CONSTRUCTORS/DESTRUCTOR === */

  /** \fn LinearDS(DynamicalSystemXML * nsdsXML)
   *  \brief xml constructor
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  LinearDS(DynamicalSystemXML *, NonSmoothDynamicalSystem* = NULL);

  /** \fn LinearDS(int number, int n, SiconosVector* x0, NSDS * nsds)
   *  \brief constructor from a set of data
   *  \param int : reference number of this DynamicalSystem
   *  \param int : dimension of this DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param string: plugin for A=jacobianX (optional)
   *  \param string: plugin for b (optional)
   *  \exception RuntimeException
   */
  LinearDS(const int&, const unsigned int&, const SiconosVector&, const std::string& = "DefaultPlugin:jacobianX",
           const std::string& = "DefaultPlugin:computeB");

  /** \fn LinearDS( const int& newNumber, const SiconosVector& newX0,
   *                const SiconosMatrix& newA)
   *  \brief constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \exception RuntimeException
   */
  LinearDS(const int&, const SiconosVector&, const SiconosMatrix&);

  /** \fn LinearDS( const int& newNumber, const SiconosVector& newX0,
   *                const SiconosMatrix& newA, const SiconosVector& newB)
   *  \brief constructor from a set of data
   *  \param int : reference number of the DynamicalSystem
   *  \param SiconosVector : the initial state of this DynamicalSystem
   *  \param SiconosMatrix : matrix A
   *  \param SiconosVector : b
   *  \exception RuntimeException
   */
  LinearDS(const int&, const SiconosVector&, const SiconosMatrix&, const SiconosVector&);

  /** \fn LinearDS(const LinearDS &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  LinearDS(const LinearDS &);

  /** \fn LinearDS(const DynamicalSystem &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  LinearDS(const DynamicalSystem &);

  /** \fn ~LinearDS()
   *  \brief destructor */
  virtual ~LinearDS();

  /** \fn bool checkDynamicalSystem()
   *  \brief check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

  /** \fn void initialize(const double& = 0, const unsigned int& = 1) ;
   *  \brief dynamical system initialization function: mainly set memory and compute value for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const double& = 0, const unsigned int& = 1) ;

  // --- getter and setter ---

  // --- A ---
  /** \fn  const SimpleMatrix getA() const
   *  \brief get the value of A
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getA() const
  {
    return *A;
  }

  /** \fn SiconosMatrix* getAPtr() const
   *  \brief get A
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getAPtr() const
  {
    return A;
  }

  /** \fn void setA (const SiconosMatrix& newValue)
   *  \brief set the value of A to newValue
   *  \param SiconosMatrix newValue
   */
  void setA(const SiconosMatrix& newValue);

  /** \fn void setAPtr(SiconosMatrix* newPtr)
   *  \brief set A to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setAPtr(SiconosMatrix *);

  /** \fn void setJacobianX (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianX to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianX(const SiconosMatrix&);

  /** \fn void setJacobianXPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianX to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianXPtr(SiconosMatrix *newPtr);

  // --- b ---

  /** \fn  const SimpleVector getB() const
   *  \brief get the value of b
   *  \return SimpleVector
   */
  inline const SimpleVector getB() const
  {
    return *b;
  }

  /** \fn SimpleVector* getBPtr() const
   *  \brief get b
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getBPtr() const
  {
    return b;
  }

  /** \fn void setB (const SimpleVector& newValue)
   *  \brief set the value of b to newValue
   *  \param SimpleVector newValue
   */
  void setB(const SimpleVector&);

  /** \fn void setBPtr(SimpleVector* newPtr)
   *  \brief set b to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setBPtr(SimpleVector *);

  // --- plugins related functions

  /** \fn  std::string getAFunctionName() const
   *  \brief get name of function that computes A = jacobianX
   *  \return a string
   */
  inline const std::string getAFunctionName() const
  {
    return computeJacobianXFunctionName;
  }

  /** \fn void setComputeAFunction(const string& libPath,const string& functionName)
   *  \brief set a specified function to compute the matrix A => same action as setComputeJacobianXFunction
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeAFunction(const std::string &, const std::string &);

  /** \fn  std::string getBFunctionName() const
  *  \brief get name of function that computes b (if b from plugin)
  *  \return a string
  */
  inline const std::string getBFunctionName() const
  {
    return bFunctionName;
  }

  /** \fn void setComputeBFunction(const string& libPath,const string& functionName);
   *  \brief set a specified function to compute the vector b
   *  \param string : the complete path to the plugin
   *  \param string : the function name to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeBFunction(const std::string &, const std::string &);

  /** \fn void setVectorFieldFunction(const string&, const string&)
   *  \brief overload corresponding function of DS -> did nothing.
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setVectorFieldFunction(const std::string & pluginPath, const std::string& functionName);

  /** \fn void setComputeJacobianXFunction(const string&, const string&)
   *  \brief set a specified function to compute jacobianX=A
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeJacobianXFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void vectorField (const double& time)
   * \brief compute the vector field Ax+b
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeVectorField(const double&);

  /** \fn void computeA(const double& time)
   *  \brief default function to compute matrix A => same action as computeJacobianX
   *  \exception RuntimeException
   */
  void computeA(const double&);

  /** \fn void computeB(const double& time)
   *  \brief default function to compute vector b
   *  \exception RuntimeException
   */
  void computeB(const double&);

  // --- xml related functions ---

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS into the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief data display on screen
   */
  virtual void display() const;

  /** \fn LinearDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the dynamical system if it is of the right type, NULL otherwise
   */
  static LinearDS* convert(DynamicalSystem* ds);

};

#endif // LINEARDS_H

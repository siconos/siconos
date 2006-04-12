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
#ifndef LAGRANGIANNLDS_H
#define LAGRANGIANNLDS_H
#include "DynamicalSystem.h"
#include "LagrangianDSXML.h"
#include "BlockMatrix.h"
#include<map>

class LagrangianDSXML;

/** \class LagrangianDS
 *  \brief Lagrangian non linear dynamical systems - Derived from DynamicalSystem -
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianDS  defines  and computes a generic ndof-dimensional
 * Lagrangian Non Linear Dynamical System of the form :
 * \f[
 * M(q) \ddot q + NNL(\dot q, q) + F_{Int}(\dot q , q , t) = F_{Ext}(t) + p
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In the particular case of Non Smooth evolution, the variable p contains the impulse and not the force.
 *    -  \f$ M(q)  \in  R^{ndof \times ndof}  \f$ is the inertia term saved in the SiconosMatrix mass.
 *    -  \f$ NNL(\dot q, q)  \in R^{ndof}\f$ is the non linear inertia term saved in the SiconosVector NNL.
 *    -  \f$ F_{Int}(\dot q , q , t)  \in R^{ndof} \f$ are the internal forces saved in the SiconosVector fInt.
 *    -  \f$ F_{Ext}(t)  \in R^{ndof}  \f$ are the external forces saved in the SiconosVector fExt.
 *
 *
 * Links with first order DynamicalSystem top-class are:
 *
 * \f$ n= 2 ndof \f$
 * \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
 *
 * The rhs is given by:
 * \f[
 * f(x,t) = \left[\begin{array}{c}
 *  \dot q  \\
 * M^{-1}(q)\left[F_{Ext}( q , t) - F_{Int}(\dot q , q , t) - NNL(\dot q, q) \right]\\
 * \end{array}\right]
 * \f]
 * Its jacobian is:
 * \f[
 * \nabla_{x}f(x,t) = \left[\begin{array}{cc}
 *  0  & I \\
 * \nabla_{q}M^{-1}(q)\left[F_{Ext}( q , t) - F_{Int}(\dot q , q , t) - NNL(\dot q, q) \right]-M^{-1}(q)\left[\nabla_q(F_{Int}(\dot q , q , t)-NNL(\dot q, q))\right] & -M^{-1}(q)\left[\nabla_{\dot q}(F_{Int}(\dot q , q , t)-NNL(\dot q, q))\right] \\
 * \end{array}\right]
 * \f]
 *  The input due to the non smooth law is:
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
 *
 * The control part (NOT IMPLEMENTED) is given by:
 *
 *  \f$ u(x,t) = u_L(q,t) \f$
 *
 * \f$ T(x) = \left[\begin{array}{c} 0_{ndof} \\ T_L(q)\end{array}\right]\f$
 *
 *  Main functionalities to handle a LagrangianDS are:
 *
 *    - Construction: the only required operator is M. All the operators can be set using the plug-in mechanism.
 *    - Initialization: compute state members and operators for time=t0 (usually done when calling strategy->initialize)
 *    - Computation at time t, thanks to "compute" functions. Any call to one of the following functions requires that the plug-in
 *      has been set properly thanks to the corresponding setPluginFunction:
 *        => computeMass     (setComputeMassFunction)
 *        => computeFInt     (setComputeFIntFunction)
 *        => computeFExt     (setComputeFExtFunction)
 *        => computeNNL      (setComputeNNLFunction)
 *        => computeJacobianQFInt         (setComputeJacobianQFIntFunction)
 *        => computeJacobianVelocityFInt  (setComputeJacobianVelocityFIntFunction)
 *        => computeJacobianQNNL          (setComputeJacobianQNNLFunction)
 *        => computeJacobianVelocityNNL   (setComputeJacobianVelocityNNLFunction)
 *        => computeRhs            (not set function)
 *        => computeJacobianXRhs   (not set function)
 *
 *  \todo: add a check function (for example check that NNL plug-in given => jacobian required  ...)
 */
class LagrangianDS : public DynamicalSystem
{
protected:

  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  unsigned int ndof;

  /** \var q
      coordinates of the system */
  SimpleVector *q;
  /** \var q0
      initial coordinates of the system */
  SimpleVector *q0;
  /** \var qFree
      free state of the system */
  SimpleVector *qFree;
  /** \var qMemory
      memory of previous coordinates of the system */
  SiconosMemory *qMemory;
  SimpleVector *velocity;/**<velocity of the system*/
  SimpleVector *velocity0; /**< initial velocity of the system*/
  SimpleVector *velocityFree;/**<free Velocity*/
  SiconosMemory *velocityMemory;/**<memory of previous velocities of the system */
  SimpleVector *p;/**<Reaction due to the non smooth law */

  /** mass of the system */
  SiconosMatrix *mass;
  /** internal strength of the system */
  SimpleVector *fInt;
  /** external strength of the system */
  SimpleVector *fExt;
  /** non-linear inertia term of the system */
  SimpleVector *NNL;
  /** jacobian/coordinates of internal strength */
  SiconosMatrix *jacobianQFInt;
  /** jacobian/velocity of internal strength */
  SiconosMatrix *jacobianVelocityFInt;
  /** jacobian/coordinates of inertia */
  SiconosMatrix *jacobianQNNL;
  /** jacobian/velocity of inertie */
  SiconosMatrix *jacobianVelocityNNL;

  /* contains the name of the plugin used to compute the mass */
  std::string  massFunctionName;
  /* contains the name of the plugin used to compute fInt */
  std::string  fIntFunctionName;
  /* contains the name of the plugin used to compute fExt */
  std::string  fExtFunctionName;
  /* contains the name of the plugin used to compute jacobianQFInt */
  std::string  jacobianQFIntFunctionName;
  /* contains the name of the plugin used to compute jacobianQNNL */
  std::string  jacobianQNNLFunctionName;
  /* contains the name of the plugin used to compute jacobianVelocityFInt */
  std::string  jacobianVelocityFIntFunctionName;
  /* contains the name of the plugin used to compute jacobianVelocityNNL */
  std::string  jacobianVelocityNNLFunctionName;
  /* contains the name of the plugin used to compute NNL */
  std::string  NNLFunctionName;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  // pointers to functions member to compute plug-in functions

  /** \fn void (*computeMassPtr)(const unsigned int& sizeOfq, double* qPtr, double* massPtr, double* param)
   *  \brief compute the mass
   *  \param unsigned int& sizeOfq : the size of the vector q
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* massPtr : the pointer to the first element of the matrix mass (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeMassPtr)(const unsigned int&, const double*, double*, double*);

  /** \fn void (*computeFIntPtr)(const unsigned int& sizeOfq, double* time, double* qPtr, double* velocityPtr, double* fIntPtr, double* param)
   *  \brief computes the internal strengths
   *  \param const unsigned int sizeOfq : the size of the vector q
   *    \param int* time : the current time
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* fIntPtr : the pointer to the first element of the vector FInt (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeFIntPtr)(const unsigned int&, const double*, const double*, const double*, double*, double*);

  /** \fn void (*computeFExtPtr)(const unsigned int& sizeOfq, const double* time, double* fExtPtr, double* param)
   *  \brief computes the external strengths
   *  \param unsigned int& sizeOfq : the size of the vector q
   *    \param double* time : the current time
   *  \param double* fExtPtr : the pointer to the first element of the vector FExt (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeFExtPtr)(const unsigned int&, const double*, double*, double*);

  /** \fn void (*computeNNLPtr)(const unsigned int& sizeOfq, const double* qPtr, const double* velocityPtr, double* NNLPtr, double* param)
   *  \brief computes the inertia
   *  \param unsigned int sizeOfq : the size of the vector q
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* NNLPtr : the pointer to the first element of the vector NNL (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeNNLPtr)(const unsigned int&, const double*, const double*, double*, double*);

  /** \fn void (*computeJacobianQFIntPtr)(const unsigned int& sizeOfq, const double* time, const double* qPtr, const double* velocityPtr, double* jacobPtr, double* param)
   *  \brief computes the gradient of the the internal strength compared to the state
   *  \param const unsigned int& sizeOfq : the size of the vector q
   *    \param double* time : the current time
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* jacobPtr : the pointer to the first element of the matrix JacobianQFInt (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeJacobianQFIntPtr)(const unsigned int&, const double*, const double*, const double*, double*, double*);

  /** \fn void (*computeJacobianVelocityFIntPtr)(const unsigned int& sizeOfq, const double* time, double* qPtr, double* velocityPtr, double* jacobPtr, double* param);
   *  \brief computes the gradient of the the internal strength compared to the velocity
   *  \param unsigned int& sizeOfq : the size of the vector q
   *    \param double* time : the current time
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* jacobPtr : the pointer to the first element of the matrix JacobianQFInt (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeJacobianVelocityFIntPtr)(const unsigned int&, const double*, const double*, const double*, double*, double*);

  /** \fn void (*computeJacobianQNNL)(const unsigned int& sizeOfq, const double* qPtr, double* velocityPtr, double* jacobPtr, double* param)
   *  \brief computes the gradient of the the external strength compared to the state
   *  \param unsigned int& sizeOfq : the size of the vector q
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* jacobPtr : the pointer to the first element of the matrix JacobianQFInt (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeJacobianQNNLPtr)(const unsigned int&, const double*, const double*, double*, double*);

  /** \fn void (*computeJacobianVelocityNNLPtr)(unsigned int* sizeOfq, const double* qPtr, double* velocityPtr, double* jacobPtr, double* param)
   *  \brief computes the gradient of the the external strength compared to the velocity
   *  \param unsigned int& sizeOfq : the size of the vector q
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* jacobPtr : the pointer to the first element of the matrix JacobianQFInt (in-out parameter)
   *    \param double*: vector of parameters
   */
  void (*computeJacobianVelocityNNLPtr)(const unsigned int&, const double*, const double*, double*, double*);


  /** Specific variables to handle links with DynamicalSystem class */
  /** \var workMatrix
   * a map of SiconosMatrix*, zero-matrix, Id-matrix, inverse of Mass or any tmp work matrix
   * No get-set functions at the time. Only used as a protected member.*/
  std::map<std::string, SiconosMatrix*> workMatrix;

  /** \fn void connectToDS()
   *  \brief set links with DS members
   */
  virtual void connectToDS();

  /** \fn initAllocationFlags(const bool& = true);
   *  \brief set all allocation flags (isAllocated map)
   *  \param bool: = if true (default) set default configuration, else set all to false
   */
  virtual void initAllocationFlags(const bool & = true);

  /** \fn initPluginFlags(const bool& val);
   *  \brief set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  virtual void initPluginFlags(const bool&);

  // -- DEFAULT CONSTRUCTOR --
  /** \fn LagrangianDS()
   *  \brief Default constructor
   */
  LagrangianDS();

public:

  // === CONSTRUCTORS - DESTRUCTOR ===

  /** \fn LagrangianDS(DynamicalSystemXML * nsdsXML)
   *  \brief constructor from an xml file
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  LagrangianDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** \fn LagrangianDS(int number, const unsigned int& ndof,
      const SimpleVector& q0, const SimpleVector& velocity0,
      const SiconosMatrix& Mass)
      *  \brief constructor from a minimum set of data
      *  \param int : the number for this DynamicalSystem
      *  \param int : the dimension of this DynamicalSystem
      *  \param SimpleVector : initial coordinates of this DynamicalSystem
      *  \param SimpleVector : initial velocity of this DynamicalSystem
      *  \param SiconosMatrix : mass matrix
      *  \exception RuntimeException
      */
  LagrangianDS(const int&, const unsigned int& ,
               const SimpleVector& , const SimpleVector& ,
               const SiconosMatrix&);

  /** \fn LagrangianDS(int number, const unsigned int& ndof,
      const SimpleVector& q0, const SimpleVector& velocity0,
      const std::string& massPluginName)
      *  \brief constructor from a minimum set of data
      *  \param int : the number for this DynamicalSystem
      *  \param int : the dimension of this DynamicalSystem
      *  \param SimpleVector : initial coordinates of this DynamicalSystem
      *  \param SimpleVector : initial velocity of this DynamicalSystem
      *  \param string: plugin path to compute mass matrix
      *  \exception RuntimeException
      */
  LagrangianDS(const int&, const unsigned int& ,
               const SimpleVector& , const SimpleVector& ,
               const std::string&  = "DefaultPlugin:computeMass");

  /** \fn LagrangianDS(const DynamicalSystem &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  LagrangianDS(const DynamicalSystem &);

  /** \fn ~LagrangianDS()
   *  \brief destructor */
  virtual ~LagrangianDS();

  /** \fn bool checkDynamicalSystem()
   *  \brief check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

  /** \fn void initialize(const double& = 0, const unsigned int& = 1) ;
   *  \brief dynamical system initialization function: mainly set memory and compute plug-in for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const double& = 0, const unsigned int& = 1) ;

  // === GETTERS AND SETTERS ===

  /** \fn const unsigned int getNdof() const
   *  \brief to get the value of ndof
   *  \return the value of ndof
   */
  inline const unsigned int getNdof() const
  {
    return ndof;
  };

  /** \fn void setNdof(const unsigned int&)
   *  \brief to set ndof
   *  \param unsigned int ndof : the value to set ndof
   */
  inline void setNdof(const unsigned int& newNdof)
  {
    ndof = newNdof;
  };

  /** \fn const unsigned int getDim(void) const;
   *  \brief return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  inline const unsigned int getDim(void) const
  {
    return ndof;
  }

  // -- q --

  /** \fn  const SimpleVector getQ() const
   *  \brief get the value of q
   *  \return SimpleVector
   */
  inline const SimpleVector getQ() const
  {
    return *q;
  }

  /** \fn SimpleVector* getQPtr() const
   *  \brief get q
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQPtr() const
  {
    return q;
  }

  /** \fn void setQ (const SimpleVector& newValue)
   *  \brief set the value of q to newValue
   *  \param SimpleVector newValue
   */
  void setQ(const SimpleVector&);

  /** \fn void setQPtr(SimpleVector* newPtr)
   *  \brief set Q to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQPtr(SimpleVector *newPtr);

  // -- q0 --

  /** \fn  const SimpleVector getQ0() const
   *  \brief get the value of q0
   *  \return SimpleVector
   */
  inline const SimpleVector getQ0() const
  {
    return *q0;
  }

  /** \fn SimpleVector* getQ0Ptr() const
   *  \brief get q0
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQ0Ptr() const
  {
    return q0;
  }

  /** \fn void setQ0 (const SimpleVector& newValue)
   *  \brief set the value of q0 to newValue
   *  \param SimpleVector newValue
   */
  void setQ0(const SimpleVector&);

  /** \fn void setQ0Ptr(SimpleVector* newPtr)
   *  \brief set Q0 to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQ0Ptr(SimpleVector *newPtr);

  // -- qFree --

  /** \fn  const SimpleVector getQFree() const
   *  \brief get the value of qFree
   *  \return SimpleVector
   */
  inline const SimpleVector getQFree() const
  {
    return *qFree;
  }

  /** \fn SimpleVector* getQFreePtr() const
   *  \brief get qFree
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQFreePtr() const
  {
    return qFree;
  }

  /** \fn void setQFree (const SimpleVector& newValue)
   *  \brief set the value of qFree to newValue
   *  \param SimpleVector newValue
   */
  void setQFree(const SimpleVector&);

  /** \fn void setQFreePtr(SimpleVector* newPtr)
   *  \brief set QFree to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQFreePtr(SimpleVector *newPtr);

  // Q memory

  /** \fn  const SiconosMemory getQMemory(void) const
   *  \brief get the value of qMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getQMemory() const
  {
    return *qMemory;
  }

  /** \fn SiconosMemory getQMemoryPtr(void) const
   *  \brief get all the values of the state vector q stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getQMemoryPtr() const
  {
    return qMemory;
  }

  /** \fn void setQMemory(const SiconosMemory &)
   *  \brief set the value of qMemory
   *  \param a ref on a SiconosMemory
   */
  void setQMemory(const SiconosMemory&);

  /** \fn void setQMemory(SiconosMemory * newPtr)
   *  \brief set qMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setQMemoryPtr(SiconosMemory *);

  // -- velocity --

  /** \fn  const SimpleVector getVelocity() const
   *  \brief get the value of velocity
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocity() const
  {
    return *velocity;
  }

  /** \fn SimpleVector* getVelocityPtr() const
   *  \brief get velocity
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocityPtr() const
  {
    return velocity;
  }

  /** \fn void setVelocity (const SimpleVector& newValue)
   *  \brief set the value of velocity to newValue
   *  \param SimpleVector newValue
   */
  void setVelocity(const SimpleVector&);

  /** \fn void setVelocityPtr(SimpleVector* newPtr)
   *  \brief set Velocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocityPtr(SimpleVector *newPtr);

  // -- velocity0 --

  /** \fn  const SimpleVector getVelocity0() const
   *  \brief get the value of velocity0
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocity0() const
  {
    return *velocity0;
  }

  /** \fn SimpleVector* getVelocity0Ptr() const
   *  \brief get velocity0
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocity0Ptr() const
  {
    return velocity0;
  }

  /** \fn void setVelocity0 (const SimpleVector& newValue)
   *  \brief set the value of velocity0 to newValue
   *  \param SimpleVector newValue
   */
  void setVelocity0(const SimpleVector&);

  /** \fn void setVelocity0Ptr(SimpleVector* newPtr)
   *  \brief set Velocity0 to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocity0Ptr(SimpleVector *newPtr) ;

  // -- velocityFree --

  /** \fn  const SimpleVector getVelocityFree() const
   *  \brief get the value of velocityFree
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocityFree() const
  {
    return *velocityFree;
  }

  /** \fn SimpleVector* getVelocityFreePtr() const
   *  \brief get velocityFree
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocityFreePtr() const
  {
    return velocityFree;
  }

  /** \fn void setVelocityFree (const SimpleVector& newValue)
   *  \brief set the value of velocityFree to newValue
   *  \param SimpleVector newValue
   */
  void setVelocityFree(const SimpleVector&);

  /** \fn void setVelocityFreePtr(SimpleVector* newPtr)
   *  \brief set VelocityFree to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocityFreePtr(SimpleVector *newPtr);

  // Velocity memory

  /** \fn  const SiconosMemory getVelocityMemory(void) const
   *  \brief get the value of velocityMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getVelocityMemory() const
  {
    return *velocityMemory;
  }

  /** \fn SiconosMemory getVelocityMemoryPtr(void) const
   *  \brief get all the values of the state vector velocity stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getVelocityMemoryPtr() const
  {
    return velocityMemory;
  }

  /** \fn void setVelocityMemory(const SiconosMemory &)
   *  \brief set the value of velocityMemory
   *  \param a ref on a SiconosMemory
   */
  void setVelocityMemory(const SiconosMemory&);

  /** \fn void setVelocityMemory(SiconosMemory * newPtr)
   *  \brief set velocityMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setVelocityMemoryPtr(SiconosMemory *);

  // -- p --

  /** \fn  const SimpleVector getP() const
   *  \brief get the value of p
   *  \return SimpleVector
   */
  inline const SimpleVector getP() const
  {
    return *p;
  }

  /** \fn SimpleVector* getPPtr() const
   *  \brief get p
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getPPtr() const
  {
    return p;
  }

  /** \fn void setP (const SimpleVector& newValue)
   *  \brief set the value of p to newValue
   *  \param SimpleVector newValue
   */
  void setP(const SimpleVector&);

  /** \fn void setPPtr(SimpleVector* newPtr)
   *  \brief set P to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setPPtr(SimpleVector *newPtr);

  // -- Mass --

  /** \fn  const SimpleMatrix getMass() const
   *  \brief get the value of Mass
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getMass() const
  {
    return *mass;
  }

  /** \fn SiconosMatrix* getMassPtr() const
   *  \brief get Mass
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getMassPtr() const
  {
    return mass;
  }

  /** \fn void setMass (const SiconosMatrix& newValue)
   *  \brief set the value of Mass to newValue
   *  \param SiconosMatrix newValue
   */
  void setMass(const SiconosMatrix&);

  /** \fn void setMassPtr(SiconosMatrix* newPtr)
   *  \brief set Mass to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setMassPtr(SiconosMatrix *newPtr);

  // -- FInt --

  /** \fn  const SimpleVector getFInt() const
   *  \brief get the value of fInt
   *  \return SimpleVector
   */
  inline const SimpleVector getFInt() const
  {
    return *fInt;
  }

  /** \fn SimpleVector* getFIntPtr() const
   *  \brief get fInt
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getFIntPtr() const
  {
    return fInt;
  }

  /** \fn void setFInt (const SimpleVector& newValue)
   *  \brief set the value of fInt to newValue
   *  \param SimpleVector newValue
   */
  void setFInt(const SimpleVector&);

  /** \fn void setFIntPtr(SimpleVector* newPtr)
   *  \brief set FInt to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setFIntPtr(SimpleVector *newPtr);

  // -- Fext --

  /** \fn  const SimpleVector getFExt() const
   *  \brief get the value of fExt
   *  \return SimpleVector
   */
  inline const SimpleVector getFExt() const
  {
    return *fExt;
  }

  /** \fn SimpleVector* getFExtPtr() const
   *  \brief get fExt
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getFExtPtr() const
  {
    return fExt;
  }

  /** \fn void setFExt (const SimpleVector& newValue)
   *  \brief set the value of fExt to newValue
   *  \param SimpleVector newValue
   */
  void setFExt(const SimpleVector&);

  /** \fn void setFExtPtr(SimpleVector* newPtr)
   *  \brief set FExt to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setFExtPtr(SimpleVector *newPtr);

  // -- NNL --

  /** \fn  const SimpleVector getNNL() const
   *  \brief get the value of NNL
   *  \return SimpleVector
   */
  inline const SimpleVector getNNL() const
  {
    return *NNL;
  }

  /** \fn SimpleVector* getNNLPtr() const
   *  \brief get NNL
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getNNLPtr() const
  {
    return NNL;
  }

  /** \fn void setNNL (const SimpleVector& newValue)
   *  \brief set the value of NNL to newValue
   *  \param SimpleVector newValue
   */
  void setNNL(const SimpleVector&);

  /** \fn void setNNLPtr(SimpleVector* newPtr)
   *  \brief set NNL to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setNNLPtr(SimpleVector *newPtr);

  // -- Jacobian Q Fint --

  /** \fn  const SimpleMatrix getJacobianQFInt() const
   *  \brief get the value of JacobianQFInt
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianQFInt() const
  {
    return *jacobianQFInt;
  }

  /** \fn SiconosMatrix* getJacobianQFIntPtr() const
   *  \brief get JacobianQFInt
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianQFIntPtr() const
  {
    return jacobianQFInt;
  }

  /** \fn void setJacobianQFInt (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianQFInt to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianQFInt(const SiconosMatrix&);

  /** \fn void setJacobianQFIntPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianQFInt to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianQFIntPtr(SiconosMatrix *newPtr);

  // -- Jacobian velocity Fint --

  /** \fn  const SimpleMatrix getJacobianVelocityFInt() const
   *  \brief get the value of JacobianVelocityFInt
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianVelocityFInt() const
  {
    return *jacobianVelocityFInt;
  }

  /** \fn SiconosMatrix* getJacobianVelocityFIntPtr() const
   *  \brief get JacobianVelocityFInt
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianVelocityFIntPtr() const
  {
    return jacobianVelocityFInt;
  }

  /** \fn void setJacobianVelocityFInt (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianVelocityFInt to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianVelocityFInt(const SiconosMatrix&);

  /** \fn void setJacobianVelocityFIntPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianVelocityFInt to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianVelocityFIntPtr(SiconosMatrix *newPtr);

  // -- Jacobian Q NNL --

  /** \fn  const SimpleMatrix getJacobianQNNL() const
   *  \brief get the value of JacobianQNNL
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianQNNL() const
  {
    return *jacobianQNNL;
  }

  /** \fn SiconosMatrix* getJacobianQNNLPtr() const
   *  \brief get JacobianQNNL
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianQNNLPtr() const
  {
    return jacobianQNNL;
  }

  /** \fn void setJacobianQNNL (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianQNNL to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianQNNL(const SiconosMatrix&);

  /** \fn void setJacobianQNNLPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianQNNL to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianQNNLPtr(SiconosMatrix *newPtr);

  // -- Jacobian velocity NNL --

  /** \fn  const SimpleMatrix getJacobianVelocityNNL() const
   *  \brief get the value of JacobianVelocityNNL
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianVelocityNNL() const
  {
    return *jacobianVelocityNNL;
  }

  /** \fn SiconosMatrix* getJacobianVelocityNNLPtr() const
   *  \brief get JacobianVelocityNNL
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianVelocityNNLPtr() const
  {
    return jacobianVelocityNNL;
  }

  /** \fn void setJacobianVelocityNNL (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianVelocityNNL to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianVelocityNNL(const SiconosMatrix&);

  /** \fn void setJacobianVelocityNNLPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianVelocityNNL to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianVelocityNNLPtr(SiconosMatrix *newPtr);

  /** \fn  std::string getMassFunctionName() const
   *  \brief get name of function that computes mass (if mass from plugin)
   *  \return a string
   */
  inline const std::string getMassFunctionName() const
  {
    return massFunctionName;
  }

  /** \fn  std::string getFIntFunctionName() const
   *  \brief get name of function that computes fInt (if fInt from plugin)
   *  \return a string
   */
  inline const std::string getFIntFunctionName() const
  {
    return fIntFunctionName;
  }

  /** \fn  std::string getFExtFunctionName() const
   *  \brief get name of function that computes fExt (if fExt from plugin)
   *  \return a string
   */
  inline const std::string getFExtFunctionName() const
  {
    return fExtFunctionName;
  }

  /** \fn  std::string getNNLFunctionName() const
   *  \brief get name of function that computes NNL (if NNL from plugin)
   *  \return a string
   */
  inline const std::string getNNLFunctionName() const
  {
    return NNLFunctionName;
  }

  /** \fn  std::string getJacobianQFIntFunctionName() const
   *  \brief get name of function that computes jacobianQFInt (if jacobianQFInt from plugin)
   *  \return a string
   */
  inline const std::string getJacobianQFIntFunctionName() const
  {
    return jacobianQFIntFunctionName;
  }

  /** \fn  std::string getJacobianVelocityFIntFunctionName() const
   *  \brief get name of function that computes jacobianVelocityFInt (if jacobianVelocityFInt from plugin)
   *  \return a string
   */
  inline const std::string getJacobianVelocityFIntFunctionName() const
  {
    return jacobianVelocityFIntFunctionName;
  }

  /** \fn  std::string getJacobianQNNLFunctionName() const
   *  \brief get name of function that computes jacobianQNNL (if jacobianQNNL from plugin)
   *  \return a string
   */
  inline const std::string getJacobianQNNLFunctionName() const
  {
    return jacobianQNNLFunctionName;
  }

  /** \fn  std::string getJacobianVelocityNNLFunctionName() const
   *  \brief get name of function that computes jacobianVelocityNNL (if jacobianVelocityNNL from plugin)
   *  \return a string
   */
  inline const std::string getJacobianVelocityNNLFunctionName() const
  {
    return jacobianVelocityNNLFunctionName;
  }

  // --- PLUGINS RELATED FUNCTIONS ---

  /** \fn void setComputeMassFunction(const string pluginPath, const string functionName&)
   *  \brief allow to set a specified function to compute the mass
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeMassFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeFIntFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeFExtFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception ICDLL_CSharedLibraryException
   */
  void setComputeFExtFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeNNLFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the inertia
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosCSharedLibraryException
   */
  void setComputeNNLFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeJacobianQFIntFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the gradient of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQFIntFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeJacobianVelocityFIntFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the internal strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityFIntFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeJacobianQNNLFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the gradient of the the external strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQNNLFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeJacobianVelocityNNLFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the external strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityNNLFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void computeMass()
   *  \brief default function to compute the mass
   *  \exception RuntimeException
   */
  void computeMass();

  /** \fn void computeMass(SimpleVector *q)
   *  \brief function to compute the mass
   *  \param double time : the current time, SimpleVector*: pointer on the state vector q
   *  \exception RuntimeException
   */
  void computeMass(SimpleVector *);

  /** \fn void computeFInt(const double &)
   *  \brief default function to compute the internal strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeFInt(const double &);

  /** \fn void computeFInt(const double &, SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the internal strengths
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  void computeFInt(const double &, SimpleVector *, SimpleVector *);

  /** \fn void computeFExt(const double &)
   *  \brief default function to compute the external strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeFExt(const double &);

  /** \fn void computeNNL();
   *  \brief default function to compute the inertia
   *  \exception RuntimeException
   */
  void computeNNL();

  /** \fn void computeNNL(SimpleVector q, SimpleVector velocity);
   *  \brief function to compute the inertia
   *  \param SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  void computeNNL(SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianQFInt(const double &)
   *  \brief default function to compute the gradient of the internal strengths compared to the state
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeJacobianQFInt(const double &);

  /** \fn void computeJacobianQFInt(const double &,SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the gradient of the internal strengths compared to state q
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  void computeJacobianQFInt(const double &, SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianVelocityFInt(double time)
   *  \brief function to compute the gradient of the internal strengths compared to velocity
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeJacobianVelocityFInt(const double &);

  /** \fn void computeJacobianVelocityFInt(const double &, SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the gradient of the internal strengths compared to velocity
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  void computeJacobianVelocityFInt(const double &, SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianQNNL(const double &)
   *  \brief function to compute the gradient of the inertia strengths compared to the state q
   *  \exception RuntimeException
   */
  void computeJacobianQNNL();

  /** \fn void computeJacobianQNNL(const double &,SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the gradient of the inertia strengths compared to the state q
   *  \param SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  void computeJacobianQNNL(SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianVelocityNNL(const double & )
   *  \brief function to compute the gradient of the inertia strengths compared to velocity
   */
  void computeJacobianVelocityNNL();

  /** \fn void computeJacobianVelocityNNL(double time, SimpleVector q, SimpleVector velocity )
   *  \brief function to compute the gradient of the inertia strengths compared to velocity
   *  \param SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   */
  void computeJacobianVelocityNNL(SimpleVector *q, SimpleVector *velocity);

  /** \fn void  computeInverseOfMass();
   *  \brief function to compute inverse of the mass matrix
   */
  void computeInverseOfMass();

  /** \fn void computeRhs(const double& time, const bool & =false)
   *  \brief Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeRhs(const double&, const bool & = false);

  /** \fn void computeJacobianXRhs(const double& time, const bool & =false)
   *  \brief Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXRhs(const double&, const bool & = false);

  // --- miscellaneous ---

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  virtual void display() const;

  /** \fn void initMemory(const unsigned int& steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  void initMemory(const unsigned int& steps);

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory();

  /** \fn LagrangianDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the system if it is of the right type, NULL otherwise
   */
  static LagrangianDS* convert(DynamicalSystem* ds);

  /** \fn double dsConvergenceIndicator()
   *  \brief compute $\frac{|\dot q_{i+1} - \dot qi|}{|\dot q_i|}$ where $\dot q_{i+1}$ represents the present state and $\dot q_i$ the previous one
   * \return a double
   */
  virtual double dsConvergenceIndicator();

};

#endif // LAGRANGIANNLDS_H

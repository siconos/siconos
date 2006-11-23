/* Siconos-Kernel version 1.3.0, Copyright INRIA 2005-2006.
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

/*! \file LagrangianDS.h
  \brief LagrangianDS class

*/

#ifndef LAGRANGIANNLDS_H
#define LAGRANGIANNLDS_H
#include "DynamicalSystem.h"
#include "LagrangianDSXML.h"
#include "BlockMatrix.h"
#include<map>

class LagrangianDSXML;

//! Lagrangian non linear dynamical systems - Derived from DynamicalSystem -
/**  \author SICONOS Development Team - copyright INRIA
 *  \version 1.3.0.
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
 * \nabla_{q}M^{-1}(q)\left[F_{Ext}( q , t) - F_{Int}(\dot q , q , t) - NNL(\dot q, q) \right]-M^{-1}(q)\left[\nabla_q(F_{Int}(\dot q , q , t)+NNL(\dot q, q))\right] & -M^{-1}(q)\left[\nabla_{\dot q}(F_{Int}(\dot q , q , t)+NNL(\dot q, q))\right] \\
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
 *    - Initialization: compute state members and operators for time=t0 (usually done when calling simulation->initialize)
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

  /** coordinates of the system */
  SimpleVector *q;
  /** initial coordinates of the system */
  SimpleVector *q0;
  /** free state of the system */
  SimpleVector *qFree;
  /** memory of previous coordinates of the system */
  SiconosMemory *qMemory;

  SimpleVector *velocity;/**<velocity of the system*/
  SimpleVector *velocity0; /**< initial velocity of the system*/
  SimpleVector *velocityFree;/**<free Velocity*/
  SiconosMemory *velocityMemory;/**<memory of previous velocities of the system */

  /** "Reaction" due to the non smooth law - The index corresponds to the dynamic level. */
  std::vector<SimpleVector*> p;

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

  /** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param[in,out] mass : pointer to the first element of mass
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*computeMassPtr)(unsigned int, const double*, double*, double*);

  /** LagrangianDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
   * @param sizeOfq : size of vector q
   * @param time : current time
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] fInt : pointer to the first element of fInt
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*computeFIntPtr)(unsigned int, double, const double*, const double*, double*, double*);

  /** LagrangianDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
   * @param sizeOfq : size of vector q
   * @param time : current time
   * @param[in,out] fExt : pointer to the first element of fExt
   * @param[in,out] param : a vector of user-defined parameters
   */
  void (*computeFExtPtr)(const unsigned int, double, double*, double*);

  /** LagrangianDS plug-in to compute \f$NNL(\dot q, q)\f$, id = "NNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] NNL : pointer to the first element of NNL
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*computeNNLPtr)(unsigned int, const double*, const double*, double*, double*);

  /** LagrangianDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianQFInt"
   * @param sizeOfq : size of vector q
   * @param time : current time
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*computeJacobianQFIntPtr)(unsigned int, double, const double*, const double*, double*, double*);

  /** LagrangianDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianVelocityFInt"
   * @param sizeOfq : size of vector q
   * @param time : current time
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*computeJacobianVelocityFIntPtr)(unsigned int, double, const double*, const double*, double*, double*);

  /** LagrangianDS plug-in to compute \f$\nabla_qNNL(\dot q, q)\f$, id = "jacobianQNNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*computeJacobianQNNLPtr)(const unsigned int, const double*, const double*, double*, double*);

  /** LagrangianDS plug-in to compute \f$\nabla_{\dot q}NNL(\dot q, q)\f$, id = "jacobianVelocityNNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param[in,out] param  : a vector of user-defined parameters
   */
  void (*computeJacobianVelocityNNLPtr)(unsigned int, const double*, const double*, double*, double*);

  /** Specific variables to handle links with DynamicalSystem class */
  /** \var workMatrix
   * a map of SiconosMatrix*, zero-matrix, Id-matrix, inverse of Mass or any tmp work matrix
   * No get-set functions at the time. Only used as a protected member.*/
  std::map<std::string, SiconosMatrix*> workMatrix;

  /** set links with DS members
   */
  virtual void connectToDS();

  /** set all allocation flags (isAllocated map)
   *  \param bool: = if true (default) set default configuration, else set all to false
   */
  virtual void initAllocationFlags(const bool  = true);

  /** set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  virtual void initPluginFlags(const bool);

  // -- DEFAULT CONSTRUCTOR --
  /** Default constructor
   */
  LagrangianDS();

public:

  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from an xml file
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  LagrangianDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** constructor from a minimum set of data
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SimpleVector : initial coordinates of this DynamicalSystem
   *  \param SimpleVector : initial velocity of this DynamicalSystem
   *  \param SiconosMatrix : mass matrix
   *  \exception RuntimeException
   */
  LagrangianDS(int, unsigned int, const SimpleVector&, const SimpleVector&, const SiconosMatrix&);

  /** constructor from a minimum set of data
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SimpleVector : initial coordinates of this DynamicalSystem
   *  \param SimpleVector : initial velocity of this DynamicalSystem
   *  \param string: plugin path to compute mass matrix
   *  \exception RuntimeException
   */
  LagrangianDS(int, unsigned int, const SimpleVector& , const SimpleVector&, const std::string& = "DefaultPlugin:computeMass");

  /** copy constructor
   *  \param a Dynamical system to copy
   */
  LagrangianDS(const DynamicalSystem &);

  /** destructor */
  virtual ~LagrangianDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

  /** initialization of qFree, vFree
   *  \param a string: the simulation type. For TimeStepping: memory allocation. For EventDriven: links (pointers) to q and velocity.
   */
  void initFreeVectors(const std::string);

  /** allocate memory for p[...] vectors
   *  \param string: simulation type
   */
  void initP(const std::string);

  /** dynamical system initialization function: mainly set memory and compute plug-in for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const std::string&, double = 0, unsigned int = 1) ;

  /** dynamical system update: mainly call compute for all time or state depending functions (mass, FInt ...).
   *  \param current time
   */
  virtual void update(const double);

  // === GETTERS AND SETTERS ===

  /** to get the value of ndof
   *  \return the value of ndof
   */
  inline const unsigned int getNdof() const
  {
    return ndof;
  };

  /** to set ndof
   *  \param unsigned int ndof : the value to set ndof
   */
  inline void setNdof(const unsigned int newNdof)
  {
    ndof = newNdof;
  };

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  inline const unsigned int getDim(void) const
  {
    return ndof;
  }

  // -- q --

  /** get the value of q
   *  \return SimpleVector
   */
  inline const SimpleVector getQ() const
  {
    return *q;
  }

  /** get q
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQPtr() const
  {
    return q;
  }

  /** set the value of q to newValue
   *  \param SimpleVector newValue
   */
  void setQ(const SimpleVector&);

  /** set Q to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQPtr(SimpleVector *newPtr);

  // -- q0 --

  /** get the value of q0
   *  \return SimpleVector
   */
  inline const SimpleVector getQ0() const
  {
    return *q0;
  }

  /** get q0
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQ0Ptr() const
  {
    return q0;
  }

  /** set the value of q0 to newValue
   *  \param SimpleVector newValue
   */
  void setQ0(const SimpleVector&);

  /** set Q0 to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQ0Ptr(SimpleVector *newPtr);

  // -- qFree --

  /** get the value of qFree
   *  \return SimpleVector
   */
  inline const SimpleVector getQFree() const
  {
    return *qFree;
  }

  /** get qFree
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQFreePtr() const
  {
    return qFree;
  }

  /** set the value of qFree to newValue
   *  \param SimpleVector newValue
   */
  void setQFree(const SimpleVector&);

  /** set QFree to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQFreePtr(SimpleVector *newPtr);

  // Q memory

  /** get the value of qMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getQMemory() const
  {
    return *qMemory;
  }

  /** get all the values of the state vector q stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getQMemoryPtr() const
  {
    return qMemory;
  }

  /** set the value of qMemory
   *  \param a ref on a SiconosMemory
   */
  void setQMemory(const SiconosMemory&);

  /** set qMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setQMemoryPtr(SiconosMemory *);

  // -- velocity --

  /** get the value of velocity
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocity() const
  {
    return *velocity;
  }

  /** get velocity
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocityPtr() const
  {
    return velocity;
  }

  /** set the value of velocity to newValue
   *  \param SimpleVector newValue
   */
  void setVelocity(const SimpleVector&);

  /** set Velocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocityPtr(SimpleVector *newPtr);

  // -- velocity0 --

  /** get the value of velocity0
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocity0() const
  {
    return *velocity0;
  }

  /** get velocity0
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocity0Ptr() const
  {
    return velocity0;
  }

  /** set the value of velocity0 to newValue
   *  \param SimpleVector newValue
   */
  void setVelocity0(const SimpleVector&);

  /** set Velocity0 to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocity0Ptr(SimpleVector *newPtr) ;

  // -- velocityFree --

  /** get the value of velocityFree
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocityFree() const
  {
    return *velocityFree;
  }

  /** get velocityFree
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocityFreePtr() const
  {
    return velocityFree;
  }

  /** set the value of velocityFree to newValue
   *  \param SimpleVector newValue
   */
  void setVelocityFree(const SimpleVector&);

  /** set VelocityFree to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocityFreePtr(SimpleVector *newPtr);

  // -- acceleration --

  /** get acceleration
   *  \return pointer on a SimpleVector
   */
  SimpleVector* getAccelerationPtr() const ;

  // Velocity memory

  /** get the value of velocityMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getVelocityMemory() const
  {
    return *velocityMemory;
  }

  /** get all the values of the state vector velocity stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getVelocityMemoryPtr() const
  {
    return velocityMemory;
  }

  /** set the value of velocityMemory
   *  \param a ref on a SiconosMemory
   */
  void setVelocityMemory(const SiconosMemory&);

  /** set velocityMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setVelocityMemoryPtr(SiconosMemory *);

  // -- p --

  /** get the value of p[index]
   *  \param unsigned int, required level for p, default = 2
   *  \return SimpleVector
   */
  inline const SimpleVector getP(const unsigned int level = 2) const
  {
    return *(p[level]);
  }

  /** get p
   *  \param unsigned int, required level for p, default = 2
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getPPtr(const unsigned int level = 2) const
  {
    return p[level];
  }

  /** set the value of p to newValue
   *  \param unsigned int, required level for p, default = 2
   *  \param SimpleVector newValue
   */
  void setP(const SimpleVector&, const unsigned int level = 2);

  /** set P to pointer newPtr
   *  \param unsigned int, required level for p, default = 2
   *  \param SimpleVector * newPtr
   */
  void setPPtr(SimpleVector *newPtr, const unsigned int level = 2);

  // -- Mass --

  /** get the value of Mass
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getMass() const
  {
    return *mass;
  }

  /** get Mass
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getMassPtr() const
  {
    return mass;
  }

  /** get inverse of Mass - Warning: in this function we do not checked that tha matrix is up to date. If M depends on q, it may require a recomputation before the get.
   *  \return pointer to a SiconosMatrix
   */
  SiconosMatrix* getInverseOfMassPtr();

  /** set the value of Mass to newValue
   *  \param SiconosMatrix newValue
   */
  void setMass(const SiconosMatrix&);

  /** set Mass to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setMassPtr(SiconosMatrix *newPtr);

  // -- FInt --

  /** get the value of fInt
   *  \return SimpleVector
   */
  inline const SimpleVector getFInt() const
  {
    return *fInt;
  }

  /** get fInt
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getFIntPtr() const
  {
    return fInt;
  }

  /** set the value of fInt to newValue
   *  \param SimpleVector newValue
   */
  void setFInt(const SimpleVector&);

  /** set FInt to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setFIntPtr(SimpleVector *newPtr);

  // -- Fext --

  /** get the value of fExt
   *  \return SimpleVector
   */
  inline const SimpleVector getFExt() const
  {
    return *fExt;
  }

  /** get fExt
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getFExtPtr() const
  {
    return fExt;
  }

  /** set the value of fExt to newValue
   *  \param SimpleVector newValue
   */
  void setFExt(const SimpleVector&);

  /** set FExt to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setFExtPtr(SimpleVector *newPtr);

  // -- NNL --

  /** get the value of NNL
   *  \return SimpleVector
   */
  inline const SimpleVector getNNL() const
  {
    return *NNL;
  }

  /** get NNL
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getNNLPtr() const
  {
    return NNL;
  }

  /** set the value of NNL to newValue
   *  \param SimpleVector newValue
   */
  void setNNL(const SimpleVector&);

  /** set NNL to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setNNLPtr(SimpleVector *newPtr);

  // -- Jacobian Q Fint --

  /** get the value of JacobianQFInt
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianQFInt() const
  {
    return *jacobianQFInt;
  }

  /** get JacobianQFInt
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianQFIntPtr() const
  {
    return jacobianQFInt;
  }

  /** set the value of JacobianQFInt to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianQFInt(const SiconosMatrix&);

  /** set JacobianQFInt to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianQFIntPtr(SiconosMatrix *newPtr);

  // -- Jacobian velocity Fint --

  /** get the value of JacobianVelocityFInt
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianVelocityFInt() const
  {
    return *jacobianVelocityFInt;
  }

  /** get JacobianVelocityFInt
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianVelocityFIntPtr() const
  {
    return jacobianVelocityFInt;
  }

  /** set the value of JacobianVelocityFInt to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianVelocityFInt(const SiconosMatrix&);

  /** set JacobianVelocityFInt to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianVelocityFIntPtr(SiconosMatrix *newPtr);

  // -- Jacobian Q NNL --

  /** get the value of JacobianQNNL
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianQNNL() const
  {
    return *jacobianQNNL;
  }

  /** get JacobianQNNL
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianQNNLPtr() const
  {
    return jacobianQNNL;
  }

  /** set the value of JacobianQNNL to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianQNNL(const SiconosMatrix&);

  /** set JacobianQNNL to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianQNNLPtr(SiconosMatrix *newPtr);

  // -- Jacobian velocity NNL --

  /** get the value of JacobianVelocityNNL
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianVelocityNNL() const
  {
    return *jacobianVelocityNNL;
  }

  /** get JacobianVelocityNNL
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianVelocityNNLPtr() const
  {
    return jacobianVelocityNNL;
  }

  /** set the value of JacobianVelocityNNL to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianVelocityNNL(const SiconosMatrix&);

  /** set JacobianVelocityNNL to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianVelocityNNLPtr(SiconosMatrix *newPtr);

  /** get name of function that computes mass (if mass from plugin)
   *  \return a string
   */
  inline const std::string getMassFunctionName() const
  {
    return massFunctionName;
  }

  /** get name of function that computes fInt (if fInt from plugin)
   *  \return a string
   */
  inline const std::string getFIntFunctionName() const
  {
    return fIntFunctionName;
  }

  /** get name of function that computes fExt (if fExt from plugin)
   *  \return a string
   */
  inline const std::string getFExtFunctionName() const
  {
    return fExtFunctionName;
  }

  /** get name of function that computes NNL (if NNL from plugin)
   *  \return a string
   */
  inline const std::string getNNLFunctionName() const
  {
    return NNLFunctionName;
  }

  /** get name of function that computes jacobianQFInt (if jacobianQFInt from plugin)
   *  \return a string
   */
  inline const std::string getJacobianQFIntFunctionName() const
  {
    return jacobianQFIntFunctionName;
  }

  /** get name of function that computes jacobianVelocityFInt (if jacobianVelocityFInt from plugin)
   *  \return a string
   */
  inline const std::string getJacobianVelocityFIntFunctionName() const
  {
    return jacobianVelocityFIntFunctionName;
  }

  /** get name of function that computes jacobianQNNL (if jacobianQNNL from plugin)
   *  \return a string
   */
  inline const std::string getJacobianQNNLFunctionName() const
  {
    return jacobianQNNLFunctionName;
  }

  /** get name of function that computes jacobianVelocityNNL (if jacobianVelocityNNL from plugin)
   *  \return a string
   */
  inline const std::string getJacobianVelocityNNLFunctionName() const
  {
    return jacobianVelocityNNLFunctionName;
  }

  // --- PLUGINS RELATED FUNCTIONS ---

  /** allow to set a specified function to compute the mass
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeMassFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception ICDLL_CSharedLibraryException
   */
  void setComputeFExtFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute the inertia
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosCSharedLibraryException
   */
  void setComputeNNLFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute the gradient of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQFIntFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute the internal strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityFIntFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute the gradient of the the external strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQNNLFunction(const std::string  pluginPath, const std::string  functionName);

  /** allow to set a specified function to compute the external strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityNNLFunction(const std::string  pluginPath, const std::string  functionName);

  /** default function to compute the mass
   *  \exception RuntimeException
   */
  void computeMass();

  /** function to compute the mass
   *  \param double time : the current time, SimpleVector*: pointer on the state vector q
   *  \exception RuntimeException
   */
  void computeMass(SimpleVector *);

  /** default function to compute the internal strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeFInt(const double);

  /** function to compute the internal strengths
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity
   *  \exception RuntimeException
   */
  void computeFInt(const double , SimpleVector *, SimpleVector *);

  /** default function to compute the external strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeFExt(const double);

  /** default function to compute the inertia
   *  \exception RuntimeException
   */
  void computeNNL();

  /** function to compute the inertia
   *  \param SimpleVector*: pointers on the state vectors q and velocity
   *  \exception RuntimeException
   */
  void computeNNL(SimpleVector *q, SimpleVector *velocity);

  /** default function to compute the gradient of the internal strengths compared to the state
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeJacobianQFInt(const double);

  /** function to compute the gradient of the internal strengths compared to state q
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity
   *  \exception RuntimeException
   */
  void computeJacobianQFInt(const double , SimpleVector *q, SimpleVector *velocity);

  /** function to compute the gradient of the internal strengths compared to velocity
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  void computeJacobianVelocityFInt(const double);

  /** function to compute the gradient of the internal strengths compared to velocity
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity
   *  \exception RuntimeException
   */
  void computeJacobianVelocityFInt(const double , SimpleVector *q, SimpleVector *velocity);

  /** function to compute the gradient of the inertia strengths compared to the state q
   *  \exception RuntimeException
   */
  void computeJacobianQNNL();

  /** function to compute the gradient of the inertia strengths compared to the state q
   *  \param SimpleVector*: pointers on the state vectors q and velocity
   *  \exception RuntimeException
   */
  void computeJacobianQNNL(SimpleVector *q, SimpleVector *velocity);

  /** function to compute the gradient of the inertia strengths compared to velocity
   */
  void computeJacobianVelocityNNL();

  /** function to compute the gradient of the inertia strengths compared to velocity
   *  \param SimpleVector*: pointers on the state vectors q and velocity
   */
  void computeJacobianVelocityNNL(SimpleVector *q, SimpleVector *velocity);

  /** function to compute inverse of the mass matrix
   */
  void computeInverseOfMass();

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

  // --- miscellaneous ---

  /** copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** print the data to the screen
   */
  virtual void display() const;

  /** initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  void initMemory(const unsigned int steps);

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the system if it is of the right type, NULL otherwise
   */
  static LagrangianDS* convert(DynamicalSystem* ds);

  /** compute \f$\frac{|\dot q_{i+1} - \dot qi|}{|\dot q_i|}\f$ where \f$\dot q_{i+1}\f$ represents the present state and \f$\dot q_i\f$ the previous one
   * \return a double
   */
  virtual double dsConvergenceIndicator();

  /** function to compute derivative number level of qFree
   *  \param double: current time
   *  \param unsigned int: derivative number
   *  \param SimpleVector*: in-out parameter, qFree
   */
  void computeQFree(const double, const unsigned int, SiconosVector*);

  /** set p[...] to zero
   */
  void resetNonSmoothPart();
};

#endif // LAGRANGIANNLDS_H

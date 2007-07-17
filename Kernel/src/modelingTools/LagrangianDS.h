/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2006.
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
  \brief LagrangianDS class - Second Order Non Linear Dynamical Systems.

*/

#ifndef LAGRANGIANNLDS_H
#define LAGRANGIANNLDS_H

#include "DynamicalSystem.h"

/** Pointer to function for plug-in. For NNL and its jacobian. */
typedef void (*FPtr5)(unsigned int, const double*, const double*, double*, unsigned int, double*);

class DynamicalSystem;

/** Lagrangian non linear dynamical systems - Derived from DynamicalSystem -
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianDS  defines  and computes a generic ndof-dimensional
 * Lagrangian Non Linear Dynamical System of the form :
 * \f[
 * M(q,z) \ddot q + NNL(\dot q, q, z) + F_{Int}(\dot q , q , t, z) = F_{Ext}(t, z) + p
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
 *    -  \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *
 *  Or:
  * \f[
 * M(q,z) \ddot q = f_L(\dot q, q, t, z) + p
 * \f]
 *
 * Links with first order DynamicalSystem top-class are:
 *
 * \f$ n= 2 ndof \f$
 * \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
 *
 * The rhs is given by:
 * \f[
 * \dot x = \left[\begin{array}{c}
 *  \dot q  \\
 * \ddot q = M^{-1}(q)\left[f_L(\dot q, q , t, z) + p \right]\\
 * \end{array}\right]
 * \f]
 * Its jacobian is:
 * \f[
 * \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
 *  0  & I \\
 * \nabla_{q}(M^{-1}(q)f_L(\dot q, q , t, z)) &  \nabla_{\dot q}(M^{-1}(q)f_L(\dot q, q , t, z)) \\
 * \end{array}\right]
 * \f]
 *  The input due to the non smooth law is:
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
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
 *        => computeRhs            (no set function)
 *        => computeJacobianXRhs   (no set function)
 *
 * About notation:
 *    - q[i] is the derivative number i of q.
 * Thus: q[0]=\f$ q \f$, global coordinates, q[1]=\f$ \dot q\f$, velocity, q[2]=\f$ \ddot q \f$, acceleration.
 *    - qFree is the state of the system when non-smooth effects are not taken into account. If required they are saved in workVector.
 *
 *
 */
class LagrangianDS : public DynamicalSystem
{
protected:

  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  unsigned int ndof;

  /** state of the system. See details on top of page. */
  VectorOfVectors q;

  /** initial coordinates of the system */
  SiconosVector* q0;

  /** initial velocity of the system */
  SiconosVector* velocity0;

  /** memory of previous coordinates of the system */
  SiconosMemory *qMemory;

  /** memory of previous velocities of the system */
  SiconosMemory *velocityMemory;

  /** "Reaction" due to the non smooth law - The index corresponds to the dynamic level. */
  std::vector<SiconosVector*> p;

  /** mass of the system */
  SiconosMatrix *mass;

  /** internal strength of the system */
  SiconosVector *fInt;

  /** jacobian/coordinates, jacobianFInt[0] and velocity, jacobianFInt[1], of internal strength */
  VectorOfMatrices jacobianFInt;

  /** external strength of the system */
  SiconosVector *fExt;

  /** non-linear inertia term of the system */
  SiconosVector *NNL;

  /** jacobian/coordinates, jacobianNNL[0] and velocity, jacobianNNL[1], of NNL */
  VectorOfMatrices jacobianNNL;

  /** fL(q[0],q[1],t)= fExt - fInt -NNL */
  SiconosVector * fL;

  /** jacobian/coordinates, jacobianFL[0], and velocity, jacobianFL[1], of fL */
  VectorOfMatrices jacobianFL;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  // pointers to functions member to compute plug-in functions

  /** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param[in,out] mass : pointer to the first element of mass
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  void (*computeMassPtr)(unsigned int, const double*, double*, unsigned int, double*);

  /** LagrangianDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] fInt : pointer to the first element of fInt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr6 computeFIntPtr;

  /** LagrangianDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param[in,out] fExt : pointer to the first element of fExt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  void (*computeFExtPtr)(double, unsigned int, double*, unsigned int, double*);

  /** LagrangianDS plug-in to compute \f$NNL(\dot q, q)\f$, id = "NNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] NNL : pointer to the first element of NNL
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr5 computeNNLPtr;

  /** LagrangianDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianFInt0" and \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianFInt1"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  std::vector<FPtr6> computeJacobianFIntPtr;

  /** LagrangianDS plug-in to compute \f$\nabla_qNNL(\dot q, q)\f$, id = "jacobianNNL0" and \f$\nabla_{\dot q}NNL(\dot q, q)\f$, id = "jacobianNNL1"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  std::vector<FPtr5> computeJacobianNNLPtr;

  /** Map that links operators names with a bool. If bool is true, the operator
   is up to date, else not. Uselful to avoid re-computation of mass, fInt ... */
  BoolMap isUp;

  /** set links with DS members
   */
  void connectToDS();

  /** set all allocation flags (isAllocated map)
   *  \param bool: = if true (default) set default configuration, else set all to false
   */
  virtual void initAllocationFlags(bool  = true);

  /** set all plug-in flags (isPlugin map) to val
   *  \param a bool
   */
  virtual void initPluginFlags(bool);

  /** Default constructor
   */
  LagrangianDS();

public:

  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from an xml file
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   */
  LagrangianDS(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** constructor from a minimum set of data
   *  \param int : the number for this DynamicalSystem
   *  \param SiconosVector : initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param SiconosMatrix : mass matrix
   */
  LagrangianDS(int, const SiconosVector&, const SiconosVector&, const SiconosMatrix&);

  /** constructor from a minimum set of data
   *  \param int : the number for this DynamicalSystem
   *  \param SiconosVector : initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param string: plugin path to compute mass matrix
   */
  LagrangianDS(int, const SiconosVector& , const SiconosVector&, const std::string& = "DefaultPlugin:computeMass");

  /** destructor */
  virtual ~LagrangianDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

  /** initialization of qFree, vFree
   *  \param a string: the simulation type. For TimeStepping: memory allocation. For EventDriven: links (pointers) to q and velocity.
   */
  void initFreeVectors(const std::string&);

  /** allocate memory for p[...] vectors
   *  \param string: simulation type
   */
  void initP(const std::string&);

  /** allocate memory for fL and its jacobians, if required.
   */
  void initFL();

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization
   */
  virtual void initRhs(double) ;

  /** dynamical system initialization function: mainly set memory and compute plug-in for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const std::string&, double = 0, unsigned int = 1) ;

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
  inline void setNdof(unsigned int newNdof)
  {
    ndof = newNdof;
  };

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  inline const unsigned int getDim() const
  {
    return ndof;
  }

  // -- q --

  /** get the value of q
   *  \return SimpleVector
   */
  inline const SimpleVector getQ() const
  {
    return *q[0];
  }

  /** get q
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getQPtr() const
  {
    return q[0];
  }

  /** set the value of q to newValue
   *  \param SiconosVector newValue
   */
  void setQ(const SiconosVector&);

  /** set Q to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setQPtr(SiconosVector *newPtr);

  // -- q0 --

  /** get the value of q0
   *  \return SimpleVector
   */
  inline const SimpleVector getQ0() const
  {
    return *q0;
  }

  /** get q0
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getQ0Ptr() const
  {
    return q0;
  }

  /** set the value of q0 to newValue
   *  \param SiconosVector newValue
   */
  void setQ0(const SiconosVector&);

  /** set Q0 to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setQ0Ptr(SiconosVector *newPtr);

  // -- qFree --

  /** get the value of qFree[0]
   *  \return SimpleVector
   */
  inline const SimpleVector getQFree() const
  {
    return *workVector.find("qFree")->second;
  }

  /** get qFree[0]
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getQFreePtr() const
  {
    return workVector.find("qFree")->second;
  }

  /** set the value of qFree[0] to newValue
   *  \param SiconosVector newValue
   */
  void setQFree(const SiconosVector&);

  /** set QFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setQFreePtr(SiconosVector *newPtr);

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
    return *q[1];
  }

  /** get velocity
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getVelocityPtr() const
  {
    return q[1];
  }

  /** set the value of velocity to newValue
   *  \param SiconosVector newValue
   */
  void setVelocity(const SiconosVector&);

  /** set Velocity to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setVelocityPtr(SiconosVector *newPtr);

  // -- velocity0 --

  /** get the value of velocity0
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocity0() const
  {
    return *velocity0;
  }

  /** get velocity0
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getVelocity0Ptr() const
  {
    return velocity0;
  }

  /** set the value of velocity0 to newValue
   *  \param SiconosVector newValue
   */
  void setVelocity0(const SiconosVector&);

  /** set Velocity0 to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setVelocity0Ptr(SiconosVector *newPtr) ;

  // -- velocityFree --

  /** get the value of velocityFree
   *  \return SimpleVector
   */
  inline const SimpleVector getVelocityFree() const
  {
    return *workVector.find("velocityFree")->second;
  }

  /** get velocityFree
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getVelocityFreePtr() const
  {
    return workVector.find("velocityFree")->second;
  }

  /** set the value of velocityFree to newValue
   *  \param SiconosVector newValue
   */
  void setVelocityFree(const SiconosVector&);

  /** set VelocityFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setVelocityFreePtr(SiconosVector *newPtr);

  // -- acceleration --

  /** get acceleration
   *  \return pointer on a SiconosVector
   */
  SiconosVector* getAccelerationPtr() const ;

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
  inline const SimpleVector getP(unsigned int level = 2) const
  {
    return *(p[level]);
  }

  /** get p
   *  \param unsigned int, required level for p, default = 2
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getPPtr(unsigned int level = 2) const
  {
    return p[level];
  }

  /** set the value of p to newValue
   *  \param unsigned int, required level for p, default = 2
   *  \param SiconosVector newValue
   */
  void setP(const SiconosVector&, unsigned int level = 2);

  /** set P to pointer newPtr
   *  \param unsigned int, required level for p, default = 2
   *  \param SiconosVector * newPtr
   */
  void setPPtr(SiconosVector *newPtr, unsigned int level = 2);

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

  /** set the value of Mass to newValue
   *  \param SiconosMatrix newValue
   */
  void setMass(const SiconosMatrix&);

  /** set Mass to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setMassPtr(SiconosMatrix *newPtr);

  /** get MassLU: a copy of the mass matrix which is LU-factorized. Temporary function?
   *  \return a pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getMassLUPtr() const
  {
    return (workMatrix.find("invMass"))->second;
  }

  // -- FInt --

  /** get the value of fInt
   *  \return SimpleVector
   */
  inline const SimpleVector getFInt() const
  {
    return *fInt;
  }

  /** get fInt
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getFIntPtr() const
  {
    return fInt;
  }

  /** set the value of fInt to newValue
   *  \param SiconosVector newValue
   */
  void setFInt(const SiconosVector&);

  /** set FInt to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setFIntPtr(SiconosVector *newPtr);

  // -- Fext --

  /** get the value of fExt
   *  \return SimpleVector
   */
  inline const SimpleVector getFExt() const
  {
    return *fExt;
  }

  /** get fExt
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getFExtPtr() const
  {
    return fExt;
  }

  /** set the value of fExt to newValue
   *  \param SiconosVector newValue
   */
  void setFExt(const SiconosVector&);

  /** set FExt to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setFExtPtr(SiconosVector *newPtr);

  // -- NNL --

  /** get the value of NNL
   *  \return SimpleVector
   */
  inline const SimpleVector getNNL() const
  {
    return *NNL;
  }

  /** get NNL
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getNNLPtr() const
  {
    return NNL;
  }

  /** set the value of NNL to newValue
   *  \param SiconosVector newValue
   */
  void setNNL(const SiconosVector&);

  /** set NNL to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setNNLPtr(SiconosVector *newPtr);

  // -- Jacobian Fint --

  /** get the value of JacobianFInt[i]
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianFInt(unsigned int i) const
  {
    return *jacobianFInt[i];
  }

  /** get JacobianFInt[i]
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianFIntPtr(unsigned int i) const
  {
    return jacobianFInt[i];
  }

  /** set the value of JacobianFInt to newValue
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param SiconosMatrix newValue
   */
  void setJacobianFInt(unsigned int, const SiconosMatrix&);

  /** set JacobianFInt to pointer newPtr
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianFIntPtr(unsigned int, SiconosMatrix *newPtr);

  // -- Jacobian NNL --

  /** get the value of JacobianNNL
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianNNL(unsigned int i) const
  {
    return *jacobianNNL[i];
  }

  /** get JacobianNNL
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianNNLPtr(unsigned int i) const
  {
    return jacobianNNL[i];
  }

  /** set the value of JacobianNNL to newValue
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param SiconosMatrix newValue
   */
  void setJacobianNNL(unsigned int, const SiconosMatrix&);

  /** set JacobianNNL to pointer newPtr
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param SiconosMatrix * newPtr
   */
  void setJacobianNNLPtr(unsigned int, SiconosMatrix *newPtr);

  // -- fL --

  /** get the value of fL
   *  \return SimpleVector
   */
  inline const SimpleVector getFL() const
  {
    return *fL;
  }

  /** get fL
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getFLPtr() const
  {
    return fL;
  }

  // -- Jacobian fL --

  /** get the value of JacobianFL
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianFL(unsigned int i) const
  {
    return *jacobianFL[i];
  }

  /** get JacobianFL
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianFLPtr(unsigned int i) const
  {
    return jacobianFL[i];
  }

  // --- PLUGINS RELATED FUNCTIONS ---

  /** allow to set a specified function to compute the mass
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeMassFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName);

  /** allow to set a specified function to compute the inertia
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeNNLFunction(const std::string& pluginPath, const std::string&  functionName);

  /** allow to set a specified function to compute the gradient of the internal strength compared to the state
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeJacobianFIntFunction(unsigned int, const std::string&  pluginPath, const std::string&  functionName);

  /** allow to set a specified function to compute the gradient of the the external strength compared to the state
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeJacobianNNLFunction(unsigned int, const std::string&  pluginPath, const std::string&  functionName);
  void setJacobianNNLFunction(unsigned int, FPtr5);

  /** default function to compute the mass
   */
  void computeMass();

  /** function to compute the mass
   *  \param double time : the current time, SiconosVector*: pointer on the state vector q
   */
  void computeMass(SiconosVector *);

  /** default function to compute the internal strengths
   *  \param double time : the current time
   */
  void computeFInt(double);

  /** function to compute the internal strengths
   *  with some specific values for q and velocity (ie not those of the current state).
   *  \param double time : the current time, SiconosVector*: pointers on the state vectors q and velocity
   */
  void computeFInt(double , SiconosVector *, SiconosVector *);

  /** default function to compute the external strengths
   *  \param double time : the current time
   */
  void computeFExt(double);

  /** default function to compute the inertia
   */
  void computeNNL();

  /** function to compute the inertia
   *  with some specific values for q and velocity (ie not those of the current state).
   *  \param SiconosVector*: pointers on the state vectors q and velocity
   */
  void computeNNL(SiconosVector *q, SiconosVector *velocity);

  /** default function to compute the gradient of the internal strengths compared to the state
   *  \param double time : the current time
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   */
  void computeJacobianFInt(unsigned int, double);

  /** function to compute the gradient of the internal strengths compared to state q
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param double time : the current time, SiconosVector*: pointers on the state vectors q and velocity
   */
  void computeJacobianFInt(unsigned int, double , SiconosVector *q, SiconosVector *velocity);

  /** function to compute the gradient of the inertia strengths compared to the state q
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   */
  void computeJacobianNNL(unsigned int);

  /** function to compute the gradient of the inertia strengths compared to the state q
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param SiconosVector*: pointers on the state vectors q and velocity
   */
  void computeJacobianNNL(unsigned int, SiconosVector *q, SiconosVector *velocity);

  /** Default function to compute the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeRhs(double, bool  = false);

  /** Default function to compute jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeJacobianXRhs(double, bool  = false);

  /** Default function to compute fL
   *  \param double, the current time
   */
  virtual void computeFL(double time);

  /** function to compute fL with some specific values for q and velocity (ie not those of the current state).
   *  \param double time : the current time
   *  \param SiconosVector*: pointers on q
   *  \param SiconosVector*: pointers on velocity
   */
  virtual void computeFL(double , SiconosVector *, SiconosVector *);

  /** Default function to compute the jacobian of fL
   *  \param index (0: \f$ \nabla_q \f$, 1: \f$ \nabla_{\dot q} \f$ )
   *  \param double, the current time
   */
  virtual void computeJacobianFL(unsigned int, double);

  // --- miscellaneous ---

  /** copy the data of the DS into the XML tree
   */
  virtual void saveSpecificDataToXML();

  /** print the data to the screen
   */
  virtual void display() const;

  /** initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int steps);

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
   *  \param SiconosVector*: in-out parameter, qFree
   */
  void computeQFree(double, unsigned int, SiconosVector*);

  /** set p[...] to zero
   */
  void resetNonSmoothPart();

  /** Computes post-impact velocity, using pre-impact velocity and impulse (p) value.
   * Used in EventDriven (Lsodar->updateState)
   */
  void computePostImpactVelocity();

};

#endif // LAGRANGIANNLDS_H

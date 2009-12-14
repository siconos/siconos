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

/*! \file LagrangianDS.h
  \brief LagrangianDS class - Second Order Non Linear Dynamical Systems.

*/

#ifndef LAGRANGIANNLDS_H
#define LAGRANGIANNLDS_H

#include "DynamicalSystem.hpp"
#include "Plugin.hpp"

class DynamicalSystem;
/** Pointer to function for plug-in. For NNL and its jacobian. */
typedef void (*FPtr5)(unsigned int, const double*, const double*, double*, unsigned int, double*);
typedef void (*FPtrMass)(unsigned int, const double*, double*, unsigned int, double*);
typedef  void (*FPtrFExt)(double, unsigned int, double*, unsigned int, double*);

/** Lagrangian non linear dynamical systems - Derived from DynamicalSystem -
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
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
 *
 *
 */
class LagrangianDS : public DynamicalSystem
{
public:

  /** List of indices used to save tmp work matrics (last one is the size of the present list) */
  enum WorkMatrixNames {invMass, jacobianXBloc10, jacobianXBloc11, zeroMatrix, idMatrix, sizeWorkMat};

protected:

  // -- MEMBERS --

  /** number of degrees of freedom of the system */
  unsigned int _ndof;

  /** state of the system. See details on top of page. */
  VectorOfVectors _q;

  /** initial coordinates of the system */
  SP::SiconosVector _q0;

  /** initial velocity of the system */
  SP::SiconosVector _velocity0;

  /** memory of previous coordinates of the system */
  SP::SiconosMemory _qMemory;

  /** memory of previous velocities of the system */
  SP::SiconosMemory _velocityMemory;

  /** "Reaction" due to the non smooth law - The index corresponds to the dynamic levels. */
  std::vector<SP::SiconosVector> _p;

  /** mass of the system */
  SP::SiconosMatrix _mass;

  /** internal strength of the system */
  SP::SiconosVector _fInt;

  /** jacobian_q Fint*/
  SP::SiconosMatrix _jacobianQFInt;
  /** jacobian_{qDot} Fint*/
  SP::SiconosMatrix _jacobianQDotFInt;

  /** external strength of the system */
  SP::SiconosVector _fExt;

  /** non-linear inertia term of the system */
  SP::SiconosVector _NNL;

  /** jacobian_q QNNL*/
  SP::SiconosMatrix _jacobianQNNL;
  /** jacobian_{qDot} QNNL*/
  SP::SiconosMatrix _jacobianQDotNNL;

  /** fL(q[0],q[1],t)= fExt - fInt -NNL */
  SP::SiconosVector _fL;

  /** jacobian_q FL*/
  SP::SiconosMatrix _jacobianQFL;
  /** jacobian_{qDot} FL*/
  SP::SiconosMatrix _jacobianQDotFL;

  /** set links with DS members
   */
  void connectToDS();

  /** Default constructor
   */
  LagrangianDS();


  // pointers to functions member to compute plug-in functions

  /** LagrangianDS plug-in to compute mass(q,t) - id = "mass"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param[in,out] mass : pointer to the first element of mass
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  //  void (*computeMassPtr)(unsigned int, const double*, double*, unsigned int, double*);
  SP::PluggedObject _pluginMass;


  /** LagrangianDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] fInt : pointer to the first element of fInt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  SP::PluggedObject _pluginFInt;
  //  FPtr6 computeFIntPtr;

  /** LagrangianDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param[in,out] fExt : pointer to the first element of fExt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  //  void (*computeFExtPtr)(double, unsigned int, double*, unsigned int, double* );
  SP::PluggedObject _pluginFExt;

  /** LagrangianDS plug-in to compute \f$NNL(\dot q, q)\f$, id = "NNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] NNL : pointer to the first element of NNL
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr5 computeNNLPtr;
  SP::PluggedObject _pluginNNL;

  /** LagrangianDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianQFInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr6 computeJacobianQFIntPtr;
  SP::PluggedObject _pluginJacQFInt;

  /** LagrangianDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianQDotFInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr6 computeJacobianQDotFIntPtr;
  SP::PluggedObject _pluginJacQDotFInt;

  /** LagrangianDS plug-in to compute \f$\nabla_qNNL(\dot q, q)\f$, id = "jacobianQNNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr5 computeJacobianQNNLPtr;
  SP::PluggedObject _pluginJacQNNL;
  /** LagrangianDS plug-in to compute \f$\nabla_{\dot q}NNL(\dot q, q)\f$, id = "jacobianQDotNNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  //  FPtr5 computeJacobianQDotNNLPtr;
  SP::PluggedObject _pluginJacQDotNNL;

  virtual void zeroPlugin();
public:

  // === CONSTRUCTORS - DESTRUCTOR ===

  /** constructor from a minimum set of data
   *  \param SiconosVector : initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   */
  LagrangianDS(SP::SiconosVector, SP::SiconosVector);

  /** constructor from an xml file
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   */
  LagrangianDS(SP::DynamicalSystemXML);

  /** constructor from a minimum set of data
   *  \param SiconosVector : initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param SiconosMatrix : mass matrix
   */
  LagrangianDS(SP::SiconosVector, SP::SiconosVector, SP::SiconosMatrix);

  /** constructor from a minimum set of data
   *  \param SiconosVector : initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param string: plugin path to compute mass matrix
   */
  LagrangianDS(SP::SiconosVector , SP::SiconosVector, const std::string&);

  /** destructor */
  virtual ~LagrangianDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

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
  void initRhs(double) ;

  /** dynamical system initialization function: mainly set memory and compute plug-in for initial state values.
   *  \param string: simulation type
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  void initialize(const std::string&, double = 0, unsigned int = 1) ;

  // === GETTERS AND SETTERS ===

  /** to get the value of ndof
   *  \return the value of ndof
   */
  inline const unsigned int getNdof() const
  {
    return _ndof;
  };

  /** to set ndof
   *  \param unsigned int ndof : the value to set ndof
   */
  inline void setNdof(unsigned int newNdof)
  {
    _ndof = newNdof;
  };

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline const unsigned int getDim() const
  {
    return _ndof;
  }

  // -- q --

  /** get q
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const
  {
    return _q[0];
  }

  /** set the value of q to newValue
   *  \param SiconosVector newValue
   */
  void setQ(const SiconosVector&);

  /** set Q to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setQPtr(SP::SiconosVector newPtr);

  // -- q0 --

  /** get q0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q0() const
  {
    return _q0;
  }

  /** set the value of q0 to newValue
   *  \param SiconosVector newValue
   */
  void setQ0(const SiconosVector&);

  /** set Q0 to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setQ0Ptr(SP::SiconosVector newPtr);

  // Q memory

  /** get all the values of the state vector q stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory qMemory() const
  {
    return _qMemory;
  }

  // -- velocity --

  /** get velocity
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity() const
  {
    return _q[1];
  }

  /** set the value of velocity to newValue
   *  \param SiconosVector newValue
   */
  void setVelocity(const SiconosVector&);

  /** set Velocity to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setVelocityPtr(SP::SiconosVector newPtr);

  // -- velocity0 --

  /** get velocity0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity0() const
  {
    return _velocity0;
  }

  /** set Velocity0 to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setVelocity0Ptr(SP::SiconosVector newPtr) ;

  // -- acceleration --

  /** get acceleration
   *  \return pointer on a SiconosVector
   */
  SP::SiconosVector acceleration() const ;

  // Velocity memory

  /** get all the values of the state vector velocity stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory velocityMemory() const
  {
    return _velocityMemory;
  }

  // -- p --

  /** get p
   *  \param unsigned int, required level for p, default = 2
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector p(unsigned int level = 2) const
  {
    return _p[level];
  }

  /** set the value of p to newValue
   *  \param unsigned int, required level for p, default = 2
   *  \param SiconosVector newValue
   */
  void setP(const SiconosVector&, unsigned int level = 2);

  /** set P to pointer newPtr
   *  \param unsigned int, required level for p, default = 2
   *  \param SP::SiconosVector newPtr
   */
  void setPPtr(SP::SiconosVector newPtr, unsigned int level = 2);

  // -- Mass --

  /** get mass
   *  \return pointer on a plugged-matrix
   */
  inline SP::SiconosMatrix mass() const
  {
    return _mass;
  }

  /** set mass to pointer newPtr
   *  \param a plugged matrix SP
   */
  inline void setMassPtr(SP::SiconosMatrix newPtr)
  {
    _mass = newPtr;
  }

  /** get MassLU: a copy of the mass matrix which is LU-factorized. Temporary function?
   *  \return a pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix massLU() const
  {
    return (_workMatrix[invMass]);
  }

  // --- fInt ---

  /** get fInt
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fInt() const
  {
    return _fInt;
  }


  /** set fInt to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setFIntPtr(SP::SiconosVector newPtr)
  {
    _fInt = newPtr;
  }

  // -- Fext --

  /** get fExt
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector fExt() const
  {
    return _fExt;
  }


  /** set fExt to pointer newPtr
   *  \param a SP to a Simple vector
   */
  inline void setFExtPtr(SP::SimpleVector newPtr)
  {
    _fExt = newPtr;
  }

  // -- NNL --

  /** get NNL
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector NNL() const
  {
    return _NNL;
  }


  /** set NNL to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setNNLPtr(SP::SiconosVector newPtr)
  {
    _NNL = newPtr;
  }


  // -- Jacobian Fint --

  /** get jacobianQFInt
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqFInt() const
  {
    return _jacobianQFInt;
  }
  /** get jacobianQDotFInt
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianqDotFInt() const
  {
    return _jacobianQDotFInt;
  }
  //  inline SP::SiconosMatrix jacobianZFInt() const { return jacobianZFInt; }


  /** set jacobianQFInt to pointer newPtr
   *  \param a SP SiconosMatrix
   */
  inline void setJacobianQFIntPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianQFInt = newPtr;
  }
  /** set jacobianQDotFInt to pointer newPtr
   *  \param a SP SiconosMatrix
   */
  inline void setJacobianQDotFIntPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianQDotFInt = newPtr;
  }

  // -- Jacobian NNL --


  /** get jacobianQNNL
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianQNNL() const
  {
    return _jacobianQNNL;
  }
  /** get jacobianQDotNNL
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianQDotNNL() const
  {
    return _jacobianQDotNNL;
  }


  /** set jacobianQNNL to pointer newPtr
   *  \param a SP SiconosMatrix
   */
  inline void setJacobianQNNLPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianQNNL = newPtr;
  }
  /** set jacobianQDotNNL to pointer newPtr
   *  \param a SP SiconosMatrix
   */
  inline void setJacobianQDotNNLPtr(SP::SiconosMatrix newPtr)
  {
    _jacobianQDotNNL = newPtr;
  }

  // -- fL --

  /** get the value of fL
   *  \return SimpleVector
   */
  inline const SimpleVector getFL() const
  {
    return *_fL;
  }

  /** get fL
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector fL() const
  {
    return _fL;
  }

  // -- Jacobian fL --


  /** get JacobianFL
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianQFL() const
  {
    return _jacobianQFL;
  }
  /** get JacobianFL
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianQDotFL() const
  {
    return _jacobianQDotFL;
  }
  //  inline SP::SiconosMatrix jacobianZFL() const { return jacobianZFL; }

  // --- PLUGINS RELATED FUNCTIONS ---

  /** allow to set a specified function to compute the mass
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeMassFunction(const std::string&  pluginPath, const std::string&  functionName)
  {
    _pluginMass->setComputeFunction(pluginPath, functionName);
    //     Plugin::setFunction(&computeMassPtr, pluginPath,functionName);
  }

  /** set a specified function to compute Mass
   *  \param a pointer on the plugin function
   */
  void setComputeMassFunction(FPtr7 fct)
  {
    _pluginMass->setComputeFunction((void*)fct);
    //    computeMassPtr=fct;
  }

  /** allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFIntFunction(const std::string&  pluginPath, const std::string&  functionName)
  {
    _pluginFInt->setComputeFunction(pluginPath, functionName);
    //    Plugin::setFunction(&computeFIntPtr, pluginPath,functionName);
  }

  /** set a specified function to compute fInt
   *  \param a pointer on the plugin function
   */
  void setComputeFIntFunction(FPtr6 fct)
  {
    _pluginFInt->setComputeFunction((void*)fct);
    //    computeFIntPtr = fct;
  }

  /** allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    _pluginFExt->setComputeFunction(pluginPath, functionName);
    //    Plugin::setFunction(&computeFExtPtr, pluginPath,functionName);
  }

  /** set a specified function to compute fExt
   *  \param a pointer on the plugin function
   */
  void setComputeFExtFunction(VectorFunctionOfTime fct)
  {
    _pluginFExt->setComputeFunction((void*)fct);
    //   computeFExtPtr = fct ;
  }

  /** allow to set a specified function to compute the inertia
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeNNLFunction(const std::string& pluginPath, const std::string&  functionName);

  /** set a specified function to compute NNL
   *  \param a pointer on the plugin function
   */
  void setComputeNNLFunction(FPtr5 fct);

  /** allow to set a specified function to compute the jacobian following q of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeJacobianQFIntFunction(const std::string&  pluginPath, const std::string&  functionName);
  /** allow to set a specified function to compute the jacobian following qDot of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeJacobianQDotFIntFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *  \param a pointer on the plugin function
   */
  void setComputeJacobianQFIntFunction(FPtr6 fct);
  /** set a specified function to compute jacobian following qDot of the FInt
   *  \param a pointer on the plugin function
   */
  void setComputeJacobianQDotFIntFunction(FPtr6 fct);

  /** allow to set a specified function to compute the jacobian following q of the the external strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeJacobianQNNLFunction(const std::string&  pluginPath, const std::string&  functionName);
  /** allow to set a specified function to compute the jacobian following qDot of the the external strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeJacobianQDotNNLFunction(const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute the jacobian following q of NNL
   *  \param a pointer on the plugin function
   */
  void setComputeJacobianQNNLFunction(FPtr5 fct);
  /** set a specified function to compute the jacobian following qDot of NNL
   *  \param a pointer on the plugin function
   */
  void setComputeJacobianQDotNNLFunction(FPtr5 fct);

  /** default function to compute the mass
   */
  virtual void computeMass();

  /** function to compute the mass
   *  \param double time : the current time, SP::SiconosVector: pointer on the state vector q
   */
  virtual void computeMass(SP::SiconosVector);

  /** default function to compute the internal strengths
   *  \param double time : the current time
   */
  virtual void computeFInt(double);

  /** function to compute the internal strengths
   *  with some specific values for q and velocity (ie not those of the current state).
   *  \param double time : the current time, SP::SiconosVector: pointers on the state vectors q and velocity
   */
  virtual void computeFInt(double , SP::SiconosVector, SP::SiconosVector);

  /** default function to compute the external strengths
   *  \param double time : the current time
   */
  virtual void computeFExt(double);

  /** default function to compute the inertia
   */
  virtual void computeNNL();

  /** function to compute the inertia
   *  with some specific values for q and velocity (ie not those of the current state).
   *  \param SP::SiconosVector: pointers on the state vectors q and velocity
   */
  virtual void computeNNL(SP::SiconosVector q, SP::SiconosVector velocity);

  /** To compute the jacobian following q of the internal strengths compared to the state
   *  \param double time : the current time
   */
  virtual void computeJacobianQFInt(double);
  /** To compute the jacobian following qDot of the internal strengths compared to the state
   *  \param double time : the current time
   */
  virtual void computeJacobianQDotFInt(double);

  /** To compute the jacobian following q of the internal strengths compared to state q
   *  \param double time : the current time, SP::SiconosVector: pointers on the state vectors q and velocity
   */
  virtual void computeJacobianQFInt(double , SP::SiconosVector q, SP::SiconosVector velocity);
  /** To compute the jacobian following qDot of the internal strengths compared to state q
   *  \param double time : the current time, SP::SiconosVector: pointers on the state vectors q and velocity
   */
  virtual void computeJacobianQDotFInt(double , SP::SiconosVector q, SP::SiconosVector velocity);

  /** function to compute the jacobian following q of the inertia strengths compared to the state q
   */
  virtual void computeJacobianQNNL();
  /** function to compute the jacobian following qDot of the inertia strengths compared to the state q
   */
  virtual void computeJacobianQDotNNL();

  /** function to compute the jacobian following q of the inertia strengths compared to the state q
   *  \param SP::SiconosVector: pointers on the state vectors q and velocity
   */
  virtual void computeJacobianQNNL(SP::SiconosVector q, SP::SiconosVector velocity);
  /** function to compute the jacobian following qDot of the inertia strengths compared to the state q
   *  \param SP::SiconosVector: pointers on the state vectors q and velocity
   */
  virtual void computeJacobianQDotNNL(SP::SiconosVector q, SP::SiconosVector velocity);

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
  virtual void computeFL(double);

  /** function to compute fL with some specific values for q and velocity (ie not those of the current state).
   *  \param double time : the current time
   *  \param SP::SiconosVector: pointers on q
   *  \param SP::SiconosVector: pointers on velocity
   */
  virtual void computeFL(double , SP::SiconosVector, SP::SiconosVector);

  /** Default function to compute the jacobian following q of fL
   *  \param double, the current time
   */
  virtual void computeJacobianQFL(double);
  /** Default function to compute the jacobian following qDot of fL
   *  \param double, the current time
   */
  virtual void computeJacobianQDotFL(double);

  // --- miscellaneous ---

  /** copy the data of the DS into the XML tree
   */
  void saveSpecificDataToXML();

  /** print the data to the screen
   */
  void display() const;

  /** initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  void initMemory(unsigned int);

  /** push the current values of x, q and r in the stored previous values
   *  xMemory, qMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory();

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param SP::DynamicalSystem : the system which must be converted
   * \return a pointer on the system if it is of the right type, NULL otherwise
   */
  static LagrangianDS* convert(DynamicalSystem* ds);

  /** To compute \f$\frac{|q_{i+1} - qi|}{|q_i|}\f$ where \f$ q_{i+1}\f$ represents the present state and \f$ q_i\f$ the previous one
   * \return a double
   */
  /*  double dsConvergenceIndicator(); */

  /** function to compute derivative number level of qFree
   *  \param double: current time
   *  \param unsigned int: derivative number
   *  \param SP::SiconosVector: in-out parameter, qFree
   */
  //  void computeQFree(double, unsigned int, SP::SiconosVector);

  /** set p[...] to zero
   */
  void resetNonSmoothPart();

  /** Computes post-impact velocity, using pre-impact velocity and impulse (p) value.
   * Used in EventDriven (Lsodar->updateState)
   */
  void computePostImpactVelocity();

};

TYPEDEF_SPTR(LagrangianDS);

#endif // LAGRANGIANNLDS_H

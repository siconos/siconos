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

/*! \file NewtonEulerDS.h
  \brief NewtonEulerDS class - Second Order Non Linear Dynamical Systems.

*/

#ifndef NEWTONEULERNLDS_H
#define NEWTONEULERNLDS_H

#include "DynamicalSystem.hpp"
#include "Plugin.hpp"

class DynamicalSystem;
/** Pointer to function for plug-in. For NNL and its jacobian. */
typedef void (*FPtr5)(unsigned int, const double*, const double*, double*, unsigned int, double*);

/** NewtonEuler non linear dynamical systems - Derived from DynamicalSystem -
 *
 *
 */
class NewtonEulerDS : public DynamicalSystem
{
protected:

  // -- MEMBERS --

  /** _v conatains the speed of the Newton Euler system:
      _v[0:2] : velocity of the mass center.
      _v[3:5] : angular velocity in the referencial attached to the object.
  */
  SP::SiconosVector _v;
  SP::SiconosVector _v0;


  /** _vMemory: to acces at previous step */
  SP::SiconosMemory _vMemory;
  SP::SiconosMemory _qMemory;
  SP::SiconosMemory _fLMemory;
  SP::SiconosMemory _dotqMemory;

  /** space dimention (1,2 or 3) */
  unsigned int _ndof;


  /** _q dimension, is not necessary _n. In our case, _qDim = 7. _n =6*/
  int _qDim;

  /** _q contains the representation of the system, in current implementation:
      _q[0:2] : the coordinates of the center of mass.
      _q[3:6] : a quaternion representing the orientation of the solid.
  */
  SP::SiconosVector _q;
  SP::SiconosVector _q0;
  /** The time derivative of q*/
  SP::SiconosVector _dotq;

  /** Inertial matrix*/
  SP::SiconosMatrix _I;
  /** mass of the system */
  double _mass;


  /** used for concatenate _I and _mass*/
  SP::SimpleMatrix _W;

  SP::SimpleMatrix _luW;

  /** Matrix depending of the meaning of x.*/
  SP::SiconosMatrix _T;


  /** "Reaction" due to the non smooth law - The index corresponds to the dynamic levels. */
  std::vector<SP::SiconosVector> _p;


  /** extarnal moment of the forces */
  SP::SiconosVector _mExt;

  /** jacobian_q */
  SP::SiconosMatrix _jacobianQmExt;
  /** jacobian_{qDot} */
  SP::SiconosMatrix _jacobianQDotmExt;

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
  NewtonEulerDS();


  // pointers to functions member to compute plug-in functions

  /** NewtonEulerDS plug-in to compute internal forces \f$F_{int}(t,q,\dot q)\f$ - id = "fInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] fInt : pointer to the first element of fInt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  FPtr6 computeFIntPtr;

  /** NewtonEulerDS plug-in to compute external forces \f$F_{Ext}(t)\f$, id = "fExt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param[in,out] fExt : pointer to the first element of fExt
   * @param  size of vector z
   * @param[in,out] z : a vector of user-defined parameters
   */
  void (*computeFExtPtr)(double, unsigned int, double*, unsigned int, double*);
  void (*computeMmtExtPtr)(double, unsigned int, double*, unsigned int, double*);

  /** NewtonEulerDS plug-in to compute \f$NNL(\dot q, q)\f$, id = "NNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] NNL : pointer to the first element of NNL
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr5 computeNNLPtr;

  /** NewtonEulerDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianQFInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr6 computeJacobianQFIntPtr;

  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianQDotFInt"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr6 computeJacobianQDotFIntPtr;

  /** NewtonEulerDS plug-in to compute \f$\nabla_qNNL(\dot q, q)\f$, id = "jacobianQNNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr5 computeJacobianQNNLPtr;
  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}NNL(\dot q, q)\f$, id = "jacobianQDotNNL"
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr5 computeJacobianQDotNNLPtr;

  void zeroPlugin();

public:

  // === CONSTRUCTORS - DESTRUCTOR ===


  /** constructor from a minimum set of data
   *  \param SiconosVector : initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param double : mass
   *  \param SiconosMatrix : inertial matrix
   */
  NewtonEulerDS(SP::SiconosVector, SP::SiconosVector, double  , SP::SiconosMatrix);



  /** destructor */
  virtual ~NewtonEulerDS();

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

  /** return the dim. of the system (n for first order). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline const unsigned int getDim() const
  {
    return _n;
  }

  // -- q --

  /** get q
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q() const
  {
    return _q;
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
  inline SP::SiconosMemory vMemory() const
  {
    return _vMemory;
  }

  // -- velocity --

  /** get velocity
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity() const
  {
    return _v;
  }



  // -- velocity0 --

  /** get velocity0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector velocity0() const
  {
    return _v0;
  }

  /** set Velocity0 to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setVelocity0Ptr(SP::SiconosVector newPtr) ;

  // -- acceleration --

  /** get acceleration
   *  \return pointer on a SiconosVector
   */
  //  SP::SiconosVector acceleration() const ;

  // Velocity memory

  /** get all the values of the state vector velocity stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory velocityMemory() const
  {
    return _vMemory;
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
  //   inline SP::SiconosMatrix mass() const { return _mass; }

  //   /** set mass to pointer newPtr
  //    *  \param a plugged matrix SP
  //    */
  //   inline void setMassPtr(SP::SiconosMatrix newPtr) {_mass = newPtr;}

  //   /** get MassLU: a copy of the mass matrix which is LU-factorized. Temporary function?
  //    *  \return a pointer on a SiconosMatrix
  //    */
  //   inline SP::SiconosMatrix massLU() const { return (_workMatrix[invMass]); }

  // --- fInt ---

  /** get fInt
   *  \return pointer on a plugged vector
   */
  //  inline SP::SiconosVector fInt() const { return _fInt; }


  /** set fInt to pointer newPtr
   *  \param a SP to plugged vector
   */
  //  inline void setFIntPtr(SP::SiconosVector newPtr) {_fInt = newPtr;}

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
  //  inline SP::SiconosMatrix jacobianqFInt() const { return _jacobianQFInt; }
  /** get jacobianQDotFInt
   *  \return pointer on a SiconosMatrix
   */
  //  inline SP::SiconosMatrix jacobianqDotFInt() const { return _jacobianQDotFInt; }
  //  inline SP::SiconosMatrix jacobianZFInt() const { return jacobianZFInt; }


  /** set jacobianQFInt to pointer newPtr
   *  \param a SP SiconosMatrix
   */
  //  inline void setJacobianQFIntPtr( SP::SiconosMatrix newPtr) {_jacobianQFInt = newPtr;}
  /** set jacobianQDotFInt to pointer newPtr
   *  \param a SP SiconosMatrix
   */
  //  inline void setJacobianQDotFIntPtr( SP::SiconosMatrix newPtr) {_jacobianQDotFInt = newPtr;}

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
  //   void setComputeMassFunction(const std::string&  pluginPath, const std::string&  functionName){
  //      Plugin::setFunction(&computeMassPtr, pluginPath,functionName);
  //   }

  //   /** set a specified function to compute Mass
  //    *  \param a pointer on the plugin function
  //    */
  //   void setComputeMassFunction(FPtr7 fct ){
  //     computeMassPtr=fct;
  //   }

  /** allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  //   void setComputeFIntFunction(const std::string&  pluginPath, const std::string&  functionName){
  //     Plugin::setFunction(&computeFIntPtr, pluginPath,functionName);
  //   }

  /** set a specified function to compute fInt
   *  \param a pointer on the plugin function
   */
  //   void setComputeFIntFunction(FPtr6 fct ){
  //     computeFIntPtr = fct;
  //   }

  /** allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    Plugin::setFunction(&computeFExtPtr, pluginPath, functionName);
  }

  /** set a specified function to compute fExt
   *  \param a pointer on the plugin function
   */
  void setComputeFExtFunction(VectorFunctionOfTime fct)
  {
    computeFExtPtr = fct ;
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
  //  void setComputeJacobianQFIntFunction( const std::string&  pluginPath, const std::string&  functionName);
  /** allow to set a specified function to compute the jacobian following qDot of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  //  void setComputeJacobianQDotFIntFunction( const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobian following q of the FInt
   *  \param a pointer on the plugin function
   */
  //  void setComputeJacobianQFIntFunction(FPtr6 fct);
  /** set a specified function to compute jacobian following qDot of the FInt
   *  \param a pointer on the plugin function
   */
  //  void setComputeJacobianQDotFIntFunction(FPtr6 fct);

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
  //  virtual void computeMass();

  /** function to compute the mass
   *  \param double time : the current time, SP::SiconosVector: pointer on the state vector q
   */
  //  virtual void computeMass(SP::SiconosVector);

  /** default function to compute the internal strengths
   *  \param double time : the current time
   */
  //  virtual void computeFInt(double );

  /** function to compute the internal strengths
   *  with some specific values for q and velocity (ie not those of the current state).
   *  \param double time : the current time, SP::SiconosVector: pointers on the state vectors q and velocity
   */
  //  virtual void computeFInt(double , SP::SiconosVector, SP::SiconosVector);

  /** default function to compute the external strengths
   *  \param double time : the current time
   */
  virtual void computeFExt(double);
  virtual void computeMmtExt(double);

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
  //  virtual void computeJacobianQFInt(double);
  /** To compute the jacobian following qDot of the internal strengths compared to the state
   *  \param double time : the current time
   */
  //  virtual void computeJacobianQDotFInt(double);

  /** To compute the jacobian following q of the internal strengths compared to state q
   *  \param double time : the current time, SP::SiconosVector: pointers on the state vectors q and velocity
   */
  //  virtual void computeJacobianQFInt( double , SP::SiconosVector q, SP::SiconosVector velocity);
  /** To compute the jacobian following qDot of the internal strengths compared to state q
   *  \param double time : the current time, SP::SiconosVector: pointers on the state vectors q and velocity
   */
  //  virtual void computeJacobianQDotFInt( double , SP::SiconosVector q, SP::SiconosVector velocity);

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
  //virtual void computeFL(double , SP::SiconosVector, SP::SiconosVector);

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
  static NewtonEulerDS* convert(DynamicalSystem* ds);

  /** To compute \f$\frac{|q_{i+1} - qi|}{|q_i|}\f$ where \f$ q_{i+1}\f$ represents the present state and \f$ q_i\f$ the previous one
   * \return a double
   */
  /*  double dsConvergenceIndicator(); */

  /** function to compute derivative number level of qFree
   *  \param double: current time
   *  \param unsigned int: derivative number
   *  \param SP::SiconosVector: in-out parameter, qFree
   */
  void computeQFree(double, unsigned int, SP::SiconosVector);

  /** set p[...] to zero
   */
  void resetNonSmoothPart();


  virtual void updateT();

  inline SP::SiconosMatrix W()
  {
    return _W;
  }
  inline SP::SiconosMatrix luW()
  {
    return _luW;
  }
  inline SP::SiconosMatrix T()
  {
    return _T;
  }
  inline SP::SiconosMemory fLMemory()
  {
    return _fLMemory;
  }
  inline SP::SiconosMemory dotqMemory()
  {
    return _dotqMemory;
  }
  inline SP::SiconosVector dotq()
  {
    return _dotq;
  }



};

TYPEDEF_SPTR(NewtonEulerDS);

#endif // NEWTONEULERNLDS_H

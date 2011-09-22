/* Siconos-Kernel, Copyright INRIA 2005-2011.
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
 * Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
 */

/*! \file NewtonEulerDS.hpp
  \brief NewtonEulerDS class - Second Order Non Linear Dynamical Systems.

  For more details, see the DevNotes.pdf, chapter NewtonEuler.
*/
#ifndef NEWTONEULERNLDS_H
#define NEWTONEULERNLDS_H

#include "DynamicalSystem.hpp"
#include "Plugin.hpp"

class DynamicalSystem;
/** Pointer to function for plug-in. */
typedef void (*FPtr5)(unsigned int, const double*, const double*, double*, unsigned int, double*);
typedef void (*Fext)(double , double*, double*, double*);
/** NewtonEuler non linear dynamical systems - Derived from DynamicalSystem -
 *
 *
 */
class NewtonEulerDS : public DynamicalSystem
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(NewtonEulerDS);

  void internalInit(SP::SiconosVector Q0, SP::SiconosVector Velocity0, double mass , SP::SiconosMatrix inertialMatrix);

  // -- MEMBERS --

  /** _v contains the speed of the Newton Euler system:
      _v[0:2] : velocity of the mass center.
      _v[3:5] : Omega, angular velocity in the referencial attached to the object.
  */
  SP::SiconosVector _v;
  SP::SiconosVector _v0;


  /** _vMemory: to acces at previous step */
  SP::SiconosMemory _vMemory;
  SP::SiconosMemory _qMemory;
  SP::SiconosMemory _forcesMemory;
  SP::SiconosMemory _dotqMemory;


  /** _q dimension, is not necessary _n. In our case, _qDim = 7. _n =6*/
  unsigned int _qDim;

  /** _q contains the representation of the system, in current implementation:
      _q[0:2] : the coordinates of the center of mass.
      _q[3:6] : a quaternion representing the orientation of the solid.
  */
  SP::SimpleVector _q;
  SP::SimpleVector _deltaq;
  SP::SiconosVector _q0;
  /** The time derivative of q*/
  SP::SiconosVector _dotq;
  /*Matrix converting the object coordinates in the absolute coordinates.*/
  SP::SimpleMatrix _MObjToAbs;

  /** Inertial matrix*/
  SP::SiconosMatrix _I;
  /** mass of the system */
  double _mass;


  /** used for concatenate _I and _mass.I_3 */
  SP::SimpleMatrix _massMatrix;

  SP::SimpleMatrix _luW;

  /** Matrix depending of the meaning of x.*/
  SP::SimpleMatrix _T;


  /** "Reaction" due to the non smooth law - The index corresponds to the dynamic levels. */
  std::vector<SP::SiconosVector> _p;


  /** external moment of the forces */
  SP::SiconosVector _mExt;

  /*The following code is commented because the jacobian of m/fext are not yet integrated in the numerical scheme.*/
  /** jacobian_q */
  //  SP::SiconosMatrix _jacobianqmExt;
  /** jacobian_{qDot} */
  //  SP::SiconosMatrix _jacobianqDotmExt;

  /** external strength of the system */
  SP::SiconosVector _fExt;

  /** forces(q[0],q[1],t)= fExt - fInt */
  SP::SiconosVector _forces;

  /** jacobian_q FL*/
  SP::SimpleMatrix _jacobianvFL;
  /** jacobian_{qDot} FL*/
  SP::SimpleMatrix _jacobianqDotForces;

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
   * @param[in] time : current time
   * @param[in] q : current dof
   * @param[in,out] fExt : pointer to the first element of fExt
   * @param[in,out] z : a vector of user-defined parameters
   */
  void (*computeFExtPtr)(double, double *, double*, double*);
  void (*computeMExtPtr)(double, double *, double*, double*);


  /** NewtonEulerDS plug-in to compute \f$\nabla_qF_{Int}(\dot q, q, t)\f$, id = "jacobianFIntq"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr6 computeJacobianFIntqPtr;

  /** NewtonEulerDS plug-in to compute \f$\nabla_{\dot q}F_{Int}(\dot q, q, t)\f$, id = "jacobianFIntqDot"
   * @param time : current time
   * @param sizeOfq : size of vector q
   * @param q : pointer to the first element of q
   * @param velocity : pointer to the first element of velocity
   * @param[in,out] jacob : pointer to the first element of the jacobian
   * @param  size of vector z
   * @param[in,out] z  : a vector of user-defined parameters
   */
  FPtr6 computeJacobianFIntqDotPtr;

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

  /** allocate memory for forces and its jacobians, if required.
   */
  void initForces();

  /** Initialization function for the rhs and its jacobian.
   *  \param time of initialization
   */
  void initRhs(double) ;

  /** dynamical system initialization function except for _p:
   *  mainly set memory and compute plug-in for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  void initialize(double = 0, unsigned int = 1) ;

  /** dynamical system initialization function for _p
   *  \param level for _p
   */
  void initializeNonSmoothInput(unsigned int level) ;

  // === GETTERS AND SETTERS ===

  /** return the dim. of the system (n for first order). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline unsigned int getDim() const
  {
    return _n;
  }
  virtual inline unsigned int getqDim() const
  {
    return _qDim;
  }

  // -- q --

  /** get q
   *  \return pointer on a SiconosVector
   */
  inline SP::SimpleVector q() const
  {
    return _q;
  }
  inline SP::SimpleVector deltaq() const
  {
    return _deltaq;
  }


  // -- q0 --

  /** get q0
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector q0() const
  {
    return _q0;
  }
  inline SP::SiconosVector v0() const
  {
    return _v0;
  }

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



  // -- Mass --

  /** get mass value
   *  \return a double
   */
  inline double massValue() const
  {
    return _mass;
  };


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



  // -- forces --

  /** get the value of fL
   *  \return SimpleVector
   */
  inline const SimpleVector getForces() const
  {
    return *_forces;
  }

  /** get forces
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector forces() const
  {
    return _forces;
  }

  // -- Jacobian fL --


  /** get JacobianFL
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SimpleMatrix jacobianvFL() const
  {
    return _jacobianvFL;
  }
  /** get JacobianFL
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SimpleMatrix jacobianqDotForces() const
  {
    return _jacobianqDotForces;
  }
  //  inline SP::SiconosMatrix jacobianZFL() const { return jacobianZFL; }

  // --- PLUGINS RELATED FUNCTIONS ---



  /** allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeFExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    Plugin::setFunction(&computeFExtPtr, pluginPath, functionName);
  }
  /** allow to set a specified function to compute Mext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   */
  void setComputeMExtFunction(const std::string&  pluginPath, const std::string& functionName)
  {
    Plugin::setFunction(&computeMExtPtr, pluginPath, functionName);
  }

  /** set a specified function to compute fExt
   *  \param a pointer on the plugin function
   */
  void setComputeFExtFunction(Fext fct)
  {
    computeFExtPtr = fct ;
  }
  void setComputeMExtFunction(Fext fct)
  {
    computeMExtPtr = fct ;
  }



  /** default function to compute the external strengths
   *  \param double time : the current time
   */
  virtual void computeFExt(double);
  virtual void computeMExt(double);



  /** Default function to compute the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeRhs(double, bool  = false);

  /** Default function to compute jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeJacobianRhsx(double, bool  = false);

  /** Default function to compute fL
   *  \param double, the current time
   */
  virtual void computeForces(double);

  /** function to compute fL with some specific values for q and velocity (ie not those of the current state).
   *  \param double time : the current time
   *  \param SP::SiconosVector: pointers on q
   *  \param SP::SiconosVector: pointers on velocity
   */
  virtual void computeForces(double , SP::SiconosVector, SP::SiconosVector);

  /** Default function to compute the jacobian following q of fL
   *  \param double, the current time
   */
  virtual void computeJacobianvFL(double);
  /** Default function to compute the jacobian following qDot of fL
   *  \param double, the current time
   */
  virtual void computeJacobianqDotForces(double);

  // --- miscellaneous ---

  /** copy the data of the DS into the XML tree
   */
  void saveSpecificDataToXML();

  /** print the data to the screen
   */
  void display() const;

  /** initialize the SiconosMemory objects with a positive size.
   *  \param int levelMin for allocation of _p
   *  \param int levelMax for allocation of _p
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


  /** set p[...] to zero
   */
  void resetNonSmoothPart();


  virtual void updateT();
  virtual void normalizeq();

  inline SP::SimpleMatrix massMatrix()
  {
    return _massMatrix;
  }
  inline SP::SimpleMatrix luW()
  {
    return _luW;
  }
  inline SP::SimpleMatrix T()
  {
    return _T;
  }
  inline SP::SiconosMemory fLMemory()
  {
    return _forcesMemory;
  }
  inline SP::SiconosMemory dotqMemory()
  {
    return _dotqMemory;
  }
  inline SP::SiconosVector dotq()
  {
    return _dotq;
  }
  /*get the matrix converting the object coordinates in the absolute coordinates.*/
  SP::SimpleMatrix MObjToAbs()
  {
    return _MObjToAbs;
  }
  /*update the _MObjToAbs from the current quaternion.*/
  void updateMObjToAbs();

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(NewtonEulerDS);

#endif // NEWTONEULERNLDS_H

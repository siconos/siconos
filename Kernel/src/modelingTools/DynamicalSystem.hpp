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

/*! \file DynamicalSystem.hpp
  \brief Abstract class - General interface for all Dynamical Systems.
*/

#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosPointers.hpp"
#include "SSLH.hpp"
#include "RuntimeException.hpp"
#include "Tools.hpp"
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#include "SiconosMemory.hpp"
#include "DynamicalSystemTypes.hpp"
#include "PluggedObject.hpp"
#include "PluginTypes.hpp"
#include "SiconosVisitor.hpp"


class NonSmoothDynamicalSystem;
class DynamicalSystemXML;
class SiconosVector;
class SiconosMatrix;
class SimpleMatrix;
class SiconosVector;
class SiconosMemory;
class SiconosSharedLibrary;


/** Pointer to function for plug-in. */
typedef void (*FPtr6)(double, unsigned int, const double*, const double*, double*, unsigned int, double*);



/** */

/**  Abstract class to handle Dynamical Systems => interface for
   derived classes (First Order or Lagrangian systems)

   \author SICONOS Development Team - copyright INRIA
   \version 3.0.0.
   \date (Creation) January 15, 2007

  This class is used to describe dynamical systems of the form :
  \f[
  g(\dot x, x, t, z) = 0
  \f]
  where

     - \f$ x \in R^{n} \f$ is the state.

     - \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic
   variables, some sort of discret state.  For example, z may be used
   to set some perturbation parameters, or to control the system (z
   will be set by some actuators) or anything else.

   with \f$ g : R^{n} \times R  \to  R^{n}   \f$ .

  Operators and the functions used to compute them:

  - g: computeg(t)

  - jacobianXG[0] = \f$ \nabla_x g(t,\dot x,x,z) \f$:
    computeJacobiang(0,...)

  - jacobianXG[1] = \f$ \nabla_{\dot x} g(t,\dot x,x,z) \f$:
    computeJacobiang(1,...)

  By default, the DynamicalSystem is considered to be an Initial Value
  Problem (IVP) and the initial conditions are given by

   \f[
   x(t_0)=x_0
  \f]

  Under some specific conditions, the system can be written as:

  \f[
  \dot x = rhs(x, t, z)
  \f]

  In that case, \f$ \nabla_{\dot x} g \f$ must be invertible.
  Right-hand side (\f$ \dot x \f$) of the equation is computed thanks
  to computeRhs(t)

  and its Jacobian according to x, named jacobianRhsx, with
  computeJacobianRhsx(t).

  <b> Those two functions (computeRhs and computeJacobianRhsx) are
  pure virtual and must be implemented in all the derived
  classes. </b>

  Dynamical System types (followed by derived classes names):

   - First Order Non Linear Dynamical Systems (FirstOrderNonLinearDS)

   - First Order Linear DS (FirstOrderLinearDS)

   - First Order Linear and Time Invariant Coefficient DS
     (FirstOrderLinearTIDS)

   - Lagrangian DS (LagrangianDS)

   - Lagrangian Linear and Time Invariant coefficients DS
     (LagrangianLinearTIDS)

  About members:

  - A DynamicalSystem is identified thanks to a number.

   - A VectorOfVectors, x, is used to saved the state: x[0]=\f$ x \f$
     and x[1]=\f$ \dot x \f$ = right-hand side.

   - number is set automatically using count static variable except in
     the case of XML loading, where number is read in the xml file,
     because it's to be given explicitely by user to set the list of
     DS in the Interactions.

   Warning:

   - At the time, nothing is implemented in simulation to proceed with
     systems written as \f$ g(...) = 0 \f$. Then use only the form \f$
     \dot x = rhs(...) \f$.

 */
class DynamicalSystem
{

public:
  /** List of indices used to save tmp work vectors
   * The last value is the size of the present list, so you HAVE to leave it at the end position.
   */
  enum WorkNames {local_buffer, qtmp, sizeWorkV};

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(DynamicalSystem);

  /** used to set ds number */
  static unsigned int count;

  /** copy constructor => private, no copy nor pass-by-value.
   */
  //  DynamicalSystem(const DynamicalSystem & );

protected:

  /** An id number for the DynamicalSystem */
  unsigned int _number;

  /** the dimension of the system (\e ie size of the state vector x) */
  unsigned int _n;

  /** initial state of the system */
  SP::SiconosVector _x0;

  /** ResiduFree  */
  SP::SiconosVector _residuFree;

  /** the input vector due to the non-smooth law \f$ r \in R^{n}\f$
   * (multiplier, force, ...)
   * \remark V.A. 17/09/2011 :
   * This should be a VectorOfVectors as for _x when higher relative degree
   * systems will be simulated
   */
  SP::SiconosVector _r;

  /** used by the relative convergence criteron*/
  double _normRef;

  /** state of the system, \f$  x \in R^{n}\f$ - With \f$ x[0]=\f$ x \f$ , x[1]= \f$ \dot x \f$ . */
  VectorOfVectors _x;

  /** jacobian according to x of the right-hand side (\f$ \dot x =
      f(x,t) + r \f$) */
  SP::SiconosMatrix _jacxRhs;

  //  SP::SiconosMatrix _jacgx;
  //  SP::SiconosMatrix _jacxDotG;
  //  SP::SiconosMatrix jacobianZG;

  /** Arbitrary algebraic values vector, z, discret state of the
      system. */
  SP::SiconosVector _z;
  SP::SiconosVector _g;





  /** DynamicalSystem plug-in to compute \f$ g(t,\dot x,x,z) \f$
   *  @param   current time
   *  @param   the size of the vector x
   *  @param   the pointer to the first element of the vector x[0]=\f$ x \f$
   *  @param   the pointer to the first element of the vector x[1]=\f$ \dot x \f$
   *  @param   the pointer to the first element of the vector g(t, ...)
   *  @param   the size of the vector z
   *  @param   a vector of parameters, z
   */
  SP::PluggedObject _pluging;

  /** Plug-in to compute jacobianG (computeJacobiangPtr[i] for jacobianG[i]).
   *  @param   current time
   *  @param   the size of the vector x
   *  @param   the pointer to the first element of the vector x[0]=\f$ x \f$
   *  @param   the pointer to the first element of the vector x[1]=\f$ \dot x \f$
   *  @param   the pointer to the first element of the vector g(t, ...)
   *  @param   the size of the vector z
   *  @param   a vector of parameters, z
   */
  SP::PluggedObject _pluginJacgx;
  SP::PluggedObject _pluginJacxDotG;


  /** the  previous state vectors stored in memory*/
  SP::SiconosMemory _xMemory;

  /** number of previous states stored in memory */
  unsigned int _stepsInMemory;

  /** A container of vectors to save temporary values (for Newton convergence computation for example)*/
  VectorOfVectors _workV;

  /** A container of matrices to save temporary values (zeroMatrix, idMatrix, inverse of Mass or any tmp work matrix ...)
   * No get-set functions at the time. Only used as a protected member.*/
  VectorOfSimpleMatrices _workMatrix;

  /** the XML object linked to the DynamicalSystem  */
  SP::DynamicalSystemXML _dsxml;

  // ===== CONSTRUCTORS =====

  /** default constructors/destructor
   */
  DynamicalSystem();
  virtual void zeroPlugin();

protected:
  /** a vector reserved to compute the freeState.*/
  SP::SiconosVector _workFree;

public:

  /*! @name Constructors */
  //@{

  /** xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   */
  DynamicalSystem(SP::DynamicalSystemXML dsXML);

  /** constructor from a set of data
   *  \param type of the system
   *  \param int : size of the system (n)
   */
  DynamicalSystem(unsigned int newN);

  /** Copy constructor
   * \param ds the DynamicalSystem to copy
   */
  DynamicalSystem(const DynamicalSystem & ds);
  // ===== DESTRUCTOR =====

  /** destructor
   */
  virtual ~DynamicalSystem() {}

  //@}

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem() = 0;

  /*! @name Members access */
  //@{

  // --- Number ---

  /** to get the number of the DynamicalSystem
   *  \return the value of number
   */
  inline int number() const
  {
    return _number;
  }

  /** function used to sort DynamicalSystem in SiconosSet<SP::DynamicalSystem>
   *  \return an int (warning: must be const, despite intel compilers warning, because of SiconosSet Cmp function arguments)
   */
  inline int getSort() const
  {
    return _number;
  }

  // --- n ---

  /** allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline unsigned int getN() const
  {
    return _n;
  }

  /** allows to set the value of n
   *  \param an integer to set the value of n
   */
  inline void setN(unsigned int newN)
  {
    _n = newN;
  }

  /** return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline unsigned int getDim() const
  {
    return _n;
  };

  // --- X0 ---

  /** get the value of x0, the initial state of the DynamicalSystem
   *  \return SiconosVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */
  inline const SiconosVector getX0() const
  {
    return *_x0;
  }

  /** get x0, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector x0() const
  {
    return _x0;
  };

  /** get _normRef
   * \return a reference to _normRef
   */
  inline const double& normRef() const
  {
    return _normRef;
  };

  // --- R ---

  /** get the value of r
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   *  \return a vector
   */
  inline const SiconosVector getR() const
  {
    return *_r;
  }

  /** get r
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector r() const
  {
    return _r;
  }

  /** set the value of r to newValue
   *  \param SiconosVector newValue
   */
  void setR(const SiconosVector&);

  /** set R to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setRPtr(SP::SiconosVector);

  // --- Residu ---

  /** get Residu,
    *  \return pointer on a SiconosVector
    */
  inline SP::SiconosVector residuFree() const
  {
    return _residuFree;
  }

  /** set the value of x0 to newValue
   *  \param SiconosVector newValue
   */
  void setX0(const SiconosVector&);

  /** set x0 to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setX0Ptr(SP::SiconosVector);

  // --- X ---

  /** get the value of \f$ x \f$, the state of the DynamicalSystem
   *  \return SiconosVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */

  inline const SiconosVector getx() const
  {
    return *(_x[0]);
  }

  /** get \f$ x \f$ (pointer), the state of the DynamicalSystem.
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector x() const
  {
    return _x[0];
  }

  /** set the value of \f$ x \f$ (ie (*x)[0]) to newValue
   *  \param SiconosVector newValue
   */
  void setX(const SiconosVector&);

  /** set \f$ x \f$ (ie (*x)[0]) to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setXPtr(SP::SiconosVector);

  // ---  Rhs ---

  /** get the value of the right-hand side, \f$ \dot x \f$, derivative of the state of the DynamicalSystem.
   *  \return SiconosVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */
  inline const SiconosVector getRhs() const
  {
    return *(_x[1]);
  }

  /** get the right-hand side, \f$ \dot x \f$, the derivative of the state of the DynamicalSystem.
   *  \return a pointer on a SiconosVector
   */
  inline SP::SiconosVector rhs() const
  {
    return _x[1];
  }

  /** set the value of the right-hand side, \f$ \dot x \f$, to newValue
   *  \param SiconosVector newValue
   */
  void setRhs(const SiconosVector&);

  /** set right-hand side, \f$ \dot x \f$, to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setRhsPtr(SP::SiconosVector);

  // --- JacobianRhsx ---

  /** get the value of the gradient according to \f$ x \f$ of the right-hand side
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianRhsx() const
  {
    return *_jacxRhs;
  }

  /** get gradient according to \f$ x \f$ of the right-hand side (pointer)
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix jacobianRhsx() const
  {
    return _jacxRhs;
  }

  /** set the value of JacobianRhsx to newValue
   *  \param SiconosMatrix newValue
   */
  void setJacobianRhsx(const SiconosMatrix&);

  /** set JacobianRhsx to pointer newPtr
   *  \param SP::SiconosMatrix  newPtr
   */
  void setJacobianRhsxPtr(SP::SiconosMatrix newPtr);

  // -- z --

  /** get the value of \f$ z \f$, the vector of algebraic parameters.
   * \return a SiconosVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SiconosVector
   */
  inline const SiconosVector getz() const
  {
    return *_z;
  }

  /** get \f$ z \f$ (pointer), the vector of algebraic parameters.
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector z() const
  {
    return _z;
  }

  /** set the value of \f$ z \f$ to newValue
   *  \param SiconosVector newValue
   */
  void setz(const SiconosVector&);

  /** set \f$ z \f$ to pointer newPtr
   *  \param SP::SiconosVector newPtr
   */
  void setzPtr(SP::SiconosVector);

  // --- g ---
  /** get the value of g
   *  \return plugged vector
  inline const PVFInt getg() const { return *g; }
   */

  /** get g
   *  \return pointer on a plugged vector
   */
  inline SP::SiconosVector g() const
  {
    return _g;
  }

  /** set the value of g to newValue
   *  \param a plugged vector

  void setG(const PVFInt&);
  */

  /** set g to pointer newPtr
   *  \param a SP to plugged vector
   */
  inline void setGPtr(SP::SiconosVector newPtr)
  {
    _g = newPtr;
  }

  // -- Jacobian g --
  /** get the value of jacobianG
    \param index of the desired jacobian
    *  \return a plugged-matrix

  inline const PMFInt getJacobianG(unsigned int i) const { return *(jacobianG[i]); }
  */

  /** get jacobianG
      \param index of the desired jacobian
      *  \return pointer on a plugged-matrix
      */
  //  inline SP::SiconosMatrix jacobianXG() const { return _jacgx; }
  //  inline SP::SiconosMatrix jacobianXDotG() const { return _jacxDotG; }
  //  inline SP::SiconosMatrix jacobianZG() const { return jacobianZG; }

  /** set the value of jacobianG to newValue
      \param index of the desired jacobian
      *  \param plugged-matrix newValue
  void setJacobianG(unsigned int, const PMFInt&);
      */

  /** set jacobianG to pointer newPtr
      \param index of the desired jacobian
      *  \param a plugged matrix SP
      */
  //  inline void setJacobianXGPtr( SP::SiconosMatrix newPtr) {_jacgx = newPtr;}
  //  inline void setJacobianXDotGPtr( SP::SiconosMatrix newPtr) {_jacxDotG = newPtr;}
  //  inline void setJacobianZGPtr( SP::SiconosMatrix newPtr) {jacobianZG = newPtr;}

  // X memory

  /** get all the values of the state vector x stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory xMemory() const
  {
    return _xMemory;
  }

  // --- Steps in memory ---

  /** get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline int getStepsInMemory() const
  {
    return _stepsInMemory;
  }

  /** set the value of stepsInMemory
   *  \param int steps : the value to set stepsInMemory
   */
  inline void setStepsInMemory(int steps)
  {
    _stepsInMemory = steps;
  }

  // --- dsxml ---

  /** get the object DynamicalSystemXML of the DynamicalSystem
   *  \return a pointer on the DynamicalSystemXML of the DynamicalSystem
   */
  inline const SP::DynamicalSystemXML dynamicalSystemXML() const
  {
    return _dsxml;
  }

  /** set the DynamicalSystemXML of the DynamicalSystem
   *  \param DynamicalSystemXML* dsxml : the address of theDynamicalSystemXML to set
   */
  inline void setDynamicalSystemXMLPtr(SP::DynamicalSystemXML newDsxml)
  {
    _dsxml = newDsxml;
  }

  // ===== WORK VECTOR =====

  /** get the vector of temporary saved vector
   *  \return a VectorOfVectors (map that links string to vectors)
   */
  inline VectorOfVectors getWorkVector() const
  {
    return _workV;
  }

  /** get a temporary saved vector, ref by id
   *  \return a SP::SiconosVector
   */
  inline SP::SiconosVector getWorkVector(const WorkNames& id) const
  {
    return _workV[id];
  }

  /** set WorkVector
   *  \param a VectorOfVectors
   */
  inline void setWorkVector(const VectorOfVectors& newVect)
  {
    _workV = newVect;
  }

  /** to add a temporary vector
   *  \param a SP::SiconosVector
   *  \param a string id
   */
  inline void addWorkVector(SP::SiconosVector newVal, const WorkNames& id)
  {
    *_workV[id] = *newVal;
  }
  /** sub a vector to a temporary one
   *  \param a SP::SiconosVector
   *  \param a string id
   */
  inline void subWorkVector(SP::SiconosVector newVal, const WorkNames& id)
  {
    *_workV[id] -= *newVal;
  }

  /** to allocate memory for a new vector in tmp map
   *  \param the id of the SiconosVector
   *  \param an int to set the size
   */
  inline void allocateWorkVector(const WorkNames& id, int size)
  {
    _workV[id].reset(new SiconosVector(size));
  }

  /** to get the free vector.
   */
  inline SP::SiconosVector workFree() const
  {
    return _workFree;
  };

  //@}

  /**
   * return true if the Dynamical system is linear.
   */
  virtual bool isLinear()
  {
    return false;
  }

  /** Initialization function for the rhs and its jacobian (including
   *  memory allocation).  \param time of initialization
   */
  virtual void initRhs(double) = 0 ;

  /** dynamical system initialization function except for _r :
   *  mainly set memory and compute value for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(double = 0, unsigned int = 1) = 0;

  /** dynamical system initialization function for NonSmoothInput _r
   *  \param level of _r.
   */
  virtual void initializeNonSmoothInput(unsigned int level) = 0;

  /** dynamical system update: mainly call compute for all time or state depending functions
   *  \param current time
   */
  void update(double);

  /*! @name Memory vectors management  */
  //@{
  /** initialize the SiconosMemory objects: reserve memory for i vectors in memory and reset all to zero.
   *  \param the size of the SiconosMemory (i)
   */
  virtual void initMemory(unsigned int);

  /** push the current values of x and r in memory (index 0 of memory is the last inserted vector)
   *  xMemory and rMemory,
   */
  virtual void swapInMemory() = 0;

  //@}
  /*! @name Plugins management  */
  //@{
  /** to set a specified function to compute g
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   */
  void setComputegFunction(const std::string&  pluginPath, const std::string& functionName);

  /** set a specified function to compute g
   *  \param a pointer on the plugin function
   */
  void setComputegFunction(FPtr6 fct);

  /** to set a specified function to compute jacobianG
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   */
  void setComputeJacobianXGFunction(const std::string&  pluginPath, const std::string&  functionName);
  void setComputeJacobianDotXGFunction(const std::string&  pluginPath, const std::string&  functionName);
  //  void setComputeJacobianZGFunction( const std::string&  pluginPath, const std::string&  functionName);

  /** set a specified function to compute jacobianG
   *  \param a pointer on the plugin function
   */
  virtual void setComputeJacobianXGFunction(FPtr6) {};
  virtual void setComputeJacobianDotXGFunction(FPtr6) {};
  //  void setComputeJacobianZGFunction( FPtr6);

  /** Default function to compute g
   *  \param double, the current time
   */
  void computeg(double);

  /** default function to compute the gradient of g
   *  \param double time : the current time
   */
  //  void computeJacobianXG(double);
  //  void computeJacobianDotXG(double);

  /**
   *default function to update the plugins functions using a new time:
   * \param double time : the current time
   */
  virtual void updatePlugins(double)
  {
    ;
  }

  //@}


  /*! @name Right-hand side computation */
  //@{
  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeRhs(double, bool  = false) = 0;

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   */
  virtual void computeJacobianRhsx(double, bool  = false) = 0;

  //@}

  /*! @name XML */
  //@{
  /** copy the data of the DS into the XML tree
   */
  virtual void saveDSToXML();

  /** copy the data specific to each system into the XML tree
   */
  virtual void saveSpecificDataToXML() = 0;
  //@}

  // ===== MISCELLANEOUS ====

  /** print the data of the dynamical system on the standard output
   */
  virtual void display() const = 0;

  /** Default function for computing an indicator of convergence
   *  \return a double when DS is a Lagrangian
   */
  virtual double dsConvergenceIndicator() ;

  /** set R to zero
   */
  virtual void resetAllNonSmoothPart() = 0;

  /** set R to zero for a given level
   */
  virtual void resetNonSmoothPart(unsigned int level) = 0;

  /**
   * overwrite these methods to do the specific work that must be done
   * at the beginning or the end of a computation step.  It could be
   * reset some buffer vector.(like the function resetNonSmoothPart).
   *
   */
  virtual void preparStep() {};

  virtual void endStep() {};

  /** Get _pluging
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginG() const
  {
    return _pluging;
  };

  /** Get _pluginJacgx
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginJacGX() const
  {
    return _pluginJacgx;
  };

  /** Get _pluginJacxDotG
   * \return a SP::PluggedObject
   */
  inline SP::PluggedObject getPluginJacXDotG() const
  {
    return _pluginJacxDotG;
  };

  /** visitors hook
   */
  VIRTUAL_ACCEPT_VISITORS(DynamicalSystem);

};

TYPEDEF_SPTR(DynamicalSystem)

#endif // DYNAMICALSYSTEM_H

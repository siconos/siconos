/* Siconos-Kernel version 1.2.0, Copyright INRIA 2005-2006.
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
#ifndef DYNAMICALSYSTEM_H
#define DYNAMICALSYSTEM_H

#include "SiconosConst.h"
#include "RuntimeException.h"
#include "check.h"

#include "SimpleMatrix.h"
#include "BlockMatrix.h"
#include "SimpleVector.h"
#include "SiconosMemory.h"
#include "SiconosSharedLibrary.h"

#include "NonSmoothDynamicalSystem.h"
#include "DSInputOutput.h"
#include "BoundaryCondition.h"
#include "DynamicalSystemXML.h"

#include <string>
#include <vector>
#include <iostream>
#include <map>

const std::string LNLDS = "LagrangianDS";
const std::string LTIDS = "LagrangianLinearTIDS";
const std::string LDS = "LinearDS";
const std::string LITIDS = "LinearTIDS";
const std::string NLDS = "NonLinearDS";

class NonSmoothDynamicalSystem;
class BoundaryCondition;
class DSInputOutput;
class DynamicalSystemXML;
class SiconosVector;
class SimpleMatrix;
class BlockMatrix;
class SiconosMemory;
class SiconosSharedLibrary;

/** \class DynamicalSystem
 *  \brief  General first order non linear dynamical systems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.2.0.
 *  \date (Creation) April 29, 2004
 *
 *
 * This class defines and computes a generic n-dimensional
 * dynamical system of the form :
 * \f[
 * \dot x = f(x,t) + T(x) u(x,t) + r,
 * \f]
 * where
 *    - \f$x \in R^{n} \f$ is the state,
 *    - \f$ r \in R^{n} \f$  the input due to the Non Smooth Interaction.
 *    - \f$ u \in R^{uSize}\f$ a "control" term
 *
 *  with \f$ f : R^{n} \times R  \mapsto  R^{n}   \f$ .
 *
 *
 * By default, the DynamicalSystem is considered to be an Initial Value Problem (IVP)
 * and the initial conditions are given by
 *  * \f[
 *  x(t_0)=x_0
 * \f]
 * To define a boundary Value Problem, the pointer on  a BoundaryCondition must be set.
 *
 * \f$ f(x,t) \f$ is a plug-in function, and can be computed using computeF(t).
 * Its Jacobian according to x is denoted jacobianXF, and computed thanks to computeJacobianXF(t).
 * f and jacobianXF can be plugged to external functions thanks to setComputeFFunction/setComputeJacobianXFFunction.
 *
 * Right-hand side of the equation is denoted rhs and is computed thanks to computeRhs(t).
 *
 * \f[ rhs(x,t) =  f(x,t) + T(x) u(x,t)\f]
 *
 * Its Jacobian according to x is jacobianXRhs:
 *
 *  \f[
 *   jacobianXRhs = \nabla_xrhs(x,t) = \nabla_xf(x,t) + \nabla_xT(x)U(x,t)
 *  \f]
 *
 * At the time, \f$\nabla_xT(x)\f$ is never taken into account.
 * Moreover, in order to avoid useless memory allocation,
 * we try to re-use the same memory location for f and rhs (and for jacobianXRhs and jacobianXF).
 *
 * \todo Add a pointer to an object Constraint .
 */
class DynamicalSystem
{
protected:

  /** Dynamical System type: General Dynamical System (NLDS) LagrangianDS (LNLDS),
      LagrangianLinearTIDS (LTIDS), LinearDS (LDS)*/
  std::string  DSType;

  /** NonSmoothDynamicalSystem owner of this DynamicalSystem */
  NonSmoothDynamicalSystem* nsds;

  /** this number defines in a single way the DynamicalSystem */
  int number;

  /** the name of the DS ("ball", "solid1254", etc.)*/
  std::string  id;

  /** the dimension of the system (i.e. size of the state vector x)*/
  unsigned int n;

  /** initial state of the system */
  SiconosVector *x0;

  /** state of the system, \f$  x \in R^{n}\f$ */
  SiconosVector *x;

  /** the  previous state vectors stored in memory*/
  SiconosMemory *xMemory;

  /** the right hand side of the equation */
  SiconosVector *rhs;

  /** jacobian according to X of this rhs */
  SiconosMatrix * jacobianXRhs;

  /** the  free state vector (state vector for r=0) */
  SiconosVector *xFree;

  /** the  input vector due to the non-smooth law \f$  r \in R^{n}\f$ (multiplier, force, ...)*/
  SiconosVector *r;

  /**  the previous r vectors */
  SiconosMemory *rMemory;

  /** f(x,t) */
  SiconosVector *f;

  /** Gradient of \f$ f(x,t) \f$ with respect to \f$ x\f$*/
  SiconosMatrix *jacobianXF;

  /** size of vector u */
  unsigned int uSize;

  /** "control" term */
  SiconosVector *u;

  /** Matrix coefficient of u */
  SiconosMatrix *T;

  /** number of previous states stored in memory */
  unsigned int stepsInMemory;

  /** A container of vectors to save temporary values (for Newton convergence computation for example)*/
  std::map<const std::string, SimpleVector*> tmpWorkVector;

  /** boundary conditions defined if the DynamicalSystem has some */
  BoundaryCondition *BC;

  /** the XML object linked to the DynamicalSystem  */
  DynamicalSystemXML *dsxml;

  /** vector of the DS Inputs-Outputs of the Dynamical System */
  std::vector<DSInputOutput*> dsioVector;

  /** Parameters list, last argument of plug-in functions. What are those parameters depends on user´s choice.
   *  This a list of pointer to SimpleVector. Each one is identified thanks to a key which is the plug-in name.
   * A flag is also added in the isAllocatedIn map to check inside-class memory allocation for this object.*/
  std::map<std::string, SimpleVector*> parametersList;

  // --- plugins ---

  /** class for plugin managing (open, close librairy...) */
  SiconosSharedLibrary cShared;

  /* the name of the plugin used to compute f(x,t) */
  std::string  computeFFunctionName;

  /* the name of the plugin used to compute the jacobian of f according to x */
  std::string  computeJacobianXFFunctionName;

  /* the name of the plugin used to compute u */
  std::string  computeUFunctionName;

  /* the name of the plugin used to compute T */
  std::string  computeTFunctionName;

  /** Flag to check if operators are plugged or not (and thus constant)
   * For example isPlugin["jacobianXF"] = false, means that jacobianXF is constant,
   * then computeJacobianXF does not change its value, and not plugged.*/
  std::map<std::string, bool> isPlugin;

  /** \fn void (*computeFPtr) (const unsigned int& sizeOfX, const double* t, const double* x, double* f, double* param)
   *  \brief pointer on function to compute f(x,t)
   *  \param unsigned int& sizeOfX : the size of the vector x
   *  \param double* time : current time
   *  \param double* x : the pointer to the first element of the vector x
   *  \param double* f : the pointer to the first element of the vector f(x,t)
   *    \param double* param   : a vector of user-defined parameters
   */
  void (*computeFPtr)(const unsigned int&, const double*, const double*, double*, double*);

  /** \fn void (*computeJacobianXFPtr)(const unsigned int& sizeOfX, const double* t, const double* x, double* jacobianXF, double* param)
   *  \brief  Pointer on function to compute the gradient of f(x,t) with respect to the state  \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param unsigned int& sizeOfX : size of vector x
   *  \param double* time : current time
   *  \param double* xPtr : pointer to the first element of x
   *  \param double* jacobianXFPtr : pointer to the first element of jacobianXF matrix (in-out parameter)
   *    \param double* param   : a vector of user-defined parameters
   */
  void (*computeJacobianXFPtr)(const unsigned int&, const double*, const double*, double*, double*);

  /** \fn void (*computeUPtr)(const unsigned int& sizeOfU, const unsigned int& sizeOfX, const double* t, const double* x, double* u, double* param)
   *  \brief  Pointer on function to compute u
   *  \param int* sizeOfU : size of vector u
   *  \param int* sizeOfX : size of vector x
   *  \param double* time : current time
   *  \param double* xPtr : pointer to the first element of x
   *  \param double* UPtr : pointer to the first element of u vector (in-out parameter)
   *    \param double* param   : a vector of user-defined parameters
   */
  void (*computeUPtr)(const unsigned int&, const unsigned int&, const double*, const double*, double*, double*);

  /** \fn void (*computeTPtr)(const unsigned int& sizeOfU, const unsigned int& sizeOfX, const double* x, double* T, double* param)
   *  \brief  Pointer on function to compute T
   *  \param unsigned int& sizeOfU : size of vector u
   *  \param unsigned int& sizeOfX : size of vector X
   *  \param double* x : pointer to the first element of X
   *  \param double* T: pointer to the first element of T matrix (in-out parameter)
   *    \param double* param   : a vector of user-defined parameters
   */
  void (*computeTPtr)(const unsigned int&, const unsigned int&, const double*, double*, double*);

  /** Flags to know if pointers have been allocated inside constructors or not */

  /** Flags to know if pointers have been allocated inside constructors or not */
  std::map<std::string, bool> isAllocatedIn;
  /** dsio */
  std::deque<bool> isDsioAllocatedIn;

  /** \fn void fillBoundaryConditionsFromXml()
   *  \brief uses the DynamicalSystemXML of the DynamicalSystem to fill BoundaryCondition fields
   *  \exception RuntimeException
   */
  virtual void fillBoundaryConditionsFromXml();

  /** \fn void fillDsioFromXml()
   *  \brief uses the DynamicalSystemXML of the DynamicalSystem to fill DSIO vector
   *  \exception RuntimeException
   */
  virtual void fillDsioFromXml();

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

  /** \fn initParameter(const string& id);
   *  \brief init parameter vector corresponding to id to a SimpleVector* of size 1
   *  \param a string, id of the plug-in
   */
  void initParameter(const std::string&);

  /** \fn DynamicalSystem();
   *  \brief default constructor
   */
  DynamicalSystem();

public:

  // ===== CONSTRUCTORS =====

  /** \fn DynamicalSystem(DynamicalSystemXML * nsdsXML, NonSmoothDynamicalSystem* =NULL)
   *  \brief xml constructor
   *  \param DynamicalSystemXML* : the XML object for this DynamicalSystem
   *  \param NonSmoothDynamicalSystem* (optional): the NSDS that owns this ds
   *  \exception RuntimeException
   */
  DynamicalSystem(DynamicalSystemXML * dsXML, NonSmoothDynamicalSystem* = NULL);

  /** \fn DynamicalSystem(DynamicalSystemXML * nsdsXML, const unsigned int& number, const unsigned int& n,
   * const SiconosVector& x0, const string& fPlugin = "DefaultPlugin:f, const string& jacobianXFPlugin = "DefaultPlugin:jacobianXF)
   *  \brief constructor from a set of data
   *  \param int : reference number for this DynamicalSystem
   *  \param int : dimension of this DynamicalSystem
   *  \param SiconosVector : initial state of this DynamicalSystem
   *  \param string : plugin name for f of this DynamicalSystem (optional)
   *  \param string : plugin name for jacobianXF of this DynamicalSystem (optional)
   *  \exception RuntimeException
   */
  DynamicalSystem(const int&, const unsigned int&, const SiconosVector&,
                  const std::string& = "DefaultPlugin:computeF", const std::string& = "DefaultPlugin:computeJacobianXF");

  /** \fn DynamicalSystem(const DynamicalSystem &)
   *  \brief copy constructor
   *  \param a Dynamical system to copy
   */
  DynamicalSystem(const DynamicalSystem &);

  // ===== DESTRUCTOR =====

  /** \fn ~DynamicalSystem();
   *  \brief destructor
   */
  virtual ~DynamicalSystem();

  /** \fn bool checkDynamicalSystem()
   *  \brief check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  virtual bool checkDynamicalSystem();

  // ===== GETTERS/SETTERS =====

  // --- type of DS ---

  /** \fn inline string getType()
   *  \brief get the type of a DynamicalSystem
   *  \return string : the type of the DynamicalSystem
   */
  inline const std::string  getType() const
  {
    return DSType;
  }

  /** \fn inline string setType()
  *  \brief set the type of a DynamicalSystem
  *  \param string : the type of the DynamicalSystem
  */
  inline void setType(const std::string newType)
  {
    DSType = newType;
  }

  // --- NonSmoothDynamicalSystem ---

  /** \fn NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr(void) const;
   *  \brief get the NonSmoothDynamicalSystem containing this DynamicalSystem
   *  \return NonSmoothDynamicalSystem*
   */
  inline NonSmoothDynamicalSystem* getNonSmoothDynamicalSystemPtr() const
  {
    return nsds;
  }

  /** \fn void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem*);
   *  \brief set the NonSmoothDynamicalSystem containing the DynamicalSystem
   *  \param NonSmoothDynamicalSystem*
   */
  inline void setNonSmoothDynamicalSystemPtr(NonSmoothDynamicalSystem *newNsds)
  {
    nsds = newNsds;
  }

  // --- Number ---

  /** \fn const int getNumber(void) const;
   *  \brief to get the number of the DynamicalSystem
   *  \return the value of number
   */
  inline const int getNumber() const
  {
    return number;
  }

  /** \fn const int getNumberForSorting(void) const;
   *  \brief same as getNumber, but return an unsigned long int, used for set<DynamicalSystem*> in OSI, NSDS ...
   *   as sorting criterion.
   *  \return the value of number
   */
  inline const unsigned long int getNumberForSorting() const
  {
    return number;
  }

  /** \fn void setNumber(const int&)
   *  \brief allows to set the value of number
   *  \param an integer to set the value of number
   */
  inline void setNumber(const int& newNumber)
  {
    number = newNumber;
  }

  // --- Id ---

  /** \fn const string getId(void) const
   *  \brief allows to get the id of the DynamicalSystem
   *  \return the value of ths id
   */
  inline const std::string  getId() const
  {
    return id;
  }

  /** \fn void setId(const string&)
   *  \brief allows to set the value of id
   *  \param a string to set the value of id
   */
  inline void setId(const std::string & newId)
  {
    id = newId;
  }

  // --- n ---

  /** \fn const unsigned int getN(void) const;
   *  \brief allow to get n, the dimension, i.e. the size of the state x of the DynamicalSystem
   *  \return the value of n
   */
  inline const unsigned int getN(void) const
  {
    return n;
  }

  /** \fn void setN(const unsigned int&)
   *  \brief allows to set the value of n
   *  \param an integer to set the value of n
   */
  inline void setN(const unsigned int& newN)
  {
    n = newN;
  }

  /** \fn const unsigned int getDim(void) const;
   *  \brief return the dim. of the system (n for first order, ndof for Lagrangian). Usefull to avoid if(typeOfDS) when size is required.
   *  \return an unsigned int.
   */
  virtual inline const unsigned int getDim(void) const
  {
    return n;
  }

  // --- X0 ---

  /** \fn  const SimpleVector getX0(void) const
   *  \brief get the value of x0, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getX0() const
  {
    return *x0;
  }

  /** \fn SiconosVector* getX0Ptr(void) const
   *  \brief get x0, the initial state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getX0Ptr() const
  {
    return x0;
  }

  /** \fn void setX0(const SiconosVector& newValue)
   *  \brief set the value of x0 to newValue
   *  \param SiconosVector newValue
   */
  void setX0(const SiconosVector&);

  /** \fn void setX0Ptr(SiconosVector* newPtr)
   *  \brief set x0 to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setX0Ptr(SiconosVector*);

  // --- X ---

  /** \fn const SimpleVector getX(void) const
   *  \brief get the value of x, the state of the DynamicalSystem
   *  \return SimpleVector
    * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getX() const
  {
    return *x;
  }

  /** \fn SiconosVector* getXPtr(void) const
   *  \brief get x, the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXPtr() const
  {
    return x;
  }

  /** \fn void setX (const SiconosVector& newValue)
   *  \brief set the value of x to newValue
   *  \param SiconosVector newValue
   */
  void setX(const SiconosVector&);

  /** \fn void setXPtr(SiconosVector* newPtr)
   *  \brief set x to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXPtr(SiconosVector *);

  // X memory

  /** \fn  const SiconosMemory getXMemory(void) const
   *  \brief get the value of xMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getXMemory() const
  {
    return *xMemory;
  }

  /** \fn SiconosMemory getXMemoryPtr(void) const
   *  \brief get all the values of the state vector x stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getXMemoryPtr() const
  {
    return xMemory;
  }

  /** \fn void setXMemory(const SiconosMemory &)
   *  \brief set the value of xMemory
   *  \param a ref on a SiconosMemory
   */
  void setXMemory(const SiconosMemory&);

  /** \fn void setXMemory(SiconosMemory * newPtr)
   *  \brief set xMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setXMemoryPtr(SiconosMemory *);

  // ---  Rhs ---

  /** \fn  const SimpleVector getRhs(void) const
   *  \brief get the value of rhs derivative of the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getRhs() const
  {
    return *rhs;
  }

  /** \fn SiconosVector* getRhsPtr(void) const
   *  \brief get rhs, the derivative of the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getRhsPtr() const
  {
    return rhs;
  }

  /** \fn void setRhs (const SiconosVector& newValue)
   *  \brief set the value of rhs to newValue
   *  \param SiconosVector newValue
   */
  void setRhs(const SiconosVector&);

  /** \fn void setRhsPtr(SiconosVector* newPtr)
   *  \brief set rhs to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setRhsPtr(SiconosVector *);

  // --- JacobianXRhs ---

  /** \fn  const SimpleMatrix getJacobianXRhs(void) const
   *  \brief get the value of JacobianXRhs
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXRhs() const
  {
    return *jacobianXRhs;
  }

  /** \fn SiconosMatrix* getJacobianXRhsPtr(void) const
   *  \brief get JacobianXRhs
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianXRhsPtr() const
  {
    return jacobianXRhs;
  }

  /** \fn void setJacobianXRhs (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianXRhs to newValue
   *  \param SiconosMatrix newValue
   */
  virtual void setJacobianXRhs(const SiconosMatrix&);

  /** \fn void setJacobianXRhsPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianXRhs to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  virtual void setJacobianXRhsPtr(SiconosMatrix *newPtr);

  // --- XFree ---

  /** \fn  const SimpleVector getXFree(void) const
   *  \brief get the value of xFree
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getXFree() const
  {
    return *xFree;
  }

  /** \fn SiconosVector* getXFreePtr(void) const
   *  \brief get xFree
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getXFreePtr() const
  {
    return xFree;
  }

  /** \fn void setXFree (const SiconosVector& newValue)
   *  \brief set the value of xFree to newValue
   *  \param SiconosVector newValue
   */
  void setXFree(const SiconosVector&);

  /** \fn void setXFreePtr(SiconosVector* newPtr)
   *  \brief set xFree to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setXFreePtr(SiconosVector *);

  // --- R ---

  /** \fn  const SiconosVector getR(void) const
   *  \brief get the value of r
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   *  \return a SiconosVector
   */
  inline const SimpleVector getR() const
  {
    return *r;
  }

  /** \fn SiconosVector* getRPtr(void) const
   *  \brief get r
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getRPtr() const
  {
    return r;
  }

  /** \fn void setR (const SiconosVector& newValue)
   *  \brief set the value of r to newValue
   *  \param SiconosVector newValue
   */
  void setR(const SiconosVector&);

  /** \fn void setRPtr(SiconosVector* newPtr)
   *  \brief set R to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setRPtr(SiconosVector *);

  // rMemory

  /** \fn  const SiconosMemory getRMemory(void) const
   *  \brief get the value of rMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getRMemory() const
  {
    return *rMemory;
  }

  /** \fn SiconosMemory getRMemoryPtr(void) const
   *  \brief get all the values of the state vector r stored in memory
   *  \return a memory
   */
  inline SiconosMemory* getRMemoryPtr() const
  {
    return rMemory;
  }

  /** \fn void setRMemory(const SiconosMemory &)
   *  \brief set the value of rMemory
   *  \param a ref on a SiconosMemory
   */
  void setRMemory(const SiconosMemory&);

  /** \fn void setRMemory(SiconosMemory * newPtr)
   *  \brief set rMemory to pointer newPtr
   *  \param a ref on a SiconosMemory
   */
  void setRMemoryPtr(SiconosMemory *);

  // ---  F ---

  /** \fn  const SimpleVector getF(void) const
   *  \brief get the value of f derivative of the state of the DynamicalSystem
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getF() const
  {
    return *f;
  }

  /** \fn SiconosVector* getFPtr(void) const
   *  \brief get f, the derivative of the state of the DynamicalSystem
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getFPtr() const
  {
    return f;
  }

  /** \fn void setF (const SiconosVector& newValue)
   *  \brief set the value of f to newValue
   *  \param SiconosVector newValue
   */
  void setF(const SiconosVector&);

  /** \fn void setFPtr(SiconosVector* newPtr)
   *  \brief set f to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setFPtr(SiconosVector *);

  // --- JacobianXF ---

  /** \fn  const SimpleMatrix getJacobianXF(void) const
   *  \brief get the value of JacobianXF
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getJacobianXF() const
  {
    return *jacobianXF;
  }

  /** \fn SiconosMatrix* getJacobianXFPtr(void) const
   *  \brief get JacobianXF
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianXFPtr() const
  {
    return jacobianXF;
  }

  /** \fn void setJacobianXF (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianXF to newValue
   *  \param SiconosMatrix newValue
   */
  virtual void setJacobianXF(const SiconosMatrix&);

  /** \fn void setJacobianXFPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianXF to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  virtual void setJacobianXFPtr(SiconosMatrix *newPtr);

  // uSize

  /** \fn const int getUSize(void) const;
   *  \brief to get uSize, size of u
   *  \return the value of uSize
   */
  inline const unsigned int getUSize(void) const
  {
    return uSize;
  }

  /** \fn void setUSize(const unsigned int&)
   *  \brief to set the value of uSize
   *  \param an integer to set the value of uSize
   */
  void setUSize(const unsigned int&);

  // ---  U ---

  /** \fn  const SimpleVector getU(void) const
   *  \brief get the value of u, control term
   *  \return SimpleVector
   * \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getU() const
  {
    return *u;
  }

  /** \fn SiconosVector* getUPtr(void) const
   *  \brief get u, the "control" term
   *  \return pointer on a SiconosVector
   */
  inline SiconosVector* getUPtr() const
  {
    return u;
  }

  /** \fn void setU (const SiconosVector& newValue)
   *  \brief set the value of u to newValue
   *  \param SiconosVector newValue
   */
  void setU(const SiconosVector&);

  /** \fn void setUPtr(SiconosVector* newPtr)
   *  \brief set u to pointer newPtr
   *  \param SiconosVector * newPtr
   */
  void setUPtr(SiconosVector *);

  // --- T ---

  /** \fn  const SimpleMatrix getT(void) const
   *  \brief get the value of T
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getT() const
  {
    return *T;
  }

  /** \fn SiconosMatrix* getTPtr(void) const
   *  \brief get T
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getTPtr() const
  {
    return T;
  }

  /** \fn void setT (const SiconosMatrix& newValue)
   *  \brief set the value of T to newValue
   *  \param SiconosMatrix newValue
   */
  void setT(const SiconosMatrix&);

  /** \fn void setTPtr(SiconosMatrix* newPtr)
   *  \brief set T to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setTPtr(SiconosMatrix *newPtr);

  // --- Steps in memory ---

  /** \fn const int getStepsInMemory(void) const
   *  \brief get the value of stepsInMemory
   *  \return the value of stepsInMemory
   */
  inline const int getStepsInMemory() const
  {
    return stepsInMemory;
  }

  /** \fn void setStepsInMemory(const int&)
   *  \brief set the value of stepsInMemory
   *  \param int steps : the value to set stepsInMemory
   */
  inline void setStepsInMemory(const int& steps)
  {
    stepsInMemory = steps;
  }

  // --- Boundary Conditions ---

  /*\todo: to be finished when BC class will be allright */
  /** \fn  const BoundaryCondition getBoundaryCondition(void) const
   *  \brief get the value of BoundaryCondition
   *  \return an object BoundaryCondition
   */
  //inline BoundaryCondition getBoundaryCondition() const { return *BC; }

  /** \fn BoundaryCondition getBoundaryConditionPtr(void) const
   *  \brief get the BoundaryCondition
   *  \return a pointer on the BoundaryCondition object
   */
  inline BoundaryCondition* getBoundaryConditionPtr() const
  {
    return BC;
  }

  /** \fn void setBoundaryCondition(const BoundaryCondition&)
   *  \brief set the Boundary Conditions
   *  \param ref on an object BoundaryCondition
   */
  //inline void setBoundaryCondition(const BoundaryCondition& newBC) {*BC = newBC; }

  /** \fn void setBoundaryConditionPtr(BoundaryCondition*)
   *  \brief set the BoundaryCondition pointer
   *  \param BoundaryCondition *bc : the BoundaryCondition to set BC
   */
  void setBoundaryConditionPtr(BoundaryCondition *newBC);

  // --- dsxml ---

  /** \fn inline const DynamicalSystemXML* getDynamicalSystemXMLPtr() const
   *  \brief get the object DynamicalSystemXML of the DynamicalSystem
   *  \return a pointer on the DynamicalSystemXML of the DynamicalSystem
   */
  inline const DynamicalSystemXML* getDynamicalSystemXMLPtr() const
  {
    return dsxml;
  }

  /** \fn inline void setDynamicalSystemXMLPtr(DynamicalSystemXML *dsxml)
   *  \brief set the DynamicalSystemXML of the DynamicalSystem
   *  \param DynamicalSystemXML* dsxml : the address of theDynamicalSystemXML to set
   */
  inline void setDynamicalSystemXMLPtr(DynamicalSystemXML *newDsxml)
  {
    dsxml = newDsxml;
  }

  // --- DS input-output ---

  /** \fn vector<DSInputOutput*> getDSInputOutputs(void)
   *  \brief allows to get all the DSInputOutput of the DynamicalSystem
   *  \return the vector of DSInputOutput
   */
  inline std::vector<DSInputOutput*> getDSInputOutputs(void)
  {
    return dsioVector;
  }

  /** \fn DSInputOutput* getDSInputOutput(int)
   *  \brief allows to get one specific DSInputOutput, with its place in the vector of DSInputOutput
   *  \param int : the place of the DSInputOutput in the vector of DSInputOutput of the DynamicalSystem
   *  \return DSInputOutput* : dsioVector[ i ] DSInputOutput
   */
  DSInputOutput* getDSInputOutput(const unsigned int&);

  /** \fn void setDSInputOutputs(vector<DSInputOutput*>)
   *  \brief allows to set all the DSInputOutputs of the DynamicalSystem
   *  \param vector<DSInputOutput*> : the vector to set
   */
  inline void setDSInputOutputs(std::vector<DSInputOutput*> newDsioVect)
  {
    dsioVector = newDsioVect;
  }

  /** \fn void addDSInputOutput(DSInputOutput*)
   *  \brief allows to add the DSInputOutput to the DynamicalSystem
   *  \param DSInputOutput* : the DSInputOutput to add
   */
  inline void addDSInputOutput(DSInputOutput* dsio)
  {
    dsioVector.push_back(dsio);
  }

  // ===== TMP WORK VECTOR =====

  /** \fn  std::map<std::string , SimpleVector*> getTmpWorkVector()
   *  \brief get the vector of temporary saved vector
   *  \return a std vector
   */
  inline std::map<const std::string , SimpleVector*> getTmpWorkVector()
  {
    return tmpWorkVector;
  }

  /** \fn  SimpleVector getTmpWorkVector(const std::string& id)
   *  \brief get a temporary saved vector, ref by id
   *  \return a std vector
   */
  inline SimpleVector* getTmpWorkVector(const std::string & id)
  {
    return tmpWorkVector[id];
  }

  /** \fn void set(map<std::string , SimpleVector*>)
   *  \brief set TmpWorkVector
   *  \param a map<std::string , SimpleVector*>
   */
  inline void setTmpWorkVector(std::map<const std::string , SimpleVector*> newVect)
  {
    tmpWorkVector = newVect;
  }

  /** \fn void addTmpWorkVector(SimpleVector*, const string&)
  *  \brief to add a temporary vector
  *  \param a SimpleVector*
  *  \param a string id
  */
  inline void addTmpWorkVector(SimpleVector* newVal, const std::string& id)
  {
    *tmpWorkVector[id] = *newVal;
  }

  /** \fn void allocateTmpWorkVector(const std::string&, const int&)
   *  \brief to allocate memory for a new vector in tmp map
   *  \param the id of the SimpleVector
   *  \param an int to set the size
   */
  inline void allocateTmpWorkVector(const std::string& id, const int& size)
  {
    tmpWorkVector[id] = new SimpleVector(size);
  }

  /** \fn freeTmpWorkVector(const std::string& )
   *  \brief to free memory in the map
   *  \param the id of the SimpleVector to free
   */
  inline void freeTmpWorkVector(const std::string& id)
  {
    delete tmpWorkVector[id];
  }

  /** \fn void initialize(const double& = 0, const unsigned int& = 1) ;
   *  \brief dynamical system initialization function: mainly set memory and compute value for initial state values.
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  virtual void initialize(const double& = 0, const unsigned int& = 1) ;

  // ===== MEMORY MANAGEMENT FUNCTIONS =====

  /** \fn void initMemory(const unsigned int& steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory
   */
  virtual void initMemory(const unsigned int&) ;

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x and r in the stored previous values
   *  xMemory and rMemory,
   */
  virtual void swapInMemory();

  // ===== COMPUTE PLUGINS FUNCTIONS =====

  // --- getters for plugin functions names ---

  /** \fn  std::string getComputeFFunctionName() const
   *  \brief get name of function that computes f (if f from plugin)
   *  \return a string
   */
  inline const std::string getComputeFFunctionName() const
  {
    return computeFFunctionName;
  }

  /** \fn  std::string getComputeJacobianXFFunctionName() const
   *  \brief get name of function that computes computeJacobianXF (if computeJacobianXF from plugin)
   *  \return a string
   */
  inline const std::string getComputeJacobianXFFunctionName() const
  {
    return computeJacobianXFFunctionName;
  }

  /** \fn  std::string getComputeUFunctionName() const
   *  \brief get name of function that computes u (if u from plugin)
   *  \return a string
   */
  inline const std::string getComputeUFunctionName() const
  {
    return computeUFunctionName;
  }

  /** \fn  std::string getComputeTFunctionName() const
   *  \brief get name of function that computes T (if T from plugin)
   *  \return a string
   */
  inline const std::string getComputeTFunctionName() const
  {
    return computeTFunctionName;
  }

  // --- setters for functions to compute plugins ---

  /** \fn void setComputeFFunction(const string&, const string&)
   *  \brief to set a specified function to compute f(x,t)
   *  \param string pluginPath : the complete path to the plugin
   *  \param string functionName : the function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeFFunction(const std::string & pluginPath, const std::string& functionName);

  /** \fn void setComputeJacobianXFFunction(const string&, const string&)
   *  \brief to set a specified function to compute jacobianXF
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  virtual void setComputeJacobianXFFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeUFunction(const string&, const string&)
   *  \brief to set a specified function to compute u
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeUFunction(const std::string & pluginPath, const std::string & functionName);

  /** \fn void setComputeTFunction(const string&, const string&)
   *  \brief to set a specified function to compute T
   *  \param string pluginPath : the complete path to the plugin
   *  \param the string functionName : function name to use in this library
   *  \exception SiconosSharedLibraryException
   */
  void setComputeTFunction(const std::string & pluginPath, const std::string & functionName);

  // -- parametersList --

  /** \fn map<string, SimpleVector*> getParameters() const
   *  \brief get the full map of parameters
   *  \return a map<string,SimpleVector*>
   */
  inline std::map<std::string, SimpleVector*> getParameters() const
  {
    return parametersList;
  };

  /** \fn  const SimpleVector getParameter(const string & id) const
   *  \brief get the vector of parameters corresponding to plug-in function named id
   *  \return a SimpleVector
   */
  inline const SimpleVector getParameter(const std::string& id)
  {
    return *(parametersList[id]);
  };

  /** \fn SimpleVector* getParameterPtr(const string& id) const
   *  \brief get the pointer to the vector of parameters corresponding to plug-in function named id
   *  \return a pointer on a SimpleVector
   */
  inline SimpleVector* getParameterPtr(const std::string& id)
  {
    return parametersList[id];
  };

  /** \fn void setParameters(const std::map<string, SimpleVector*>& newMap)
   *  \brief set the map for parameters
   *  \param a map<string, SimpleVector*>
   */
  void setParameters(const std::map<std::string, SimpleVector*>&);

  /** \fn void setParameter(const SimpleVector& newValue, const string& id)
   *  \brief set vector corresponding to plug-in function named id to newValue
   *  \param a SimpleVector
   *  \param a string
   */
  void setParameter(const SimpleVector&, const std::string&);

  /** \fn void setParameterPtr(SimpleVector* newPtr, const string& id)
   *  \brief set vector corresponding to plug-in function named id to newPtr (!! pointer link !!)
   *  \param a pointer to SimpleVector
   *  \param a string
   */
  void setParameterPtr(SimpleVector *, const std::string&);

  // --- compute plugin functions ---

  /** \fn void computeF(const double& time)
   * \brief Default function to compute \f$ f: (x,t)\f$
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeF(const double&);

  /** \fn static void computeJacobianXF (const double& time, const bool & =false)
   *  \brief Default function to compute \f$ \nabla_x f: (x,t) \in R^{n} \times R  \mapsto  R^{n \times n} \f$
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  virtual void computeJacobianXF(const double&, const bool & = false);

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

  /** \fn static void computeU (const double&)
   *  \brief Default function to compute u
   * \param double time : current time
   *  \exception RuntimeException
   */
  virtual void computeU(const double&);

  /** \fn static void computeU (const double&)
   *  \brief function to compute u when x is not those of the current object.
   *  \param double time : current time
   *  \param SiconosVector* : pointer to a x value
   *  \exception RuntimeException
   */
  virtual void computeU(const double&,  SiconosVector* xx);

  /** \fn static void computeT ()
   *  \brief Default function to compute T
   *  \exception RuntimeException
   */
  virtual void computeT();

  // --- isPlugin ---

  /** \fn const std::map<std::string, bool> getIsPlugin() const
   *  \brief get isPlugin, map of flags to check if operators are plugged or not
   *  \return a map of bool
   */
  inline const std::map<std::string, bool> getIsPlugin() const
  {
    return isPlugin;
  }

  /** \fn const bool isPlugged(const std::string& name) const
   *  \brief return true if "name" is plugged, else false (ie name is constant)
   *  \return a map of bool
   */
  inline const bool isPlugged(const std::string& name)
  {
    return isPlugin[name];
  }

  // ===== XML MANAGEMENT FUNCTIONS =====

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** \fn void saveDSDataToXML()
   *  \brief copy the data common to each system in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSDataToXML();

  /** \fn void saveBCToXML()
   *  \brief copy the Boundary Conditions data in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveBCToXML();

  /** \fn void saveDSIOToXML()
   *  \brief copy the DS Input-Output data in the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSIOToXML();

  // ===== MISCELLANEOUS ====

  /** \fn void display()
   *  \brief print the data of the dynamical system on the standard output
   */
  virtual void display() const;

  /** \fn void dsConvergenceIndicator()
   *  \brief Default function for computing an indicator of convergence
   *  \return a double when DS is a Lagrangian
   */
  virtual double dsConvergenceIndicator();

};

#endif // DYNAMICALSYSTEM_H



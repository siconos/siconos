#ifndef LAGRANGIANNLDS_H
#define LAGRANGIANNLDS_H

#include "DynamicalSystem.h"
#include "LagrangianDSXML.h"

#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include <iostream>
#include <vector>

#include "SiconosSharedLibrary.h"

//using namespace std;


class LagrangianDSXML;

/** \class LagrangianDS
 *  \brief main class of Lagrangian dynamic systems
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianDS  defines  and computes a generic ndof-dimensional
 * Lagrangian Non Linear Dynamical System of the form :
 * \f[
 * M(q) \ddot q + Q(\dot q, q) = F_{Int}(\dot q , q , t)+F_{Ext}(t) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In the particular case of Non Smooth evolution, the variable p contains the impulse and not the force.
 *    -  \f$ M(q)  \in  R^{ndof \times ndof}  \f$ is the inertia term saved in the SiconosMatrix mass.
 *    -  \f$ Q(\dot q, q)  \in R^{ndof}\f$ is the non linear inertia term saved in the SiconosVector QNLInertia.
 *    -  \f$ F_{Int}(\dot q , q , t)  \in R^{ndof} \f$ are the internal forces saved in the SiconosVector fInt.
 *    -  \f$ F_{Ext}(t)  \in R^{ndof}  \f$ are the external forces saved in the SiconosVector fExt.
 *
 *
 * One word on the initial condition.
 *
 * One word on the bilateral constraint
 *
 * The state of the master class DynamicalSystem is defined by \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$ and then \f$ n= 2 ndof \f$ and the VectorField
 * is specified as :
 * \f[
 * f(x,t) = \left[\begin{array}{cc}
 *  0_{ndof \times ndof} & I_{ndof \times ndof} \\
 * M^{-1}(q)\left[   F_{Int}(\dot q , q , t)+F_{Ext}( q , t) -  Q(\dot q, q) \right]\\
 * \end{array}\right]
 * \f]
 *  and the input due to the non smooth law by
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
 *
 *
 * \todo Automatically, specify the function of DynamicalSystem such as
 *          VectorField.
 * \todo The specification of the master class DynamicalSystem must be made :
 *     -  Automatically,  in order to specify the function of DynamicalSystem such as
 *          VectorField.
 *     -  On the control of the OneStep Integrator and/or the user for the variable x, xDot and r and the memories.
 * we should find a right way ti do that without systematically copying the value of the state.
 * \todo Default plugin mechanism (for example for fExt FInt, M, ....) has to be enhanced following these steps :
 *   -# A Plugin is defined (in the XML or via the API) : Plug this Plugin. If a value is also given  (if it is possible),
 * this value is set as a current value in this object.
 *   -# A contant vector is given : Create a function which returns the constant value and plug it.
 *   -# Nothing is given : Plug the Basic Plugin and print a help message.
 *
 * \warning :  A constant Mass Matrix was loaded from the XML. For the moment ,  constant Mass Matrix
 * is dedicated to LagrangianLinearTIDS
 */
class LagrangianDS : public DynamicalSystem
{
public:

  /** \fn LagrangianDS(DSXML * nsdsXML)
   *  \brief constructor from an xml file
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \exception RuntimeException
   */
  LagrangianDS(DSXML * dsXML);

  /** \fn LagrangianDS(int number, int ndof,
      SiconosVector* q0, SiconosVector* velocity0,
      string fInt, string fExt,
      string jacobianQFInt, string jacobianVelocityFInt,
      string jacobianQQNLInertia, string jacobianVelocityQNLInertia,
      NSDS * nsds)
      *  \brief constructor from a minimum set of data
      *  \param int : the number for this DynamicalSystem
      *  \param int : the dimension of this DynamicalSystem
      *  \param SiconosVector* : initial coordinates of this DynamicalSystem
      *  \param SiconosVector* : initial velocity of this DynamicalSystem
      *  \param NSDS * : The NSDS which contains this DynamicalSystem
      *  \param string : fInt plugin name and location
      *  \param string : fExt plugin name and location
      *  \param string : jacobianQFInt plugin name and location
      *  \param string : jacobianVelocityFInt plugin name and location
      *  \param string : jacobianQQNLInertia plugin name and location
      *  \param string : jacobianVelocityQNLInertia plugin name and location
      *  \param NSDS * : The NSDS which contains this DynamicalSystem
      *  \exception RuntimeException
      */
  LagrangianDS(int number, int ndof,
               SiconosVector* q0, SiconosVector* velocity0,
               string mass = "BasicPlugin:computeMass",
               string fInt = "BasicPlugin:computeFInt", string fExt = "BasicPlugin:computeFExt",
               string jacobianQFInt = "BasicPlugin:computeJacobianQFInt",
               string jacobianVelocityFInt = "BasicPlugin:computeJacobianVelocityFInt",
               string jacobianQQNLInertia = "BasicPlugin:computeJacobianQQNLInertia",
               string jacobianVelocityQNLInertia = "BasicPlugin:computeJacobianVelocityQNLInertia",
               string QNLlInertia = "BasicPlugin:computeQNLInertia");

  virtual ~LagrangianDS();

  // --- GETTERS AND SETTERS ---

  /** \fn const int getNdof() const
   *  \brief allows to get the value of ndof
   *  \return the value of ndof
   */
  inline const int getNdof() const
  {
    return ndof;
  };

  /** \fn void setNdof(const int&)
   *  \brief allows to set ndof
   *  \param int ndof : the value to set ndof
   */
  inline void setNdof(const int& newNdof)
  {
    ndof = newNdof;
  };

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
  inline void setQ(const SimpleVector& newValue)
  {
    *q = newValue;
  }

  /** \fn void setQPtr(SimpleVector* newPtr)
   *  \brief set Q to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setQPtr(SimpleVector *newPtr)
  {
    delete q;
    q = newPtr;
  }

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
  inline void setQ0(const SimpleVector& newValue)
  {
    *q0 = newValue;
  }

  /** \fn void setQ0Ptr(SimpleVector* newPtr)
   *  \brief set Q0 to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setQ0Ptr(SimpleVector *newPtr)
  {
    delete q0;
    q0 = newPtr;
  }

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
  inline void setQFree(const SimpleVector& newValue)
  {
    *qFree = newValue;
  }

  /** \fn void setQFreePtr(SimpleVector* newPtr)
   *  \brief set QFree to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setQFreePtr(SimpleVector *newPtr)
  {
    delete qFree;
    qFree = newPtr;
  }

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
  inline void setVelocity(const SimpleVector& newValue)
  {
    *velocity = newValue;
  }

  /** \fn void setVelocityPtr(SimpleVector* newPtr)
   *  \brief set Velocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setVelocityPtr(SimpleVector *newPtr)
  {
    delete velocity;
    velocity = newPtr;
  }

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
  inline void setVelocity0(const SimpleVector& newValue)
  {
    *velocity0 = newValue;
  }

  /** \fn void setVelocity0Ptr(SimpleVector* newPtr)
   *  \brief set Velocity0 to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setVelocity0Ptr(SimpleVector *newPtr)
  {
    delete velocity0;
    velocity0 = newPtr;
  }

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
  inline void setVelocityFree(const SimpleVector& newValue)
  {
    *velocityFree = newValue;
  }

  /** \fn void setVelocityFreePtr(SimpleVector* newPtr)
   *  \brief set VelocityFree to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setVelocityFreePtr(SimpleVector *newPtr)
  {
    delete velocityFree;
    velocityFree = newPtr;
  }

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
  inline void setP(const SimpleVector& newValue)
  {
    *p = newValue;
  }

  /** \fn void setPPtr(SimpleVector* newPtr)
   *  \brief set P to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setPPtr(SimpleVector *newPtr)
  {
    delete p;
    p = newPtr;
  }

  // --- Memory ---

  // -- q memory --

  /** \fn  const SiconosMemory getQMemory() const
   *  \brief get the value of qMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getQMemory() const
  {
    return qMemory;
  }

  /** \fn SiconosMemory getQMemoryPtr() const
   *  \brief get all the values of the state vector q stored in memory
   *  \return the memory object which stores previous values of q
   */
  inline SiconosMemory* getQMemoryPtr()
  {
    return &qMemory;
  }

  /** \fn void setQMemory(const SiconosMemory &)
   *  \brief set the value of qMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setQMemory(const SiconosMemory& qMem)
  {
    qMemory = qMem;
  }

  // -- velocity memory --

  /** \fn  const SiconosMemory getVelocityMemory() const
   *  \brief get the value of velocityMemory
   *  \return a SiconosMemory
   */
  inline const SiconosMemory getVelocityMemory() const
  {
    return velocityMemory;
  }

  /** \fn SiconosMemory getVelocityMemoryPtr() const
   *  \brief get all the values of the state vector velocity stored in memory
   *  \return the memory object which stores previous values of velocity
   */
  inline SiconosMemory* getVelocityMemoryPtr()
  {
    return &velocityMemory;
  }

  /** \fn void setVelocityMemory(const SiconosMemory &)
   *  \brief set the value of velocityMemory
   *  \param a ref on a SiconosMemory
   */
  inline void setVelocityMemory(const SiconosMemory& velocityMem)
  {
    velocityMemory = velocityMem;
  }

  // -- Mass --

  /** \fn  const SiconosMatrix getMass() const
   *  \brief get the value of Mass
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getMass() const
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
  inline void setMass(const SiconosMatrix& newValue)
  {
    *mass = newValue;
  }

  /** \fn void setMassPtr(SiconosMatrix* newPtr)
   *  \brief set Mass to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setMassPtr(SiconosMatrix *newPtr)
  {
    delete mass;
    mass = newPtr;
  }

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
  inline void setFInt(const SimpleVector& newValue)
  {
    *fInt = newValue;
  }

  /** \fn void setFIntPtr(SimpleVector* newPtr)
   *  \brief set FInt to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setFIntPtr(SimpleVector *newPtr)
  {
    delete fInt;
    fInt = newPtr;
  }

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
  inline void setFExt(const SimpleVector& newValue)
  {
    *fExt = newValue;
  }

  /** \fn void setFExtPtr(SimpleVector* newPtr)
   *  \brief set FExt to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setFExtPtr(SimpleVector *newPtr)
  {
    delete fExt;
    fExt = newPtr;
  }

  // -- QNLInertia --

  /** \fn  const SimpleVector getQNLInertia() const
   *  \brief get the value of QNLInertia
   *  \return SimpleVector
   */
  inline const SimpleVector getQNLInertia() const
  {
    return *QNLInertia;
  }

  /** \fn SimpleVector* getQNLInertiaPtr() const
   *  \brief get QNLInertia
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQNLInertiaPtr() const
  {
    return QNLInertia;
  }

  /** \fn void setQNLInertia (const SimpleVector& newValue)
   *  \brief set the value of QNLInertia to newValue
   *  \param SimpleVector newValue
   */
  inline void setQNLInertia(const SimpleVector& newValue)
  {
    *QNLInertia = newValue;
  }

  /** \fn void setQNLInertiaPtr(SimpleVector* newPtr)
   *  \brief set QNLInertia to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  inline void setQNLInertiaPtr(SimpleVector *newPtr)
  {
    delete QNLInertia;
    QNLInertia = newPtr;
  }

  // -- Jacobian Q Fint --

  /** \fn  const SiconosMatrix getJacobianQFInt() const
   *  \brief get the value of JacobianQFInt
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianQFInt() const
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
  inline void setJacobianQFInt(const SiconosMatrix& newValue)
  {
    *jacobianQFInt = newValue;
  }

  /** \fn void setJacobianQFIntPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianQFInt to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setJacobianQFIntPtr(SiconosMatrix *newPtr)
  {
    delete jacobianQFInt;
    jacobianQFInt = newPtr;
  }

  // -- Jacobian velocity Fint --

  /** \fn  const SiconosMatrix getJacobianVelocityFInt() const
   *  \brief get the value of JacobianVelocityFInt
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianVelocityFInt() const
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
  inline void setJacobianVelocityFInt(const SiconosMatrix& newValue)
  {
    *jacobianVelocityFInt = newValue;
  }

  /** \fn void setJacobianVelocityFIntPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianVelocityFInt to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setJacobianVelocityFIntPtr(SiconosMatrix *newPtr)
  {
    delete jacobianVelocityFInt;
    jacobianVelocityFInt = newPtr;
  }

  // -- Jacobian Q QNLInertia --

  /** \fn  const SiconosMatrix getJacobianQQNLInertia() const
   *  \brief get the value of JacobianQQNLInertia
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianQQNLInertia() const
  {
    return *jacobianQQNLInertia;
  }

  /** \fn SiconosMatrix* getJacobianQQNLInertiaPtr() const
   *  \brief get JacobianQQNLInertia
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianQQNLInertiaPtr() const
  {
    return jacobianQQNLInertia;
  }

  /** \fn void setJacobianQQNLInertia (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianQQNLInertia to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setJacobianQQNLInertia(const SiconosMatrix& newValue)
  {
    *jacobianQQNLInertia = newValue;
  }

  /** \fn void setJacobianQQNLInertiaPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianQQNLInertia to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setJacobianQQNLInertiaPtr(SiconosMatrix *newPtr)
  {
    delete jacobianQQNLInertia;
    jacobianQQNLInertia = 0;
    jacobianQQNLInertia = newPtr;
  }

  // -- Jacobian velocity QNLInertia --

  /** \fn  const SiconosMatrix getJacobianVelocityQNLInertia() const
   *  \brief get the value of JacobianVelocityQNLInertia
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getJacobianVelocityQNLInertia() const
  {
    return *jacobianVelocityQNLInertia;
  }

  /** \fn SiconosMatrix* getJacobianVelocityQNLInertiaPtr() const
   *  \brief get JacobianVelocityQNLInertia
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getJacobianVelocityQNLInertiaPtr() const
  {
    return jacobianVelocityQNLInertia;
  }

  /** \fn void setJacobianVelocityQNLInertia (const SiconosMatrix& newValue)
   *  \brief set the value of JacobianVelocityQNLInertia to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setJacobianVelocityQNLInertia(const SiconosMatrix& newValue)
  {
    *jacobianVelocityQNLInertia = newValue;
  }

  /** \fn void setJacobianVelocityQNLInertiaPtr(SiconosMatrix* newPtr)
   *  \brief set JacobianVelocityQNLInertia to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  inline void setJacobianVelocityQNLInertiaPtr(SiconosMatrix *newPtr)
  {
    delete jacobianVelocityQNLInertia;
    jacobianVelocityQNLInertia = newPtr;
  }

  // --- PLUGINS RELATED FUNCTIONS ---

  /** \fn void computeMass(double time)
   *  \brief default function to compute the mass
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeMass(double time);

  /** \fn void computeMass(double time, SimpleVector *q)
   *  \brief function to compute the mass
   *  \param double time : the current time, SimpleVector*: pointer on the state vector q
   *  \exception RuntimeException
   */
  virtual void computeMass(double time, SimpleVector *q);

  /** \fn void computeFInt(double time)
   *  \brief default function to compute the internal strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeFInt(double time);

  /** \fn void computeFInt(double time, SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the internal strengths
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  virtual void computeFInt(double time, SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeFExt(double time)
   *  \brief default function to compute the external strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */

  virtual void computeFExt(double time);

  /** \fn void computeQNLInertia();
   *  \brief default function to compute the inertia
   *  \exception RuntimeException
   */
  virtual void computeQNLInertia();

  /** \fn void computeQNLInertia(SimpleVector q, SimpleVector velocity);
   *  \brief function to compute the inertia
   *  \param SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  virtual void computeQNLInertia(SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianQFInt(double time)
   *  \brief default function to compute the gradient of the internal strengths compared to the state
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianQFInt(double time);

  /** \fn void computeJacobianQFInt(double time,SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the gradient of the internal strengths compared to state q
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  virtual void computeJacobianQFInt(double time, SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianVelocityFInt(double time)
   *  \brief function to compute the gradient of the internal strengths compared to velocity
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianVelocityFInt(double time);

  /** \fn void computeJacobianVelocityFInt(double time, SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the gradient of the internal strengths compared to velocity
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  virtual void computeJacobianVelocityFInt(double time, SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianQQNLInertia(double time)
   *  \brief function to compute the gradient of the inertia strengths compared to the state q
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianQQNLInertia(double time);

  /** \fn void computeJacobianQQNLInertia(double time,SimpleVector q, SimpleVector velocity)
   *  \brief function to compute the gradient of the inertia strengths compared to the state q
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  virtual void computeJacobianQQNLInertia(double time, SimpleVector *q, SimpleVector *velocity);

  /** \fn void computeJacobianVelocityQNLInertia(double time )
   *  \brief function to compute the gradient of the inertia strengths compared to velocity
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianVelocityQNLInertia(double time);
  /** \fn void computeJacobianVelocityQNLInertia(double time, SimpleVector q, SimpleVector velocity )
   *  \brief function to compute the gradient of the inertia strengths compared to velocity
   *  \param double time : the current time, SimpleVector*: pointers on the state vectors q and velocity (\dot q)
   *  \exception RuntimeException
   */
  virtual void computeJacobianVelocityQNLInertia(double time, SimpleVector *q, SimpleVector *velocity);

  /** \fn void setComputeMassFunction(const string pluginPath, const string functionName&)
   *  \brief allow to set a specified function to compute the mass
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeMassFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeFIntFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFIntFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeFExtFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception ICDLL_CSharedLibraryException
   */
  void setComputeFExtFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeQNLInertiaFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the inertia
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosCSharedLibraryException
   */
  void setComputeQNLInertiaFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeJacobianQFIntFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the gradient of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQFIntFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeJacobianVelocityFIntFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the internal strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityFIntFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeJacobianQQNLInertiaFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the gradient of the the external strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQQNLInertiaFunction(const string& pluginPath, const string& functionName);

  /** \fn void setComputeJacobianVelocityQNLInertiaFunction(const string& pluginPath, const string& functionName)
   *  \brief allow to set a specified function to compute the external strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityQNLInertiaFunction(const string& pluginPath, const string& functionName);

  // --- miscellaneous ---

  /** \fn void SiconosVectorSizeInit()
   *  \brief Initialisation of all the vector of the dynamical system with the right size
   *  \exception RuntimeException
   */
  //virtual void SiconosVectorSizeInit();

  /** \fn void CompositeVectorInit()
   *  \brief Initialisation of all the composite vector of the dynamical system with different SiconosVector composing each composite
   *  \exception RuntimeException
   */
  //virtual void CompositeVectorInit();

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  virtual void display() const;

  /** \fn void initMemory(const int& steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  void initMemory(const int& steps);

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, xDot and r in the stored previous values
   *  xMemory, xDotMemory, rMemory,
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
  double LagrangianDS::dsConvergenceIndicator() const ;

protected:

  // -- DEFAULT CONSTRUCTOR --
  /** \fn LagrangianDS()
   *  \brief Default constructor
   */
  LagrangianDS();

  // -- MEMBERS --
  /** number of degrees of freedom of the system */
  int ndof;
  /** coordinates of the system */
  SimpleVector *q;
  /** initial coordinates of the system */
  SimpleVector *q0;
  /** free coordinate */
  SimpleVector *qFree;
  /** memory of previous coordinates of the system */
  SiconosMemory qMemory;
  /** velocity of the system */
  SimpleVector *velocity;
  /** initial velocity of the system */
  SimpleVector *velocity0;
  /** free Velocity */
  SimpleVector *velocityFree;
  /** memory of previous velocity of the system */
  SiconosMemory velocityMemory;
  /** Reaction due to the non smooth law */
  SimpleVector *p;

  /** mass of the system */
  SiconosMatrix *mass;
  /** internal strength of the system */
  SimpleVector *fInt;
  /** external strength of the system */
  SimpleVector *fExt;
  /* contains the name of the plugin used to compute fInt */
  string fIntFunctionName;
  /* contains the name of the plugin used to compute fExt */
  string fExtFunctionName;
  /* contains the name of the plugin used to compute the mass */
  string massFunctionName;
  /* contains the name of the plugin used to compute jacobianQFInt */
  string jacobianQFIntFunctionName;
  /* contains the name of the plugin used to compute jacobianQQNLInertia */
  string jacobianQQNLInertiaFunctionName;
  /* contains the name of the plugin used to compute jacobianVelocityFInt */
  string jacobianVelocityFIntFunctionName;
  /* contains the name of the plugin used to compute jacobianVelocityQNLInertia */
  string jacobianVelocityQNLInertiaFunctionName;
  /* contains the name of the plugin used to compute QNLInertia */
  string QNLInertiaFunctionName;
  /** non-linear inertia term of the system */
  SimpleVector *QNLInertia;
  /** jacobian/coordinates of internal strength */
  SiconosMatrix *jacobianQFInt;
  /** jacobian/velocity of internal strength */
  SiconosMatrix *jacobianVelocityFInt;
  /** jacobian/coordinates of inertia */
  SiconosMatrix *jacobianQQNLInertia;
  /** jacobian/velocity of inertie */
  SiconosMatrix *jacobianVelocityQNLInertia;

  /** class for manage plugin (open, close librairy...) */
  SiconosSharedLibrary cShared;

  //////////////////////////////////////

  /** \fn void (*computeMassPtr)(double time, double* qPtr, int sizeOfq, double* massPtr)
   *  \brief compute the mass
   *  \param int* sizeOfq : the size of the vector q
   *  \param double* time : the time for the computation
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* massPtr : the pointer to the first element of the matrix mass (in-out parameter)
   */
  void (*computeMassPtr)(int* sizeOfq, double* time, double* qPtr, double* massPtr);

  /** \fn void (*computeFIntPtr)(int* sizeOfq, double* time, double* qPtr, double* velocityPtr, double* fIntPtr)
   *  \brief computes the internal strengths
   *  \param int* sizeOfq : the size of the vector q
   *  \param int* time : the current time
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* fIntPtr : the pointer to the first element of the vector FInt (in-out parameter)
   */
  void (*computeFIntPtr)(int* sizeOfq, double* time, double* qPtr, double* velocityPtr, double* fIntPtr);

  /** \fn void (*computeFExtPtr)(double* time, double* qPtr, int* sizeOfq, double* fExtPtr)
   *  \brief computes the external strengths
   *  \param int* sizeOfq : the size of the vector q
   *      \param double* time : the current time
   *  \param double* fExtPtr : the pointer to the first element of the vector FInt (in-out parameter)
   */
  //void (*computeFExtPtr)(int* sizeOfq, double* time, double* qPtr, double* fExtPtr);
  void (*computeFExtPtr)(int* sizeOfq, double* time, double* fExtPtr);

  /** \fn void (*computeQNLInertiaPtr)(int* sizeOfq, double* qPtr, double* velocityPtr, double* QNLInertiaPtr)
   *  \brief computes the inertia
   *  \param int* sizeOfq : the size of the vector q
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* QNLInertiaPtr : the pointer to the first element of the vector QNLInertia (in-out parameter)
   */
  void (*computeQNLInertiaPtr)(int* sizeOfq, double* qPtr, double* velocityPtr, double* QNLInertiaPtr);

  /** \fn void (*computeJacobianQFIntPtr)(int* sizeOfq, double* time, double* qPtr, double* velocityPtr, double* jacobPtr)
   *  \brief computes the gradient of the the internal strength compared to the state
   *  \param int* sizeOfq : the size of the vector q
   *  \param double* time : the current time
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* velocityPtr : the pointer to the first element of the vector velocity
   *  \param double* jacobPtr : the pointer to the first element of the matrix JacobianQFInt (in-out parameter)
   */
  void (*computeJacobianQFIntPtr)(int* sizeOfq, double* time, double* qPtr, double* velocityPtr, double* jacobPtr);

  /** \fn void (*computeJacobianVelocityFIntPtr)(double time)
   *  \brief computes the gradient of the the internal strength compared to the velocity
   *  \param to be defined
   */
  void (*computeJacobianVelocityFIntPtr)(int* sizeOfq, double* time, double* qPtr, double* velocityPtr, double* jacobPtr);

  /** \fn void (*computeJacobianQQNLInertia)(double time)
   *  \brief computes the gradient of the the external strength compared to the state
   *  \param to be defined
   */
  void (*computeJacobianQQNLInertiaPtr)(int* sizeOfq, double* qPtr, double* velocityPtr, double* jacobPtr);

  /** \fn void (*computeJacobianVelocityQNLInertiaPtr)(double time)
   *  \brief computes the gradient of the the external strength compared to the velocity
   *  \param to be defined
   */
  void (*computeJacobianVelocityQNLInertiaPtr)(int* sizeOfq, double* qPtr, double* velocityPtr, double* jacobPtr);

};

#endif // LAGRANGIANNLDS_H

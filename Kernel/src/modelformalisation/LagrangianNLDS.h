#ifndef LAGRANGIANNLDS_H
#define LAGRANGIANNLDS_H

#include "DynamicalSystem.h"
#include "LagrangianNLDSXML.h"

#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include <iostream>
#include <vector>

#include "SiconosSharedLibrary.h"

using namespace std;


class LagrangianNLDSXML;

/** \class LagrangianNLDS
 *  \brief main class of Lagrangian dynamic systems
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianNLDS  allows to define  and compute a generic ndof-dimensional
 * Lagrangian Non Linear Dynamical System of the form :
 * \f[
 * M(q) \ddot q + Q(\dot q, q) = F_{Int}(\dot q , q , t)+F_{Ext}( q , t) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,  the variable p stored the impulse and not the force.
 *    -  \f$ M(q)  \in  R^{ndof \times ndof}  \f$ is the inertia term stored in the SiconosMatrix mass.
 *    -  \f$ Q(\dot q, q)  \in R^{ndof}\f$ the non linear inertia term stored in the SiconosVector QNLInertia.
 *    -  \f$ F_{Int}(\dot q , q , t)  \in R^{ndof} \f$ the internal forces stored in the SiconosVector fInt.
 *    -  \f$ F_{Ext}( q , t)  \in R^{ndof}  \f$ the external forces stored in the SiconosVector fExt.
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
 * is dedicated to LagrangianTIDS
 */
class LagrangianNLDS : public DynamicalSystem
{
public:

  /** \fn LagrangianNLDS()
   *  \brief Default constructor
   */
  LagrangianNLDS();

  /** \fn LagrangianNLDS(DSXML*)
   *  \brief constructor with XML object of the LagrangianNLDS
   *  \param DSXML* : the XML object corresponding
   */
  LagrangianNLDS(DSXML*);

  virtual ~LagrangianNLDS();

  /** \fn void initMemory(int steps) ;
   *  \brief initialize the SiconosMemory objects with a positive size.
   *  \param the size of the SiconosMemory. must be >= 0
   */
  void initMemory(int steps);

  /** \fn virtual void swapInMemory(void);
   * \brief push the current values of x, xDot and r in the stored previous values
   *  xMemory, xDotMemory, rMemory,
   * \todo Modify the function swapIn Memory with the new Object Memory
   */
  void swapInMemory(void);

  // getter/setter
  /** \fn int getNdof(void)
   *  \brief allows to get the value of ndof
   *  \return the value of ndof
   */
  inline int getNdof(void) const
  {
    return this->ndof;
  };

  /** \fn SiconosMatrix getMass(void)
   *  \brief allows to get the SiconosMatrix mass
   *  \return the SiconosMatrix mass
   */
  inline SiconosMatrix getMass(void) const
  {
    return this->mass;
  };

  /** \fn SimpleVector getQ(void)
   *  \brief allows to get the SiconosVector q
   *  \return SimpleVector : value of q
   */
  inline SimpleVector const getQ(void)
  {
    return this->q;
  };

  /** \fn SimpleVector getQ0(void)
   *  \brief allows to get the SiconosVector q0
   *  \return SimpleVector : value of q0
   */
  inline SimpleVector const getQ0(void)
  {
    return this->q0;
  };

  /** \fn SiconosMatrix* getMassPtr(void)
   *  \brief allows to get the SiconosMatrix* mass
   *  \return the SiconosMatrix* mass
   */
  SiconosMatrix* getMassPtr(void);

  /** \fn SimpleVector* getQPtr(void)
   *  \brief allows to get the SiconosVector* q
   *  \return SimpleVector* : pointer on q
   */
  SimpleVector* getQPtr(void);

  /** \fn SimpleVector* getQ0Ptr(void)
   *  \brief allows to get the SiconosVector* q0
   *  \return SimpleVector* : pointer on q0
   */
  SimpleVector* getQ0Ptr(void);

  /** \fn SiconosMemory* getQMemories(void)
   *  \brief allows to get all the value of old q
   *  \return the memory object containing previous values of q
   */
  SiconosMemory* getQMemories(void);

  /** \fn SimpleVector getVelocity(void)
   *  \brief allows to get the SiconosVector velocity
   *  \return SimpleVector : value of velocity
   */
  inline SimpleVector const getVelocity(void)
  {
    return this->velocity;
  };

  /** \fn SimpleVector getVelocity0(void)
   *  \brief allows to get the initial velocity
   *  \return SimpleVector : initial value of velocity
   */
  inline SimpleVector const getVelocity0(void)
  {
    return this->velocity0;
  };

  /** \fn SimpleVector* getVelocityPtr(void)
   *  \brief allows to get the SiconosVector* velocity
   *  \return SimpleVector* : pointer on velocity
   */
  SimpleVector* getVelocityPtr(void);

  /** \fn SimpleVector* getVelocity0Ptr(void)
   *  \brief allows to get the SiconosVector velocity0 Pointer
   *  \return SimpleVector* : pointer on initial velocity
   */
  SimpleVector* getVelocity0Ptr(void);

  /** \fn SiconosMemory* getVelocityMemories(void)
   *  \brief allows to get all the value of old velocity Pointer
   *  \return the memory object containing previous values of velocity
   */
  SiconosMemory* getVelocityMemories(void);

  /** \fn SimpleVector* getVelocityFreePtr(void)
   *  \brief get the SimpleVector* velocityFree
   *  \return SimpleVector* : pointer on velocityFree
   */
  inline SimpleVector* getVelocityFreePtr(void)
  {
    return (&(this->velocityFree));
  }

  /** \fn SimpleVector getVelocityFree(void)
   *  \brief get the free velocity
   *  \return SimpleVector : value of velocityFree
   */
  inline SimpleVector getVelocityFree(void) const
  {
    return this->velocityFree;
  }

  /** \fn void setVelocityFree(const SimpleVector& velocityFree)
   *  \brief set the free velocity
   *  \param SimpleVector& : new value of velocityFree
   */
  inline void setVelocityFree(const SimpleVector& velocityFree)
  {
    this->velocityFree = velocityFree;
  }

  /** \fn SimpleVector* getQFreePtr(void)
   *  \brief get the free state
   *  \return SimpleVector* : pointer on qFree
   */
  inline SimpleVector* getQFreePtr(void)
  {
    return (&(this->qFree));
  }

  /** \fn SimpleVector getQFree(void)
   *  \brief get the free state qFree
   *  \return SimpleVector : value of qFree
   */
  inline SimpleVector getQFree(void)
  {
    return this->qFree;
  }

  /** \fn void setQFree(const SimpleVector& qFree)
   *  \brief set the free state qFree
   *  \param SimpleVector& : new value of qFree
   */
  inline void setQFree(const SimpleVector &qFree)
  {
    this->qFree = qFree;
  }

  /** \fn SimpleVector* getPPtr(void)
   *  \brief get the vector p
   *  \return SimpleVector* : pointer on p
   */
  inline SimpleVector* getPPtr(void)
  {
    return (&(this->p));
  }

  /** \fn SimpleVector getP(void)
   *  \brief get the vector p
   *  \return SimpleVector : value of p
   */
  inline SimpleVector getP(void)
  {
    return this->p;
  }

  /** \fn void setP(const SimpleVector& p)
   *  \brief set the vector p
   *  \param SimpleVector&  : new value of p
   */
  inline void setP(const SimpleVector& p)
  {
    this->p = p;
  }

  /** \fn SimpleVector getFInt(void)
   *  \brief get vector of internal forces
   *  \return SimpleVector : value of fInt
   */
  inline SimpleVector getFInt(void)
  {
    return this->fInt;
  };

  /** \fn SimpleVector* getFIntPtr(void)
   *  \brief get vector of internal forces
   *  \return SimpleVector* : pointer on fInt
   */
  SimpleVector* getFIntPtr(void);

  /** \fn SimpleVector getFExt(void)
   *  \brief get vector of external forces
   *  \return SimpleVector : value of fExt
   */
  inline SimpleVector getFExt(void) const
  {
    return this->fExt;
  };

  /** \fn SimpleVector* getFExtPtr(void)
   *  \brief get vector of external forces
   *  \return SimpleVector* : pointer on fExt
   */
  SimpleVector* getFExtPtr(void);

  /** \fn SimpleVector getQNLInertia(void)
   *  \brief get the inertia
   *  \return SimpleVector : value of QNLInertia
   */
  inline SimpleVector getQNLInertia(void) const
  {
    return this->QNLInertia;
  };

  /** \fn SimpleVector* getQNLInertiaPtr(void)
   *  \brief get the inertia
   *  \return SimpleVector* : pointer on QNLInertia
   */
  SimpleVector* getQNLInertiaPtr(void);

  /** \fn SiconosMatrix getJacobianQFInt(void)
   *  \brief allows to get the SiconosMatrix jacobianQFIntMat
   *  \return the SiconosMatrix jacobianQFIntMat
   */
  inline SiconosMatrix getJacobianQFInt(void) const
  {
    return this->jacobianQFInt;
  };

  /** \fn SiconosMatrix getJacobianVelocityFInt(void)
   *  \brief allows to get  the SiconosMatrix jacobianVelocityFIntMat
   *  \return the SiconosMatrix jacobianVelocityFIntMat
   */
  SiconosMatrix getJacobianVelocityFInt(void) const
  {
    return this->jacobianQFInt;
  };

  /** \fn SiconosMatrix getJacobianQQNLInertia(void)
   *  \brief allows to get the SiconosMatrix JacobianQQNLInertiaMat
   *  \return the SiconosMatrix JacobianQQNLInertiaMat
   */
  SiconosMatrix getJacobianQQNLInertia(void) const
  {
    return this->jacobianQQNLInertia;
  };

  /** \fn SiconosMatrix getJacobianVelocityQNLInertia(void)
   *  \brief allows to get the SiconosMatrix jacobianVelocityQNLInertiaMat
   *  \return the SiconosMatrix jacobianVelocityQNLInertiaMat
   */
  SiconosMatrix getJacobianVelocityQNLInertia(void) const
  {
    return this->jacobianVelocityQNLInertia;
  };

  /** \fn SiconosMatrix* getJacobianQFIntPtr(void)
   *  \brief allows to get the SiconosMatrix* jacobianQFIntMat
   *  \return the SiconosMatrix* jacobianQFIntMat
   */
  SiconosMatrix* getJacobianQFIntPtr(void);

  /** \fn SiconosMatrix* getJacobianVelocityFIntPtr(void)
   *  \brief allows to get  the SiconosMatrix* jacobianVelocityFIntMat
   *  \return the SiconosMatrix* jacobianVelocityFIntMat
   */
  SiconosMatrix* getJacobianVelocityFIntPtr(void);

  /** \fn SiconosMatrix* getJacobianQQNLInertiaPtr(void)
   *  \brief allows to get the SiconosMatrix* JacobianQQNLInertiaMat
   *  \return the SiconosMatrix* JacobianQQNLInertiaMat
   */
  SiconosMatrix* getJacobianQQNLInertiaPtr(void);

  /** \fn SiconosMatrix* getJacobianVelocityQNLInertiaPtr(void)
   *  \brief allows to get the SiconosMatrix* jacobianVelocityQNLInertiaMat
   *  \return the SiconosMatrix* jacobianVelocityQNLInertiaMat
   */
  SiconosMatrix* getJacobianVelocityQNLInertiaPtr(void);

  /** \fn void setNdof(int)
   *  \brief allows to set ndof
   *  \param int ndof : the value to set ndof
   */
  inline void setNdof(const int ndof)
  {
    this->ndof = ndof;
  };

  /** \fn void setMass(SiconosMatrix)
   *  \brief allows to set Mass
   *  \param SiconosMatrix mass : the SiconosMatrix to set Mass
   */
  inline void setMass(const SiconosMatrix &mass)
  {
    this->mass = mass;
  };

  /** \fn void setQ(SimpleVector&)
   *  \brief set the state q
   *  \param SimpleVector& : new value of q
   */
  inline void setQ(const SimpleVector& q)
  {
    this->q = q;
  };

  /** \fn void setQ0(SimpleVector&)
   *  \brief set the initial state q0
   *  \param SimpleVector& : value of initial state q0
   */
  inline void setQ0(const SimpleVector& q0)
  {
    this->q0 = q0;
  };

  /** \fn void setQMemories(SiconosMemory&)
   *  \brief set the memory object of previous values of state q
   *  \param SiconosMemory& : memory object
   */
  inline void setQMemories(const SiconosMemory &mem)
  {
    this->qMemory = qMemory;
  };

  /** \fn void setVelocity(SimpleVector&)
   *  \brief set velocity
   *  \param SimpleVector& : new value of velocity
   */
  inline void setVelocity(const SimpleVector& velocity)
  {
    this->velocity = velocity;
  };

  /** \fn void setVelocity0(SimpleVector&)
   *  \brief set initial velocity
   *  \param SimpleVector& : value of velocity0
   */
  inline void setVelocity0(const SimpleVector& velocity0)
  {
    this->velocity0 = velocity0;
  };

  /** \fn void setVelocityMemories(SiconosMemory&)
   *  \brief set the memory object of previous values of velocity
   *  \param SiconosMemory& : memory object
   */
  inline void setVelocityMemories(const SiconosMemory& mem)
  {
    this->velocityMemory = mem;
  };


  /** \fn void setFInt(SimpleVector&)
   *  \brief set internal forces fInt
   *  \param SimpleVector& : new value of fint
   */
  inline void setFInt(const SimpleVector& fint)
  {
    this->fInt = fint;
  };

  /** \fn void setFExt(SimpleVector&)
   *  \brief set external forces fExt
   *  \param SimpleVector& : new value of fext
   */
  inline void setFExt(const SimpleVector& fext)
  {
    this->fExt = fext;
  };

  /** \fn void setQNLInertia(SimpleVector&)
   *  \brief set inertia
   *  \param SimpleVector& : new value of QNLInertia
   */
  inline void setQNLInertia(const SimpleVector& QNLInertia)
  {
    this->QNLInertia = QNLInertia;
  };


  /** \fn void setJacobianQFInt(SiconosMatrix)
   *  \brief allows to set jacobianQFIntMat
   *  \param SiconosMatrix jacob : the SiconosMatrix to set jacobianQFIntMat
   */
  inline void setJacobianQFInt(const SiconosMatrix &jacob)
  {
    this->jacobianQFInt = jacob;
  };

  /** \fn void setJacobianVelocityFInt(SiconosMatrix)
   *  \brief allows to set jacobianVelocityFIntMat
   *  \param SiconosMatrix jacob : the SiconosMatrix to set jacobianVelocityFIntMat
   */
  inline void setJacobianVelocityFInt(const SiconosMatrix &jacob)
  {
    this->jacobianVelocityFInt = jacob;
  };

  /** \fn void setJacobianVelocityFInt(SiconosMatrix)
   *  \brief allows to set JacobianQQNLInertiaMat
   *  \param SiconosMatrix jacob : the SiconosMatrix to set JacobianQQNLInertiaMat
   */
  inline void setJacobianQQNLInertia(const SiconosMatrix &jacob)
  {
    this->jacobianQQNLInertia = jacob;
  };

  /** \fn void setJacobianVelocityQNLInertia(SiconosMatrix)
   *  \brief allows to set jacobianVelocityQNLInertiaMat
   *  \param SiconosMatrix jacob : the SiconosMatrix to set jacobianVelocityQNLInertiaMat
   */
  inline void setJacobianVelocityQNLInertia(const SiconosMatrix &jacob)
  {
    this->jacobianVelocityQNLInertia = jacob;
  };


  /** \fn void fillDSWithDSXML()
   *  \brief overload of the function for a LagrangianNLDS
   *  \exception RuntimeException
   */
  virtual void fillDSWithDSXML();

  /** \fn void SiconosVectorSizeInit()
   *  \brief Initialisation of all the vector of the dynamical system with the right size
   *  \exception RuntimeException
   */
  virtual void SiconosVectorSizeInit();

  /** \fn void CompositeVectorInit()
   *  \brief Initialisation of all the composite vector of the dynamical system with different SiconosVector composing each composite
   *  \exception RuntimeException
   */
  virtual void CompositeVectorInit();


  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  virtual void saveDSToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  virtual void display() const;

  //////////////////////////////////////


  /** \fn void computeMass(double time)
   *  \brief default function to compute the mass
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeMass(double time);

  /** \fn void computeFInt(double time)
   *  \brief default function to compute the internal strengths
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeFInt(double time);

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

  /** \fn void computeJacobianQFInt(double time)
   *  \brief default function to compute the gradient of the internal strength compared to the state
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianQFInt(double time);

  /** \fn void computeJacobianVelocityFInt(double time)
   *  \brief default function to compute the gradient of the internal strength compared to the velocity
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianVelocityFInt(double time);

  /** \fn void computeJacobianQQNLInertia(double time)
   *  \brief default function to compute the gradient of the external strength compared to the state
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianQQNLInertia(double time);

  /** \fn void computeJacobianVelocityQNLInertia(double time)
   *  \brief default function to compute the gradient of the external strength compared to the velocity
   *  \param double time : the current time
   *  \exception RuntimeException
   */
  virtual void computeJacobianVelocityQNLInertia(double time);

  ////////////////////////////////////////////

  /** \fn void setComputeMassFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute the mass
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeMassFunction(string pluginPath, string functionName);

  /** \fn void setComputeFIntFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute Fint
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeFIntFunction(string pluginPath, string functionName);

  /** \fn void setComputeFExtFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute Fext
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception ICDLL_CSharedLibraryException
   */
  void setComputeFExtFunction(string pluginPath, string functionName);

  /** \fn void setComputeQNLInertiaFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute the inertia
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosCSharedLibraryException
   */
  void setComputeQNLInertiaFunction(string pluginPath, string functionName);

  /** \fn void setComputeJacobianQFIntFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute the gradient of the internal strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQFIntFunction(string pluginPath, string functionName);

  /** \fn void setComputeJacobianVelocityFIntFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute the internal strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityFIntFunction(string pluginPath, string functionName);

  /** \fn void setComputeJacobianQQNLInertiaFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute the gradient of the the external strength compared to the state
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianQQNLInertiaFunction(string pluginPath, string functionName);

  /** \fn void setComputeJacobianVelocityQNLInertiaFunction(string pluginPath, string functionName)
   *  \brief allow to set a specified function to compute the external strength compared to the velocity
   *  \param string : the complete path to the plugin
   *  \param string : the name of the function to use in this plugin
   *  \exception SiconosSharedLibraryException
   */
  void setComputeJacobianVelocityQNLInertiaFunction(string pluginPath, string functionName);

  ////////////////////////////////////////

  /** \fn void createDynamicalSystem(DSXML * nsdsXML, int number, int ndof,
          SiconosVector* q0, SiconosVector* velocity0,
          string fInt, string fExt,
          string jacobianQFInt, string jacobianVelocityFInt,
          string jacobianQQNLInertia, string jacobianVelocityQNLInertia,
          NSDS * nsds)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial coordinates of this DynamicalSystem
   *  \param SiconosVector* : the initial velocity of this DynamicalSystem
   *  \param NSDS * : The NSDS which contains this DynamicalSystem
   *  \param string : the indiaction needed to locate and use the fInt plugin
   *  \param string : the indiaction needed to locate and use the fExt plugin
   *  \param string : the indiaction needed to locate and use the jacobianQFInt plugin
   *  \param string : the indiaction needed to locate and use the jacobianVelocityFInt plugin
   *  \param string : the indiaction needed to locate and use the jacobianQQNLInertia plugin
   *  \param string : the indiaction needed to locate and use the jacobianVelocityQNLInertia plugin
   *  \param NSDS * : The NSDS which contains this DynamicalSystem
   *  \exception RuntimeException
   */
  void createDynamicalSystem(DSXML * dsXML, int number = -1, int ndof = -1,
                             SiconosVector* q0 = NULL, SiconosVector* velocity0 = NULL,
                             string mass = "BasicPlugin:computeMass",
                             string fInt = "BasicPlugin:computeFInt", string fExt = "BasicPlugin:computeFExt",
                             string jacobianQFInt = "BasicPlugin:computeJacobianQFInt",
                             string jacobianVelocityFInt = "BasicPlugin:computeJacobianVelocityFInt",
                             string jacobianQQNLInertia = "BasicPlugin:computeJacobianQQNLInertia",
                             string jacobianVelocityQNLInertia = "BasicPlugin:computeJacobianVelocityQNLInertia",
                             string QNLlInertia = "BasicPlugin:computeQNLInertia");
  //,NSDS * nsds = NULL);

  /** \fn LagrangianNLDS* convert (DynamicalSystem* ds)
   *  \brief encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param DynamicalSystem* : the system which must be converted
   * \return a pointer on the system if it is of the right type, NULL otherwise
   */
  static LagrangianNLDS* convert(DynamicalSystem* ds);

protected:

  /** \fn void init()
   *  \brief initialise value of a Lagrangian NLDS
   */
  virtual void init();


  /** number of degrees of freedom of the system */
  int ndof;
  /** mass of the system */
  SiconosMatrix mass;
  /** coordinates of the system */
  SimpleVector q;
  /** initial coordinates of the system */
  SimpleVector q0;
  /** free coordinate */
  SimpleVector qFree;


  /** memory of previous coordinates of the system */
  SiconosMemory qMemory;
  /** velocity of the system */
  SimpleVector velocity;
  /** initial velocity of the system */
  SimpleVector velocity0;
  /** memory of previous velocity of the system */
  SiconosMemory velocityMemory;
  /** free Velocity */
  SimpleVector velocityFree;
  /** Reaction due to the non smooth law */
  SimpleVector p;


  /** internal strength of the system */
  SimpleVector fInt;
  /** external strength of the system */
  SimpleVector fExt;

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
  SimpleVector QNLInertia;

  /** jacobian/coordinates of internal strength */
  SiconosMatrix jacobianQFInt;
  /** jacobian/velocity of internal strength */
  SiconosMatrix jacobianVelocityFInt;
  /** jacobian/coordinates of inertia */
  SiconosMatrix jacobianQQNLInertia;
  /** jacobian/velocity of inertie */
  SiconosMatrix jacobianVelocityQNLInertia;

  //  /**  */
  //  LagrangianNLDSXML *lnldsxml;

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
   *  \param double* time : the current time
   *  \param double* qPtr : the pointer to the first element of the vector q
   *  \param double* fExtPtr : the pointer to the first element of the vector FInt (in-out parameter)
   */
  void (*computeFExtPtr)(int* sizeOfq, double* time, double* qPtr, double* fExtPtr);

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

#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.h"
#include "LagrangianLinearTIDSXML.h"

#include "SiconosMatrix.h"
//#include "SiconosVector.h"
#include "NewSiconosVector.h"

#include <iostream>
#include <vector>

using namespace std;


class LagrangianLinearTIDSXML;

/** \class LagrangianLinearTIDS
 *  \brief class of Lagrangian invariant time systems, inherited of LagrangianDS
*  \author SICONOS Development Team - copyright INRIA
 *  \version 1.0
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianLinearTIDS  allows to define  and compute a generic ndof-dimensional
 * Lagrangian Linear Time Invariant Dynamical System of the form :
 * \f[
 * M \ddot q + C \dot q + K q =  F_{Ext}( q , t) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,  the variable p stored the impulse and not the force.
 *    -  \f$ M \in  R^{ndof \times ndof} \f$ is Mass matrix stored in the SiconosMatrix mass.
 *    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix stored in the SiconosMatrix K.
 *    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix stored in the SiconosMatrix C.
 *
 *
 *
 *
  *    -  \f$ F_{Int}(\dot q , q , t)\f$ the internal forces stored in the SiconosVector fExt.
 *    -  \f$ F_{Ext}( q , t) \f$ the external forces stored in the SiconosVector fInt.
 *
 *
 * One word on the initial condition.
 *
 * One word on the bilateral constraint
 *
 *
 * The master Class LagrangianDS is specified as follows :
 *    -  \f$ M(q) = M q \f$
 *    -  \f$ Q(\dot q, q) = 0 \f$
 *    -  \f$ F_{Int}(\dot q , q , t) = -C \dot q - K q \f$
 *
 *
 *
 * As for the master Class LagrangianDS, the state of the master class DynamicalSystem is defined by \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$ and then \f$ n= 2 ndof \f$ and the VectorField
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
 * \todo Automatically, specify the function of LagrangianDS such as
 *          Mass, QNL Inertia , FInt = K q +c velocity,
 */
class LagrangianLinearTIDS : public LagrangianDS
{
public:

  /** \fn LagrangianLinearTIDS()
   *  \brief default constructor
   */
  LagrangianLinearTIDS();

  /** \fn LagrangianLinearTIDS(DSXML*)
   *  \brief constructor with XML object of the LagrangianLinearTIDS
   *  \param DSXML* : the XML object corresponding
   */
  LagrangianLinearTIDS(DSXML*);

  ~LagrangianLinearTIDS();

  // getter/setter
  /** \fn SiconosMatrix getK(void)
   *  \brief allow to get the SiconosMatrix K
   *  \return the SiconosMatrix K
   */
  inline SiconosMatrix getK(void) const
  {
    return this->K;
  };

  /** \fn SiconosMatrix getC(void)
   *  \brief allow to get the SiconosMatrix C
   *  \return the SiconosMatrix C
   */
  inline SiconosMatrix getC(void) const
  {
    return this->C;
  };

  /** \fn SiconosMatrix* getKPtr(void)
   *  \brief allow to get the SiconosMatrix K
   *  \return the SiconosMatrix* K
   */
  SiconosMatrix* getKPtr(void);

  /** \fn SiconosMatrix* getCPtr(void)
   *  \brief allow to get the SiconosMatrix C
   *  \return the SiconosMatrix* C
   */
  SiconosMatrix* getCPtr(void);

  /** \fn void setK(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix K
   *  \param the SiconosMatrix to set K
   */
  inline void setK(const SiconosMatrix &K)
  {
    this->K = K;
  };

  /** \fn void setC(SiconosMatrix)
   *  \brief allow to set the SiconosMatrix C
   *  \param the SiconosMatrix to set C
   */
  inline void setC(const SiconosMatrix &C)
  {
    this->C = C;
  };


  //////////////////////////

  /** \fn void init()
   *  \brief initialise value of a Lagrangian TIDS
   */
  virtual void init();

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief print the data to the screen
   */
  void display() const;


  /** \fn void createDynamicalSystem(DSXML * dsXML, int number, int ndof,
        SiconosVector* q0, SiconosVector* velocity0, SiconosMatrix* mass,
        string fExt,SiconosMatrix* K, SiconosMatrix* C)
   *  \brief allows to create the DynamicalSystem with an xml file, or the needed data
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \param int : the number for this DynamicalSystem
   *  \param int : the dimension of this DynamicalSystem
   *  \param SiconosVector* : the initial coordinates of this DynamicalSystem
   *  \param SiconosVector* : the initial velocity of this DynamicalSystem
   *  \param SiconosMatrix* : the mass of this DynamicalSystem
   *  \param string : the indiaction needed to locate and use the fExt plugin
   *  \param SiconosMatrix* : the matrix K of this DynamicalSystem
   *  \param SiconosMatrix* : the matrix C of this DynamicalSystem
  //  *  \param NSDS * : The NSDS which contains this DynamicalSystem
   *  \exception RuntimeException
   */
  void createDynamicalSystem(DSXML * dsXML, int number = -1, int ndof = -1,
                             SiconosVector* q0 = NULL, SiconosVector* velocity0 = NULL,
                             SiconosMatrix* mass = NULL,
                             string fExt = "BasicPlugin:computeFExt",
                             //        string jacobianQFInt="BasicPlugin:computeJacobianQFInt",
                             //        string jacobianVelocityFInt="BasicPlugin:computeJacobianVelocityFInt",
                             //        string jacobianQQNLInertia="BasicPlugin:computeJacobianQQNLInertia",
                             //        string jacobianVelocityQNLInertia="BasicPlugin:computeJacobianVelocityQNLInertia",
                             //        string QNLlInertia="BasicPlugin:computeQNLInertia",
                             SiconosMatrix* K = NULL, SiconosMatrix* C = NULL); //, NSDS * nsds = NULL);

  static LagrangianLinearTIDS* convert(DynamicalSystem* ds);

protected :
  /** \fn void fillDSWithDSXML()
   *  \brief overload of the function for a LagrangianLinearTIDS
   *  \exception RuntimeException
   */
  void fillDSWithDSXML();


private:
  /** specific matrix for a LagrangianLinearTIDS */
  SiconosMatrix K;

  /** specific matrix for a LagrangianLinearTIDS */
  SiconosMatrix C;

  /*
   * fake of plugin ^^
   */
  friend void LTIDSComputeFExt(int* sizeOfq, double* time, double* qPtr, double* fExt);
};

#endif // LAGRANGIANTIDS_H

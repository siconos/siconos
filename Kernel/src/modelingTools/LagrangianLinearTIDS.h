#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.h"
#include "LagrangianLinearTIDSXML.h"

#include "SiconosMatrix.h"
#include "NewSiconosVector.h"
#include "check.h"

#include <iostream>
#include <vector>

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
 * M \ddot q + C \dot q + K q =  F_{Ext}(t) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,  the variable p contains the impulse and not the force.
 *    -  \f$ M \in  R^{ndof \times ndof} \f$ is Mass matrix stored in the SiconosMatrix mass.
 *    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix stored in the SiconosMatrix K.
 *    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix stored in the SiconosMatrix C.
 *
 *
 *
 *
 *    -  \f$ F_{Int}(\dot q , q , t)\f$ the internal forces stored in the SimpleVector fExt.
 *    -  \f$ F_{Ext}(t) \f$ the external forces stored in the SimpleVector fInt.
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

  /** \fn LagrangianLinearTIDS(DSXML * dsXML)
   *  \brief constructor from an xml file
   *  \param DSXML * : the XML object for this DynamicalSystem
   *  \exception RuntimeException
   */
  LagrangianLinearTIDS(DSXML * dsXML);

  /** \fn LagrangianLinearTIDS(int number, int ndof,
      SiconosVector* q0, SiconosVector* velocity0, SiconosMatrix* mass,
      string fExt,SiconosMatrix* K, SiconosMatrix* C)
      *  \brief constructor from a minimum set of data
      *  \param int : the number for this DynamicalSystem
      *  \param int : dimension of this DynamicalSystem
      *  \param SiconosVector* : initial coordinates of this DynamicalSystem
      *  \param SiconosVector* : initial velocity of this DynamicalSystem
      *  \param SiconosMatrix* : mass of this DynamicalSystem
      *  \param string : fExt plugin name and location
      *  \param SiconosMatrix* : matrix K of this DynamicalSystem
      *  \param SiconosMatrix* : matrix C of this DynamicalSystem
      *  \exception RuntimeException
      */
  LagrangianLinearTIDS(int number, int ndof,
                       SiconosVector* q0, SiconosVector* velocity0,
                       SiconosMatrix* mass,
                       std::string  fExt,
                       SiconosMatrix* K, SiconosMatrix* C);

  // destructor
  ~LagrangianLinearTIDS();

  // getter/setter
  /** \fn SiconosMatrix getK()
   *  \brief allow to get the SiconosMatrix K
   *  \return the SiconosMatrix K
   */

  // --- GETTERS AND SETTERS ---

  // -- K --
  /** \fn  const SiconosMatrix getK() const
   *  \brief get the value of K
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getK() const
  {
    return *K;
  }

  /** \fn SiconosMatrix* getKPtr() const
   *  \brief get K
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getKPtr() const
  {
    return K;
  }

  /** \fn void setK (const SiconosMatrix& newValue)
   *  \brief set the value of K to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setK(const SiconosMatrix& newValue)
  {
    *K = newValue;
  }

  /** \fn void setKPtr(SiconosMatrix* newPtr)
   *  \brief set K to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setKPtr(SiconosMatrix *newPtr);

  // -- C --
  /** \fn  const SiconosMatrix getC() const
   *  \brief get the value of C
   *  \return SiconosMatrix
   */
  inline const SiconosMatrix getC() const
  {
    return *C;
  }

  /** \fn SiconosMatrix* getCPtr() const
   *  \brief get C
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getCPtr() const
  {
    return C;
  }

  /** \fn void setC (const SiconosMatrix& newValue)
   *  \brief set the value of C to newValue
   *  \param SiconosMatrix newValue
   */
  inline void setC(const SiconosMatrix& newValue)
  {
    *C = newValue;
  }

  /** \fn void setCPtr(SiconosMatrix* newPtr)
   *  \brief set C to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setCPtr(SiconosMatrix *newPtr) ;

  // --- Miscellaneous ---

  /** \fn void saveDSToXML()
   *  \brief copy the data of the DS to the XML tree
   *  \exception RuntimeException
   */
  void saveDSToXML();

  /** \fn void display()
   *  \brief print the data onto the screen
   */
  void display() const;

  static LagrangianLinearTIDS* convert(DynamicalSystem* ds);

  double dsConvergenceIndicator()
  {
    return 0.0;
  }

private:

  /** \fn LagrangianLinearTIDS()
   *  \brief default constructor
   */
  LagrangianLinearTIDS();

  /** specific matrix for a LagrangianLinearTIDS */
  SiconosMatrix *K;

  /** specific matrix for a LagrangianLinearTIDS */
  SiconosMatrix *C;

  /** Flags to know if pointers have been allocated inside constructors or not */
  bool isKAllocatedIn;
  bool isCAllocatedIn;
};

#endif // LAGRANGIANTIDS_H

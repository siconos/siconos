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

/*! \file LagrangianLinearTIDS.hpp

*/
#ifndef LAGRANGIANTIDS_H
#define LAGRANGIANTIDS_H

#include "LagrangianDS.hpp"

class LagrangianDS;

/** Lagrangian Linear Systems with time invariant coefficients - Derived from LagrangianDS
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 *
 * The class LagrangianLinearTIDS  allows to define  and compute a generic ndof-dimensional
 * Lagrangian Linear Time Invariant Dynamical System of the form :
 * where
 * \f[
 * M \ddot q + C \dot q + K q =  F_{Ext}(t,z) + p,
 * \f]
 * where
 *    - \f$q \in R^{ndof} \f$ is the set of the generalized coordinates,
 *    - \f$ \dot q  \in R^{ndof} \f$  the velocity, i. e. the time derivative of the  generalized coordinates.
 *    - \f$ \ddot q  \in R^{ndof} \f$  the acceleration, i. e. the second time derivative of the  generalized coordinates.
 *    - \f$ p  \in R^{ndof} \f$  the forces due to the Non Smooth Interaction. In particular case of Non Smooth evolution,
 *       the variable p contains the impulse and not the force.
 *    -  \f$ M \in  R^{ndof \times ndof} \f$ is Mass matrix stored in the SiconosMatrix mass.
 *    -  \f$ K \in  R^{ndof \times ndof} \f$ is the stiffness matrix stored in the SiconosMatrix K.
 *    -  \f$ C \in  R^{ndof \times ndof} \f$ is the viscosity matrix stored in the SiconosMatrix C.
 *    -  \f$ z \in R^{zSize}\f$ is a vector of arbitrary algebraic variables, some sort of discret state.
 *
 *
 *
 * Links with first order DynamicalSystem top-class are:
 *
 * \f$ n= 2 ndof \f$
 * \f$ x = \left[\begin{array}{c}q \\ \dot q\end{array}\right]\f$
 *
 * The rhs is given by:
 * \f[
 * rhs(x,t,z) = \left[\begin{array}{c}
 *  \dot q  \\
 * \ddot q = M^{-1}(q)\left[F_{Ext}( q , t, z) - C \dot q - K q   \right]\\
 * \end{array}\right]
 * \f]
 * Its jacobian is:
 * \f[
 * \nabla_{x}rhs(x,t) = \left[\begin{array}{cc}
 *  0  & I \\
 * -M^{-1}K & -M^{-1}C) \\
 * \end{array}\right]
 * \f]
 *  The input due to the non smooth law is:
 * \f[
 * r = \left[\begin{array}{c}0 \\ p \end{array}\right]
 * \f]
 */
class LagrangianLinearTIDS : public LagrangianDS
{
protected:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(LagrangianLinearTIDS);


  /** default constructor */
  LagrangianLinearTIDS() : LagrangianDS() {};

  /** specific matrix for a LagrangianLinearTIDS */
  SP::SiconosMatrix _K;

  /** specific matrix for a LagrangianLinearTIDS */
  SP::SiconosMatrix _C;

public:

  /** constructor from an xml file
   *  \param DynamicalSystemXML * : the XML object for this DynamicalSystem
   */
  LagrangianLinearTIDS(SP::DynamicalSystemXML);

  /** constructor from a set of data
   *  \param SiconosVector: initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param SiconosMatrix : mass matrix of this DynamicalSystem
   *  \param SiconosMatrix : matrix K of this DynamicalSystem
   *  \param SiconosMatrix : matrix C of this DynamicalSystem
   */
  //  LagrangianLinearTIDS(const SiconosVector&, const SiconosVector&, const SiconosMatrix&, const SiconosMatrix&, const SiconosMatrix&);
  LagrangianLinearTIDS(SP::SiconosVector, SP::SiconosVector, SP::SiconosMatrix, SP::SiconosMatrix, SP::SiconosMatrix);

  /** constructor from a set of data
   *  \param SiconosVector: initial coordinates of this DynamicalSystem
   *  \param SiconosVector : initial velocity of this DynamicalSystem
   *  \param SiconosMatrix : mass matrix of this DynamicalSystem
   */
  LagrangianLinearTIDS(SP::SiconosVector, SP::SiconosVector, SP::SiconosMatrix);
  //  LagrangianLinearTIDS(const SiconosVector&, const SiconosVector&, const SiconosMatrix&);

  /** destructor */
  ~LagrangianLinearTIDS();

  /** check that the system is complete (ie all required data are well set)
   * \return a bool
   */
  bool checkDynamicalSystem();

  /**
   * return true if the Dynamical system is linear.
   */
  virtual bool isLinear()
  {
    return true;
  }



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
   *  \param time of initialisation, default value = 0
   *  \param the size of the memory, default size = 1.
   */
  //void initializeNonSmoothInput(double = 0, unsigned int = 1) ;


  // --- GETTERS AND SETTERS ---

  // -- K --
  /** get the value of K
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getK() const
  {
    return *_K;
  }

  /** get K
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix K() const
  {
    return _K;
  }

  /** set the value of K to newValue
   *  \param SiconosMatrix newValue
   */
  void setK(const SiconosMatrix&);

  /** set K to pointer newPtr
   *  \param SP::SiconosMatrix  newPtr
   */
  void setKPtr(SP::SiconosMatrix newPtr);

  // -- C --
  /** get the value of C
   *  \return SimpleMatrix

  inline const SimpleMatrix getC() const { return *C; }
  */
  /** get C
   *  \return pointer on a SiconosMatrix
   */
  inline SP::SiconosMatrix C() const
  {
    return _C;
  }

  /** set the value of C to newValue
   *  \param SiconosMatrix newValue

  void setC(const SiconosMatrix&);
  */

  /** set C to pointer newPtr
   *  \param SP::SiconosMatrix  newPtr
   */
  void setCPtr(SP::SiconosMatrix newPtr) ;

  /** Default function to the right-hand side term
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  void computeRhs(double, bool  = false);

  /** Default function to jacobian of the right-hand side term according to x
   *  \param double time : current time
   *  \param bool isDSup : flag to avoid recomputation of operators
   *  \exception RuntimeException
   */
  void computeJacobianRhsx(double, bool  = false);

  // --- Miscellaneous ---

  /** copy the data of the DS into the XML tree
   */
  void saveSpecificDataToXML();

  /** print the data onto the screen
   */
  void display() const;

  static LagrangianLinearTIDS* convert(DynamicalSystem* ds);

  /** overload LagrangianDS corresponding function
   * \return a double, always zero.
   */
  double dsConvergenceIndicator()
  {
    return 0.0;
  }

  ACCEPT_STD_VISITORS();

};

TYPEDEF_SPTR(LagrangianLinearTIDS)

#endif // LAGRANGIANTIDS_H

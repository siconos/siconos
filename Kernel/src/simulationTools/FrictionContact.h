/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
/*! \file
  Fricton-Contact Non-Smooth Problem
*/
#ifndef FrictionContact_H
#define FrictionContact_H

#include "OneStepNSProblem.h"
#include "SimpleVector.h"
#include "SimpleMatrix.h"

/** Formalization and Resolution of a Friction-Contact Problem
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) Dec 15, 2005
 *
 *  WARNING !!! This is a virtual class -> derived class = FrictionContact2D and 3D !!!
 *
 * This class is devoted to the formalization and the resolution of
 * friction contact problems defined by :
 * \f[
 * velocity =  q + M reaction
 * \f]
 * \f[
 * velocity \geq 0, reaction \geq 0,  reaction^{T} velocity =0
 * \f]
 * and a Coulomb friction law.
 *
 * With:
 *    - \f$velocity \in R^{n} \f$  and \f$reaction \in R^{n} \f$ are the unknown,
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *
 * The present formulation corresponds to pfc2D and 3D of Numerics package.
 *
 * \todo Correct the computation of M with a correct concatenation process
 */
class FrictionContact : public OneStepNSProblem
{
protected:

  /** contains the vector velocity of a FrictionContact system */
  SimpleVector *velocity;

  /** contains the vector reaction of a FrictionContact system */
  SimpleVector *reaction;

  /** contains the matrix M of a FrictionContact system */
  SiconosMatrix *M;

  /** contains the vector q of a FrictionContact system */
  SimpleVector *q;

  /** contains the vector q of a FrictionContact system */
  SiconosVector *mu;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isVelocityAllocatedIn;
  bool isReactionAllocatedIn;
  bool isMAllocatedIn;
  bool isQAllocatedIn;

  /** Specific structure required when a (Numerics) solver block is used */
  SparseBlockStructuredMatrix *Mspbl;

private:

  /** default constructor (private)
   */
  FrictionContact(const std::string pbType = "FrictionContact2D");

public:

  /** xml constructor
   *  \param string: problem type
   *  \param OneStepNSProblemXML* : the XML linked-object
   *  \param Simulation *: the simulation that owns the problem
   */
  FrictionContact(const std::string, OneStepNSProblemXML*, Simulation*);

  /** constructor from data
   *  \param string: problem type
   *  \param Simulation *: the simulation that owns this problem
   *  \param string: id
   */
  FrictionContact(const std::string pbType, Simulation *, const std::string);

  /** constructor from data
   *  \param string: problem type
   *  \param Solver* : pointer to object that contains solver algorithm and formulation
   *  \param Simulation *: the simulation that owns this problem
   *  \param string: id
   */
  FrictionContact(const std::string pbType, Solver*, Simulation *, const std::string);

  /** destructor
   */
  virtual ~FrictionContact();

  // GETTERS/SETTERS

  // --- Velocity ---
  /** get the value of velocity, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getVelocity() const
  {
    return *velocity;
  }

  /** get velocity, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getVelocityPtr() const
  {
    return velocity;
  }

  /** set the value of velocity to newValue
   *  \param SimpleVector newValue
   */
  void setVelocity(const SimpleVector&);

  /** set velocity to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setVelocityPtr(SimpleVector*);

  // --- Reaction ---
  /** get the value of reaction, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getReaction() const
  {
    return *reaction;
  }

  /** get reaction, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getReactionPtr() const
  {
    return reaction;
  }

  /** set the value of reaction to newValue
   *  \param SimpleVector newValue
   */
  void setReaction(const SimpleVector&);

  /** set reaction to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setReactionPtr(SimpleVector*) ;

  // --- M ---

  /** get the value of M
   *  \return SimpleMatrix
   */
  inline const SimpleMatrix getM() const
  {
    return *M;
  }

  /** get M
   *  \return pointer on a SiconosMatrix
   */
  inline SiconosMatrix* getMPtr() const
  {
    return M;
  }

  /** set the value of M to newValue
   *  \param SiconosMatrix newValue
   */
  void setM(const SiconosMatrix&);

  /** set M to pointer newPtr
   *  \param SiconosMatrix * newPtr
   */
  void setMPtr(SiconosMatrix *);

  /** get the structure used to save M as a list of blocks
   *  \return a SparseBlockStructuredMatrix
   */
  inline SparseBlockStructuredMatrix* getMspblPtr() const
  {
    return Mspbl;
  }

  // --- Q ---
  /** get the value of q, the initial state of the DynamicalSystem
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const SimpleVector getQ() const
  {
    return *q;
  }

  /** get q, the initial state of the DynamicalSystem
   *  \return pointer on a SimpleVector
   */
  inline SimpleVector* getQPtr() const
  {
    return q;
  }

  /** set the value of q to newValue
   *  \param SimpleVector newValue
   */
  void setQ(const SimpleVector&);

  /** set q to pointer newPtr
   *  \param SimpleVector * newPtr
   */
  void setQPtr(SimpleVector*);

  // --- Mu ---
  /** get the vector mu, list of the friction coefficients
   *  \return SimpleVector
   */
  inline const SimpleVector getMu() const
  {
    return *mu;
  }

  /** get a pointer to mu, the list of the friction coefficients
   *  \return pointer on a SimpleVector
   */
  inline SiconosVector* getMuPtr() const
  {
    return mu;
  }

  /** get the value of the component number i of mu, the vector of the friction coefficients
   *  \return SimpleVector
   *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
   */
  inline const double getMu(unsigned int i) const
  {
    return (*mu)(i);
  }

  // --- Others functions ---

  /** initialize the FrictionContact problem(compute topology ...)
   */
  void initialize();

  /** computes extra diagonal block-matrix that corresponds to UR1 and UR2
   *  Move this to Unitary Relation class?
   *  \param a pointer to UnitaryRelation
   *  \param a pointer to UnitaryRelation
   */
  void computeBlock(UnitaryRelation*, UnitaryRelation*);

  /** built matrix M using already computed blocks
   */
  void assembleM();

  /** compute vector q
   *  \param double : current time
   */
  void computeQ(const double time);

  /** pre-treatment for LCP
   *  \param double : current time
   *  \return void
   */
  void preCompute(const double time);

  /** Compute the unknown reaction and velocity and update the Interaction (y and lambda )
   *  \param double : current time
   *  \return int, information about the solver convergence.
   */
  virtual int compute(const double time) = 0;

  /** post-treatment for LCP
   */
  void postCompute();

  /** copy the data of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveNSProblemToXML();

  /** copy the matrix M of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveMToXML();

  /** copy the vector q of the OneStepNSProblem to the XML tree
   *  \exception RuntimeException
   */
  void saveQToXML();

  /** print the data to the screen
   */
  void display() const;

  /** encapsulates an operation of dynamic casting. Needed by Python interface.
   *  \param OneStepNSProblem* : the one step problem which must be converted
   * \return a pointer on the problem if it is of the right type, NULL otherwise
   */
  static FrictionContact* convert(OneStepNSProblem* osnsp);
};

#endif // FrictionContact_H

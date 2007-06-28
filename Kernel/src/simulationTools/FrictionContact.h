/* Siconos-Kernel version 2.1.0, Copyright INRIA 2005-2006.
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
 *  \version 2.1.0.
 *  \date (Creation) Dec 15, 2005
 *
 *  WARNING !!! This is a virtual class -> derived class = FrictionContact2D and 3D !!!
 *
 * This class is devoted to the formalization and the resolution of
 * friction contact problems defined by :
 * \f[
 * w =  q + M z
 * \f]
 * \f[
 * w \geq 0, z \geq 0,  z^{T} w =0
 * \f]
 * and a Coulomb friction law.
 *
 * With:
 *    - \f$w \in R^{n} \f$  and \f$z \in R^{n} \f$ are the unknown,
 *    - \f$M \in R^{n \times n } \f$  and \f$q \in R^{n} \f$
 *
 * The present formulation corresponds to pfc2D and 3D of Numerics package.
 *
 * \todo Correct the computation of M with a correct concatenation process
 */
class FrictionContact : public OneStepNSProblem
{
protected:

  /** contains the vector w of a FrictionContact system */
  SimpleVector *w;

  /** contains the vector z of a FrictionContact system */
  SimpleVector *z;

  /** contains the matrix M of a FrictionContact system */
  SiconosMatrix *M;

  /** contains the vector q of a FrictionContact system */
  SimpleVector *q;

  /** Flags to check wheter pointers were allocated in class constructors or not */
  bool isWAllocatedIn;
  bool isZAllocatedIn;
  bool isMAllocatedIn;
  bool isQAllocatedIn;

public:

  /** default constructor
  */
  FrictionContact(const std::string pbType = "FrictionContact2D");

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

  // --- W ---
  /** get the value of w, the initial state of the DynamicalSystem
  *  \return SimpleVector
  *  \warning: SiconosVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getW() const
  {
    return *w;
  }

  /** get w, the initial state of the DynamicalSystem
  *  \return pointer on a SimpleVector
  */
  inline SimpleVector* getWPtr() const
  {
    return w;
  }

  /** set the value of w to newValue
  *  \param SimpleVector newValue
  */
  void setW(const SimpleVector&);

  /** set w to pointer newPtr
  *  \param SimpleVector * newPtr
  */
  void setWPtr(SimpleVector*);

  // --- Z ---
  /** get the value of z, the initial state of the DynamicalSystem
  *  \return SimpleVector
  *  \warning: SimpleVector is an abstract class => can not be an lvalue => return SimpleVector
  */
  inline const SimpleVector getZ() const
  {
    return *z;
  }

  /** get z, the initial state of the DynamicalSystem
  *  \return pointer on a SimpleVector
  */
  inline SimpleVector* getZPtr() const
  {
    return z;
  }

  /** set the value of z to newValue
  *  \param SimpleVector newValue
  */
  void setZ(const SimpleVector&);

  /** set z to pointer newPtr
  *  \param SimpleVector * newPtr
  */
  void setZPtr(SimpleVector*) ;

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

  /** Compute the unknown z and w and update the Interaction (y and lambda )
  *  \param double : current time
  *  \return void
  */
  virtual void compute(const double time) = 0;

  /** post-treatment for LCP
  *  \param 2 pointers to SiconosVector: output of LCP solver
  */
  void postCompute(SiconosVector*, SiconosVector*) ;

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

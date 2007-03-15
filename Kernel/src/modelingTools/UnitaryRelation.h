/* Siconos-Kernel version 2.0.1, Copyright INRIA 2005-2006.
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
/*! \file UnitaryRelation.h

*/
#ifndef UNITARYRELATION_H
#define UNITARYRELATION_H

#include "Tools.h"
#include "SiconosAlgebra.h"
#include "DynamicalSystemsSet.h"
#include "Interaction.h"

/** Interface to single relations from Interactions
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.0.1.
 *  \date (Creation) June 06, 2006
 *
 * Remind that each Interaction is composed with one (pointer to) relation and a non smooth law (size nsLawSize). The number of relations
 * is sizeOfInteraction/nsLawSize. A UnitaryRelation is one of the relations of the Interaction. Actually, this class only provides an interface
 * to handle single relations, this for IndexSets used in Topology.
 * Each UnitaryRelation has a pointer to its "mother" Interaction  and methods to compute y, lambda and so on.
 *
 * - Blocks computation: in OneStepNSProblem, some operators/matrices are required to compute the blocks matrices (used for example for the assembly
 * of Mlcp matrix).
 * In the present class, three methods are available to get the required blocks:
 *  getLeftBlockForDS, getRightBlockForDS and getExtraBlock, with the general model for block computation:
 *
 * block = getExtraBlock  +  getLeftBlockForDS * W * getRightBlockForDS
 *
 * Examples:
 *   => LinearTIR, block = D + h*theta*C*W*B  (D != NULL only for blocks on the diagonal of the full-assembled matrix)
 *    and thus getExtraBlock = D, getLeftBlockForDS = C, getRightBlockForDS = B.
 *
 *   => Lagrangian, block = G* W* transpose(G)  (G=H for lagrangian linear)
 *     ie getExtraBlock = NULL, getLeftBlockForDS = getRightBlockForDS = G (transpose is done using matMultTranspose)
 *
 * Moreover, for, G, B, etc ... we only get the part corresponding to a specific DynamicalSystem (which belongs to the UnitaryRelation)
 *
 */
class UnitaryRelation
{

  // === PRIVATE MEMBERS ===

private:

  /** link to Interaction that owns this relation **/
  Interaction * mainInteraction;

  /** relative position of the present relation in the Interaction - For example if the present relation takes place from index 2 to 4 in y vector
   of mainInteraction, thus the relative position is equal to 2. */
  unsigned int relativePosition;

  /** number of the relation, ie the number of the corresponding block vector in the main Interaction.*/
  unsigned int number;

  /** default constructor
  */
  UnitaryRelation();

public:

  // === PUBLIC FUNCTIONS ===

  // === CONSTRUCTORS/DESTRUCTOR ===

  /** copy constructor
  *  \param UnitaryRelation* : the object to copy
  */
  UnitaryRelation(const UnitaryRelation&);

  /** constructor from a pointer to Interaction
  *  \param Interaction * : Interaction object from which a list of relation will be "extracted"
  *  \param unsigned int: gives the relative position of the relation inside the y vector of the interaction
  *  \param unsigned int: gives the number of the block in y vector of the interaction that corresponds to the present unitary relation.
  */
  UnitaryRelation(Interaction* , const unsigned int, const unsigned int);

  /** destructor
  */
  ~UnitaryRelation();

  /** get main interaction of this unitary relation
  *  \return a pointer to Interaction
  */
  inline Interaction* getInteractionPtr()
  {
    return mainInteraction;
  } ;

  /** get relative position of the Unitary Relation
  *  \return an unsigned int
  */
  inline const unsigned int getRelativePosition() const
  {
    return relativePosition;
  } ;

  /** get number of the Unitary Relation
  *  \return an unsigned int
  */
  inline const unsigned int getNumber() const
  {
    return number;
  };

  /** get vector of output derivatives
  *  \return a VectorOfVectors
  */
  //  const std::vector< SimpleVector* > getY() const;
  const VectorOfVectors getY() const;

  /** get y[i], derivative number i of output
  *  \return pointer on a SimpleVector
  */
  SiconosVector* getYPtr(const unsigned int) const;

  /** get yOld[i], derivative number i of output
  *  \return pointer on a SimpleVector
  */
  SiconosVector* getYOldPtr(const unsigned int) const;

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  const  VectorOfVectors getLambda() const;

  /** get lambda[i], derivative number i of input
  *  \return pointer on a SimpleVector
  */
  SiconosVector* getLambdaPtr(const unsigned int) const;

  /** get y[i], derivative number i of output, value used to compute indexSets
  *  \return a double
  */
  const double getYRef(const unsigned int) const;

  /** get lambda[i], derivative number i of output, value used to compute indexSets
  *  \return a double
  */
  const double getLambdaRef(const unsigned int) const;

  /** returns the size of the embedded non smooth law
  *  \return an unsigned int
  */
  const unsigned int getNonSmoothLawSize() const;

  /** returns the type of the embedded non smooth law
  *  \return a string
  */
  const std::string getNonSmoothLawType() const;

  /** returns the type of the embedded relation.
  *  \return a string
  */
  const std::string getRelationType() const;

  /** returns the subtype of the embedded relation.
   *  \return a string
   */
  const std::string getRelationSubType() const;

  /** gets an iterator to the first element of the DynamicalSystems set.
   *  \return a DSIterator.
   */
  inline DSIterator dynamicalSystemsBegin()
  {
    return mainInteraction->dynamicalSystemsBegin();
  };

  /** gets an iterator equal to DynamicalSystems.end().
   *  \return a DSIterator.
   */
  inline DSIterator dynamicalSystemsEnd()
  {
    return mainInteraction->dynamicalSystemsEnd();
  };

  /** gets a const iterator to the first element of the  DynamicalSystems set.
   *  \return a ConstDSIterator.
   */
  inline ConstDSIterator dynamicalSystemsBegin() const
  {
    return mainInteraction->dynamicalSystemsBegin();
  };

  /** gets a const iterator equal to DynamicalSystems.end().
   *  \return a ConstDSIterator.
   */
  inline ConstDSIterator dynamicalSystemsEnd() const
  {
    return mainInteraction->dynamicalSystemsEnd();
  };

  /** gets the DynamicalSystemsSet
   *  \return a DynamicalSystemsSet
   */
  DynamicalSystemsSet getDynamicalSystems() const;

  /** gets a pointer to the DynamicalSystemsSet
   *  \return a DynamicalSystemsSet*
   */
  DynamicalSystemsSet * getDynamicalSystemsPtr() ;

  /** gets the matrix used in block computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting block matrix
   *  \param unsigned int (optional): gradient index (see relations for details)
   */
  void getLeftBlockForDS(DynamicalSystem *, SiconosMatrix*, const unsigned int = 0) const;

  /** gets the matrix used in block computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting block matrix
   *  \param unsigned int (optional): gradient index (see relations for details)
   */
  void getRightBlockForDS(DynamicalSystem *, SiconosMatrix*, const unsigned int = 0) const;

  /** gets extra block corresponding to the present UR (see the top of this files for extra block meaning)
   *  \param a pointer to a SiconosMatrix (in-out parameter)
   */
  void getExtraBlock(SiconosMatrix *) const;

  /** computes free value of y
   *  \param a double (current time)
   *  \param unsigned int: derivative order for y
   *  \param string: simulation type
   *  \param pointer to SiconosVector (in-out parameter), result
   */
  void computeEquivalentY(const double, const unsigned int, const std::string, SiconosVector*);
};

#endif // UNITARYRELATION_H

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
/*! \file UnitaryRelation.hpp

*/
#ifndef UNITARYRELATION_H
#define UNITARYRELATION_H

#include "Tools.hpp"
#include "SiconosAlgebra.hpp"
#include "Interaction.hpp"
#include "RelationNamespace.hpp"

#include "SiconosPointers.hpp"

/** Interface to single relations from Interactions
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) June 06, 2006
 *
 * Remind that each Interaction is composed with one (pointer to)
 * relation and a non smooth law (size nsLawSize). The number of
 * relations is sizeOfInteraction/nsLawSize. A UnitaryRelation is one
 * of the relations of the Interaction. Actually, this class only
 * provides an interface to handle single relations, this for
 * IndexSets used in Topology.  Each UnitaryRelation has a pointer to
 * its "mother" Interaction and methods to compute y, lambda and so
 * on.
 *
 * - UnitaryBlocks computation: in OneStepNSProblem, some
 * operators/matrices are required to compute the unitaryBlocks
 * matrices (used for example for the assembly of Mlcp matrix).  In
 * the present class, three methods are available to get the required
 * unitaryBlocks: getLeftUnitaryBlockForDS, getRightUnitaryBlockForDS
 * and getExtraUnitaryBlock, with the general model for unitaryBlock
 * computation:
 *
 * unitaryBlock = getExtraUnitaryBlock  +  getLeftUnitaryBlockForDS * W * getRightUnitaryBlockForDS
 *
 * Examples:
 *
 *   => LinearTIR, unitaryBlock = D + h*theta*C*W*B (D != NULL only
 *      for unitaryBlocks on the diagonal of the full-assembled
 *      matrix) and thus getExtraUnitaryBlock = D,
 *      getLeftUnitaryBlockForDS = C, getRightUnitaryBlockForDS = B.
 *
 *   => Lagrangian, unitaryBlock = G* W* transpose(G) (G=H for
 *      xslagrangian linear) ie getExtraUnitaryBlock = NULL,
 *      getLeftUnitaryBlockForDS = getRightUnitaryBlockForDS = G
 *      (transpose is done using matMultTranspose)
 *
 * Moreover, for, G, B, etc ... we only get the part corresponding to
 * a specific DynamicalSystem (which belongs to the UnitaryRelation)
 *
 */
class UnitaryRelation
{

  // === PRIVATE MEMBERS ===

private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(UnitaryRelation);


  /** link to Interaction that owns this relation **/
  SP::Interaction _mainInteraction;

  // /** relative position of the present relation in the Interaction -
  //  For example if the present relation takes place from index 2 to 4
  //  in y vector of mainInteraction, the relative position is equal to
  //  2. */
  // unsigned int _relativePosition;

  // /** number of the relation, ie the number of the corresponding
  //     unitaryBlock vector in the main Interaction.*/
  // unsigned int _number;

  // /** Absolute position in the "global" vector of constraints (for
  //     example, the one handled by lsodar) */
  // unsigned int _absolutePosition;
  // /** Absolute position in the "global" vector of constraints for the proj formulation. */

  // unsigned int _absolutePositionProj;

  // /** work vector to save pointers to state-related data of the
  //     dynamical systems involved in the UR.*/
  // SP::SiconosVector _workX;
  // SP::SiconosVector _workXq;
  // SP::SiconosVector _workFree;

  // SP::SiconosVector _workYp;

  // /** work vector to save pointers to z data of the dynamical systems
  //     involved in the UR.*/
  // SP::SiconosVector _workZ;



  /** default constructor
   */
  UnitaryRelation() {};

  /** copy constructor: private, no copy cor pass-by value
   */
  UnitaryRelation(const UnitaryRelation&);

public:

  /** constructor from a pointer to Interaction
  *  \param SP::Interaction : Interaction object from which a list of
  *  relation will be "extracted"
  *  \param unsigned int: gives the relative position of the relation
  *  inside the y vector of the interaction
  *  \param unsigned int: gives the number of the unitaryBlock in y
  *  vector of the interaction that corresponds to the present
  *  unitary relation.
  */
  UnitaryRelation(SP::Interaction inter, unsigned int, unsigned int): _mainInteraction(inter)
  {};

  /** destructor
  */
  ~UnitaryRelation() {};

  /** get main interaction of this unitary relation
  *  \return a pointer to Interaction
  */
  inline SP::Interaction interaction() const
  {
    assert(_mainInteraction);
    return _mainInteraction;
  }

  /** get relative position of the Unitary Relation
  *  \return an unsigned int
  */
  inline unsigned int getRelativePosition() const
  {
    //return _relativePosition;
    return 0;
  } ;

  /** get id of the parent interaction
   *  \return a string
   */
  inline const std::string getId() const
  {
    return _mainInteraction->getId();
  };

  /** get number of the Unitary Relation
   *  \return an unsigned int
   */
  inline unsigned int number() const
  {
    //return _number;
    assert(0);
    return 1;
  };

  /** get y[i], derivative number i of output
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector y(unsigned int i) const
  {
    // i is the derivative number.
    return (interaction()->y(i));
  };

  /** get yOld[i], derivative number i of output
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector yOld(unsigned int i) const
  {
    // i is the derivative number.
    return (interaction()->yOld(i));
  };

  /* get y_k[i]
   *    \return pointer on a SiconosVector
   */
  inline SP::SiconosVector y_k(unsigned int i) const
  {
    //i is the derivative number.
    return (interaction()->y_k(i));
  };

  // /* get yMemory[i][j]
  //  *    \return pointer on a SiconosVector
  //  * i is the derivative number.
  //  * j is the depth in time
  //  */
  // SP::SiconosVector yMemory(unsigned int,unsigned int) const;

  /* get yMemory[i]
   *    \return pointer on a SiconosMemory
   * \param unsigned int i is the derivative number.
   */
  inline  SP::SiconosMemory yMemory(unsigned int i) const
  {
    //i is the derivative number.
    return interaction()->yMemory(i);
  };

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  inline const  VectorOfVectors getLambda() const
  {
    // A new object of type VectorOfVectors is created but it handles
    // pointers to BlockVectors, thus there is no copy of the "basic"
    // SiconosVectors.
    return  interaction()->getLambda();
  };

  /** get lambda[i], derivative number i of input
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector lambda(unsigned int i) const
  {
    // i is the derivative number.
    return ((interaction()->lambda(i)));
  };

  /** get y[i], derivative number i of output, value used to compute indexSets
  *  \return a double
  */
  double getYRef(unsigned int i) const
  {
    return ((interaction()->getYRef(i)));
  };

  /** get lambda[i], derivative number i of output, value used to compute indexSets
  *  \return a double
  */
  double getLambdaRef(unsigned int i) const
  {
    return ((interaction()->getLambdaRef(i)));
  };

  /** returns the size of the embedded non smooth law
  *  \return an unsigned int
  */
  inline unsigned int getNonSmoothLawSize() const
  {
    return interaction()->nonSmoothLaw()->size();
  };

  unsigned int absolutePosition()
  {
    return _mainInteraction->absolutePosition();
  };
  void setAbsolutePosition(unsigned int v)
  {
    _mainInteraction->setAbsolutePosition(v);
  };
  unsigned int absolutePositionProj()
  {
    return _mainInteraction->absolutePositionProj();
  };
  void setAbsolutePositionProj(unsigned int v)
  {
    _mainInteraction->setAbsolutePositionProj(v);
  };

  /** temporary visitor to get type
      must be removed */
  struct GetNSLType;
  friend class GetNSLType;

  /** returns the type of the embedded relation.
   */
  inline RELATION::TYPES getRelationType() const
  {
    return interaction()->relation()->getType();
  };

  /** returns the subtype of the embedded relation.
   */
  inline RELATION::SUBTYPES getRelationSubType() const
  {
    return interaction()->relation()->getSubType();
  } ;

  /** To initialize the UR: mainly to set work vectors.
   */
  void initialize(const std::string&)
  {
    RuntimeException::selfThrow("UnitaryRelation::initialize(simulationType) - Obsolete - should not be called");
  };

  /* to set workX content.
   \param a SP::SiconosVector to be inserted into workX
  */
  inline void insertInWorkX(SP::SiconosVector newX)
  {
    assert(_mainInteraction->workX()) ;
    _mainInteraction->workX()->insertPtr(newX);
  };
  /* to set _workFree content.
   \param a SP::SiconosVector to be inserted into workFree
  */
  inline void insertInWorkFree(SP::SiconosVector newX)
  {
    assert(_mainInteraction->workFree()) ;
    _mainInteraction->workFree()->insertPtr(newX);
  };

  /** Get a pointer to workX */
  inline SP::SiconosVector workx()
  {
    return _mainInteraction->workX();
  };
  inline SP::SiconosVector xq()
  {
    return _mainInteraction->workXq();
  };
  inline SP::SiconosVector workFree()
  {
    return _mainInteraction->workFree();
  };

  inline SP::SiconosVector yp()
  {
    return _mainInteraction->yp();
  };


  /** Get a pointer to workZ */
  inline SP::SiconosVector workz()
  {
    return _mainInteraction->workZ();
  };


  /** gets the matrix used in unitaryBlock computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting unitaryBlock matrix
   */
  inline void getLeftUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
  {
    interaction()->getLeftUnitaryBlockForDS(ds, UnitaryBlock);
  };
  /** gets the matrix used in unitaryBlock computation. Used only for the formulation projecting on the constraints.
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting unitaryBlock matrix
   */
  inline void getLeftUnitaryBlockForDSProjectOnConstraints(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
  {
    interaction()->getLeftUnitaryBlockForDSProjectOnConstraints(ds, UnitaryBlock);
  };
  /** gets the matrix used in unitaryBlock computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting unitaryBlock matrix
   */
  inline void getRightUnitaryBlockForDS(SP::DynamicalSystem ds , SP::SiconosMatrix UnitaryBlock) const
  {
    interaction()->getRightUnitaryBlockForDS(ds, UnitaryBlock);
  };

  /** gets extra unitaryBlock corresponding to the present UR (see the
   *  top of this files for extra unitaryBlock meaning)
   * \param a pointer to a SiconosMatrix (in-out parameter)
   */
  inline void getExtraUnitaryBlock(SP::SiconosMatrix UnitaryBlock) const
  {
    interaction()->getExtraUnitaryBlock(UnitaryBlock);
  };

};

TYPEDEF_SPTR(UnitaryRelation);
#endif // UNITARYRELATION_H

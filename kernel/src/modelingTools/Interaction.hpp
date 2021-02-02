/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

/*! \file Interaction.hpp
  \brief Interaction class and related typedef
*/


#ifndef INTERACTION_H
#define INTERACTION_H

#include "RelationNamespace.hpp"
#include "SiconosPointers.hpp"
#include "SiconosVector.hpp"
#include "SiconosMemory.hpp"
#include "SiconosAlgebraTypeDef.hpp"
#include <vector>

/** Description of a non-smooth interaction

    The object Interaction is used to defined a "link" between one or two DynamicalSystem,
    like unilateral constraints and some nonsmooth law (e.g. complementarity).

    It holds two vectors of "local" variables, \f$y\f$ and \f$\lambda\f$
    and their derivatives, which are related to the state variables and the inputs of the
    DynamicalSystem (x,R) through constraints defined in a Relation and completed with
    a NonSmoothLaw involving those variables.

    Remarks:

    - one and only one Relation (access: relation()) per Interaction
    - one and only one NonSmoothLaw (access: nonSmoothLaw()) per Interaction
    - dimension() is the size of the interaction and so the size of vectors y, lambda
      and their derivatives.
    - output: y(i), to get derivative i of y
    - input: lambda(i), to get derivative i of lambda
 
 */
class Interaction : public std::enable_shared_from_this<Interaction >
{
private:
  /* serialization hooks */
  ACCEPT_SERIALIZATION(Interaction);

  /* internal counter used to set interaction number */
  static size_t __count;

  /** Interaction id */
  size_t _number;

  /** Minimum required 'level' for output y
   *  y will be initialized from
   *  y[_lowerLevelForOutput] to y[_upperLevelForOutput]
   */
  unsigned int _lowerLevelForOutput = 0;

  /** Maximum required 'level' for output y
    *  y will be initialized from
    *  y[_lowerLevelForOutput] to y[_upperLevelForOutput]
    */
  unsigned int _upperLevelForOutput = 0;

  /** Minimum required 'level' for input lambda
   *  lambda will be initialized from
   *  lambda[_lowerLevelForIntput] to lambda[_upperLevelForInput]
   */
  unsigned int _lowerLevelForInput = 0;

  /** Maximum required 'level' for input lambda
   *  lambda will be initialized from
   *  lambda[_lowerLevelForIntput] to lambda[_upperLevelForInput]
   */
  unsigned int _upperLevelForInput = 0;

  /** size of the interaction, ie size of y[i] and _lambda[i] */
  unsigned int _interactionSize = 0;

  /** sum of all DS sizes, for DS involved in the interaction */
  unsigned int _sizeOfDS = 0;

  /** Bool to check the number of DS concerned by this interaction
      (1 or 2 indeed)
      True if 2 DS.
      Note FP : usefull in NewtonEuler jacobians computation.
  */
  bool _has2Bodies = false;

  /** relation between constrained variables and states variables
   * vector of output derivatives
   * y[0] is y, y[1] is yDot and so on
   */
  VectorOfVectors _y;

  /** previous value of Newton iteration for y
   * \warning : VA 24/05/2013 this has to be put into the workspace vector
   *   or we have to use the _yMemory storage
   */
  VectorOfVectors _yOld;

  /** value of the previous time-step */
  VectorOfVectors _y_k;

  /** memory of previous coordinates of the system */
  VectorOfMemories _yMemory;

  /** memory of previous coordinates of the system */
  VectorOfMemories _lambdaMemory;

  /** result of the computeInput function */
  VectorOfVectors _lambda;

  /** previous step values for _lambda
   * \warning : VA 24/05/2013 this has to be put into the workspace vector
   *   or we have to use the _yMemory storage
   */
  VectorOfVectors _lambdaOld;

  /** the Non-smooth Law of the interaction*/
  SP::NonSmoothLaw _nslaw;

  /** the type of Relation of the interaction */
  SP::Relation _relation;

  /** pointer links to DS variables needed for computation,
   *  mostly used in Relations (computeOutput and computeInput)
   * and OneStepIntegrator classes. */
  VectorOfBlockVectors _linkToDSVariables;

  VectorOfSMatrices _relationMatrices;

  VectorOfVectors _relationVectors;

  struct _setLevels;
  friend struct Interaction::_setLevels;

  // === PRIVATE FUNCTIONS ===

  /* forbid default constructor, copy and assignment */
  Interaction(const Interaction& inter) = delete;
  Interaction& operator=(const Interaction&) = delete;

protected:
  Interaction() {} /* for serialization only */

private:

  /*! @name DSlink graph property management */
  //@{

  /** update DSlink property content with dynamical systems members, FirstOrder relations case.
      \param DSlink the container of the link to DynamicalSystem attributes
      \param ds1 first ds linked to this Interaction
      \param ds2 second ds linked to this Interaction. ds1 == ds2 is allowed.
   */
  void __initDataFirstOrder(VectorOfBlockVectors& DSlink, DynamicalSystem& ds1, DynamicalSystem& ds2);

  /** Initialize the link with the DynamicalSystem, FirstOrderDS variant
   * \param ds a DynamicalSystem concerned by this Interaction
   * \param DSlink the container of the link to DynamicalSystem attributes
   */
  void __initDSDataFirstOrder(DynamicalSystem& ds, VectorOfBlockVectors& DSlink);

  /** update DSlink property content with dynamical systems members, Lagrangian relations case.
      \param DSlink the container of the link to DynamicalSystem attributes
      \param ds1 first ds linked to this Interaction
      \param ds2 second ds linked to this Interaction. ds1 == ds2 is allowed.
   */
  void __initDataLagrangian(VectorOfBlockVectors& DSlink, DynamicalSystem& ds1, DynamicalSystem& ds2);

  /** Initialize the link with the DynamicalSystem, LagrangianDS variant
   * \param ds a DynamicalSystem concerned by this Interaction
   * \param DSlink the container of the link to DynamicalSystem attributes
   */
  void __initDSDataLagrangian(DynamicalSystem& ds, VectorOfBlockVectors& DSlink);

  /** update DSlink property content with dynamical systems members, NewtonEuler relations case.
      \param DSlink the container of the link to DynamicalSystem attributes
      \param ds1 first ds linked to this Interaction
      \param ds2 second ds linked to this Interaction. ds1 == ds2 is allowed.
  */
  void __initDataNewtonEuler(VectorOfBlockVectors& DSlink, DynamicalSystem& ds1, DynamicalSystem& ds2);

  /** Initialize the link with the DynamicalSystem, NewtonEulerDS variant
   * \param ds a DynamicalSystem concerned by this Interaction
   * \param DSlink the container of the link to DynamicalSystem attributes
   */
  void __initDSDataNewtonEuler(DynamicalSystem& ds, VectorOfBlockVectors& DSlink);
  ///@}

public:

  /** Interaction constructor
      \param NSL pointer object describing the nonsmooth law;
      the interaction size if infered from the size of this law.
      \param rel a pointer object describing the functions used to compute the constraints 
  */
  Interaction(SP::NonSmoothLaw NSL, SP::Relation rel);

  /** destructor  */
  ~Interaction() {};

  /** Update interactions attributes.
      Must be called when levels have been modified.
  */
  void reset();
  
  /** set the links to the DynamicalSystem(s) and allocate the required workspaces
   *  \param interProp the InteractionProperties of this Interaction
      \param ds1 first ds linked to this Interaction (i.e IG->vertex.source)
      \param workV1 work vectors of ds1
      \param ds2 second ds linked to this Interaction (i.e IG->vertex.target) ds1 == ds2 is allowed.
      \param workV2 work vectors of ds2
   */
  //void setDSLinkAndWorkspace(InteractionProperties& interProp, DynamicalSystem& ds1, VectorOfVectors& workV1, DynamicalSystem& ds2, VectorOfVectors& workV2);

  /*! @name DSlink graph property management */
  //@{

  /** set the links  between the interaction and the DynamicalSystem(s) members.
      \param ds1 first ds linked to this Interaction (i.e IG->vertex.source)
      \param ds2 second ds linked to this Interaction (i.e IG->vertex.target) ds1 == ds2 is allowed.
  */
  void initializeLinkToDsVariables(DynamicalSystem& ds1,
                                   DynamicalSystem& ds2);
  ///@}
  
  /** set all lambda to zero */
  void resetAllLambda() ;

  /** set lambda to zero for a given level
   * \param level
   */
  void resetLambda(unsigned int level);

  /** build memories vectors for y and \f$\lambda\f$
   * \param computeResiduY true if interaction should compute extra residu value
   * \param steps number of required memories (depends on the OSI)
  */
  void initializeMemory(unsigned int steps);

  // === GETTERS/SETTERS ===
  /** \return the id of the interaction */
  inline size_t number() const
  {
    return _number;
  }

  /** Set the lower level for output y.
   * \param newVal : an unsigned int
   */
  inline void setLowerLevelForOutput(const unsigned int newVal)
  {
    _lowerLevelForOutput = newVal;
  };

  /** Set the upper level for output y.
   * \param newVal : an unsigned int
   */
  inline void setUpperLevelForOutput(const unsigned int newVal)
  {
    _upperLevelForOutput = newVal;
  };

  /** Get the lower level for output y.
      \return an unsigned int.
   */
  inline unsigned int lowerLevelForOutput()
  {
    return _lowerLevelForOutput;
  };

  /** Get the upper level for output y.
      \return an unsigned int.
   */
  inline unsigned int  upperLevelForOutput()
  {
    return _upperLevelForOutput;
  };

  /** Set the lower level for input Lambda.
   * \param newVal : an unsigned int
   */
  inline void setLowerLevelForInput(const unsigned int newVal)
  {
    _lowerLevelForInput = newVal;
  };

  /** Set the upper level for input Lambda.
   * \param newVal : an unsigned int.
   */
  inline void setUpperLevelForInput(const unsigned int newVal)
  {
    _upperLevelForInput = newVal;
  };


  /** Get the lower level for input Lambda.
      \return an unsigned int.
   */
  inline unsigned int lowerLevelForInput()
  {
    return _lowerLevelForInput ;
  };

  /** Get the upper level for input Lambda.
      \return an unsigned int.
   */
  inline unsigned int upperLevelForInput()
  {
    return _upperLevelForInput;
  };

  /** returns dimension (i.e. nslaw size == y and lambda size) */
  inline unsigned int dimension() const
  {
    return _interactionSize;
  }
  
  /** Get the sum of DS sizes, for DS involved in interaction.
   *  \return an unsigned int
   */
  inline unsigned int getSizeOfDS() const
  {
    return _sizeOfDS;
  }

  /** Set the number of dynamical systems concerned by
      this interaction. Warning FP: this function is supposed
      to be called only during topology->link(inter, ds1, ds2) call.
      \param val : true if two ds, else false
   */
  void setHas2Bodies(bool val) {_has2Bodies = val;}

  /** Check the number of dynamical systems concerned by
      this interaction.
      \return bool : true if two ds, else false
   */
  bool has2Bodies() const {return _has2Bodies;}

  // -- y --

  /** Get y[i], derivative number i of output.
   * \param i : the derivative number.
   *  \return BlockVector
   */
  inline const SiconosVector getCopyOfy(const unsigned int i) const
  {
    assert(_y[i] && "_y[i]");
    return *(_y[i]);
  }

  /** get y[i], derivative number i of output
   * \param i derivative number i of output
   * \return BlockVector
   */
  inline const SiconosVector getCopyOfyOld(const unsigned int i) const
  {
    assert(_yOld[i]);
    return *(_yOld[i]);
  }


  /** get vector of output derivatives
   *  \return a VectorOfVectors
   */
  inline const VectorOfVectors y() const
  {
    return _y;
  }


  /** get y[i], derivative number i of output
   * \param i derivative number i of output
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector y(const unsigned int i) const
  {
    assert(_y[i]);
    return _y[i];
  }

  /** set the output vector y to newVector with copy of the y[i] (ie
      memory allocation)
  *  \param v VectorOfVectors
  */
  void setY(const VectorOfVectors& v);

  /** set the output vector y to newVector with direct pointer
  *  equality for the y[i]
  * \param v VectorOfVectors
  */
  void setYPtr(const VectorOfVectors& v);

  /** set y[i] to newValue
   * \param i derivative number i of output
   * \param v a SiconosVector and an unsigned int
   */
  void setY(const unsigned int i, const SiconosVector& v);

  /** set y[i] to pointer newPtr
   * \param i derivative number i of output
   * \param v a SP::SiconosVector  and an unsigned int
   */
  void setYPtr(const unsigned int i, SP::SiconosVector v);

  // -- yOld --

  /** get vector of output derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getYOld() const
  {
    return _yOld;
  }

  /** get yOld[i], derivative number i of output
   * \param i derivative number i of output
   * \return BlockVector
   */
  inline const SiconosVector getYOld(const unsigned int i) const
  {
    assert(_yOld[i]);
    return *(_yOld[i]);
  }

  /** get yOld[i], derivative number i of output
   * \param i derivative number i of output
   * \return pointer on a SiconosVector
   */
  inline SP::SiconosVector yOld(const unsigned int i) const
  {
    assert(_yOld[i]);
    return _yOld[i];
  }
  /* get y_k[i]
   * \param i derivative number i of output
   * \return pointer on a SiconosVector
   */
  inline SP::SiconosVector y_k(const unsigned int i) const
  {
    assert(_y_k[i]);
    return _y_k[i];
  }

  /** set the output vector yOld to newVector
   * \param  v VectorOfVectors
   */
  void setYOld(const VectorOfVectors& v);

  /** set vector yOld to newVector with direct pointer equality for
   *  the yOld[i]
   * \param  v VectorOfVectors
   */
  void setYOldPtr(const VectorOfVectors& v);

  /** set yOld[i] to newValue
   * \param i derivative number i of output
   * \param v a SiconosVector and an unsigned int
   */
  void setYOld(const unsigned int i, const SiconosVector& v);

  /** set yOld[i] to pointer newPtr
   * \param i derivative number i of output
   * \param v a SP::SiconosVector  and an unsigned int
   */
  void setYOldPtr(const unsigned int i, SP::SiconosVector v);


  /** get all the values of the state vector y stored in memory
   * \param level
   * \return a memory
   */
  inline SiconosMemory& yMemory(unsigned int level)
  {
    return _yMemory[level];
  }

  /** get all the values of the multiplier lambda stored in memory
   * \param level
   * \return a memory
   */
  inline SiconosMemory& lambdaMemory(unsigned int level)
  {
    return _lambdaMemory[level];
  }

  // -- _lambda --

  /** get vector of input derivatives
   *  \return a VectorOfVectors
   */
  inline const VectorOfVectors getLambda() const
  {
    return _lambda;
  }

  /** get _lambda[i], derivative number i of input
   * \param i derivative number i of output
   * \return SiconosVector
   */
  inline const SiconosVector getLambda(const unsigned int i) const
  {
    assert(_lambda[i]);
    return *(_lambda[i]);
  }

  /** get _lambda[i], derivative number i of input
   * \param i derivative number i of output
   * \return pointer on a SiconosVector
   */
  inline SP::SiconosVector lambda(const unsigned int i) const
  {
    assert(_lambda[i]);
    return _lambda[i];
  }

  /** set the input vector _lambda to newVector
   *  \param v VectorOfVectors
   */
  void setLambda(const VectorOfVectors& v);

  /** set vector _lambda to newVector with direct pointer equality for the _lambda[i]
   *  \param v VectorOfVectors
   */
  void setLambdaPtr(const VectorOfVectors& v);

  /** set _lambda[i] to newValue
   * \param i derivative number i of output
   *  \param newValue a SiconosVector 
   */
  void setLambda(const unsigned int i, const SiconosVector& newValue);

  /** set _lambda[i] to pointer newPtr
   * \param i derivative number i of output
   * \param newPtr a SP::SiconosVector
   */
  void setLambdaPtr(const unsigned int i, SP::SiconosVector newPtr);

  // -- _lambdaOld --

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getLambdaOld() const
  {
    return _lambdaOld;
  }

  /** get _lambdaOld[i], derivative number i of input
   * \param i derivative number i of output
   *  \return SiconosVector
   */
  inline const SiconosVector getLambdaOld(const unsigned int i) const
  {
    assert(_lambdaOld[i]);
    return *(_lambdaOld[i]);
  }

  /** get _lambdaOld[i], derivative number i of input
   * \param i derivative number i of output
   *  \return pointer on a SiconosVector
   */
  inline SP::SiconosVector lambdaOld(const unsigned int i) const
  {
    assert(_lambdaOld[i]);
    return _lambdaOld[i];
  }

  /** set the input vector _lambdaOld to newVector
   * \param v VectorOfVectors
   */
  void setLambdaOld(const VectorOfVectors& v);

  /** set vector _lambdaOld to newVector with direct pointer equality for the _lambdaOld
  *  \param v VectorOfVectors
  */
  void setLambdaOldPtr(const VectorOfVectors& v);

  /** set _lambdaOld[i] to newValue
   * \param i derivative number i of output
   * \param v a SiconosVector
   */
  void setLambdaOld(const unsigned int i, const SiconosVector& v );

  /** set _lambdaOld[i] to pointer newPtr
   * \param i derivative number i of output
   * \param newPtr a SP::SiconosVector
   */
  void setLambdaOldPtr(const unsigned int i, SP::SiconosVector newPtr);

  /** get the Relation of this Interaction
   *  \return a pointer on this Relation
   */
  inline SP::Relation relation() const
  {
    return _relation;
  }

  /** get the NonSmoothLaw of this Interaction
   *  \return a pointer on this NonSmoothLaw
   */
  inline SP::NonSmoothLaw nonSmoothLaw() const
  {
    return _nslaw;
  }


  inline VectorOfBlockVectors & linkToDSVariables()
  {
    return _linkToDSVariables;
  };

  inline VectorOfVectors & relationVectors()
  {
    return _relationVectors;
  };

  inline VectorOfSMatrices & relationMatrices()
  {
    return _relationMatrices;
  };
  
  // --- OTHER FUNCTIONS ---

  /** set interaction 'ds-dimension', i.e. sum of all sizes of the dynamical systems linked
   *  by the current interaction. This must be done by topology during call to link(inter, ds, ...).
   * \param s1 int sum of ds sizes
  */
  inline void setDSSizes(unsigned int s1)
  {
    _sizeOfDS = s1;
  }

  /**   put values of y into yOld, the same for _lambda
  */
  void swapInOldVariables();

  /** Must be call to fill _y_k. (after convergence of the Newton iterations)
   */
  void swapInMemory();

  /** print the data to the screen
  */
  void display(bool brief = true) const;

  /** reset the global Interaction counter (for ids)
   *  \return the previous value of count
   */
  static inline int resetCount(int new_count=0)
  {
    int old_count = __count;
    __count = new_count;
    return old_count;
  };
  
  /** Computes output y; depends on the relation type.
   *  \param time current time
   *  \param derivativeNumber number of the derivative to compute,
   *  optional, default = 0.
   */
  void computeOutput(double time, unsigned int derivativeNumber = 0);
  /** Compute input r of all Dynamical Systems involved in the present
   *   Interaction.
   *  \param time current time
   *  \param level order of _lambda used to compute input.
   */
  void computeInput(double time, unsigned int level = 0);

  /** gets the matrix used in interactionBlock computation, (left * W * right), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to one ds.
   * \param pos int, relative position of the beginning of the required block in relation matrix.
   * \param size int, size(0) of the block
   * \param sizeDS int, size(1) of the block
   *  \return InteractionBlock a pointer to SiconosMatrix (in-out parameter): the resulting interactionBlock matrix
   */
  SP::SiconosMatrix getLeftInteractionBlockForDS(unsigned int pos, unsigned int size,  unsigned int sizeDS) const;

  /** gets the matrix used in interactionBlock computation. Used only for the formulation projecting on the constraints.
   *         We get only the part corresponding to ds.
   *  \param pos int, relative position of the beginning of the required block in relation matrix.
   *  \param InteractionBlock a pointer to SiconosMatrix (in-out parameter): the resulting interactionBlock matrix
   */
  void getLeftInteractionBlockForDSProjectOnConstraints(unsigned int pos, SP::SiconosMatrix InteractionBlock) const;

  /** gets the matrix used in interactionBlock computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   * \param pos int, relative position of the beginning of the required block in relation matrix.
   * \param sizeDS int, size(0) of the block
   * \param size int, size(1) of the block
   * \return InteractionBlock a pointer to SiconosMatrix (in-out parameter): the resulting interactionBlock matrix
   */
  SP::SiconosMatrix  getRightInteractionBlockForDS(unsigned int pos, unsigned int sizeDS, unsigned size) const;

  /** gets extra interactionBlock corresponding to the present Interaction (see the
   *  top of this files for extra interactionBlock meaning)
   * \param[in,out] InteractionBlock SP::SiconosMatrix
   */
  void getExtraInteractionBlock(SP::SiconosMatrix InteractionBlock) const;

};

#endif // INTERACTION_H

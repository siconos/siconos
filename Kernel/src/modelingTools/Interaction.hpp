
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

/*! \file Interaction.hpp
  \brief Interaction class and related typedef
*/


#ifndef INTERACTION_H
#define INTERACTION_H

// const
#include "BlockVector.hpp"
#include "DynamicalSystemsSet.hpp"
#include "Tools.hpp"
#include "NonSmoothLaw.hpp"
#include "RelationNamespace.hpp"
#include "Relation.hpp"
#include "SiconosPointers.hpp"
#include "LagrangianDS.hpp"
#include "FirstOrderLinearTIDS.hpp"

class DynamicalSystem;
class BlockVector;

typedef std::vector<SP::SiconosMemory> VectorOfMemories;

/**  An Interaction describes the non-smooth interactions between some
 *  Dynamical Systems.
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) Apr 29, 2004
 *
 * An interaction represents the "link" between a set of Dynamical
 * Systems (var: involvedDS) that interact through some relations
 * (between state variables (x,R) and local variables (y,lambda))
 * completed by a non-smooth law.
 *
 * Thus, the interaction main members are:
 *
 * - a set of Dynamical Systems (from 1 to ...) that interacts, named
     involvedDS.
 *
 * - relation: a pointer to a Relation object that determines the type
 *   of relation and so the way it is computed. Warning: there is only
 *   one Relation object (ie only one type of relation for an
 *   interaction) but there can be several "relations", in the sense
 *   of constraints equations between (y,lambda) and (x,r).
 *
 * - nslaw: the non smooth law
 *
 * - the local variable y  (its size is interactionSize).
 *   STL vectors are used and y[i] usually represents the
 *   i-eme derivative of variable y.
 *
 * - the local variable lambda  (its size is interactionSize).
 *   STL vectors are used and lambda[i] represents various level
 *   of multiplier involved in the nonsmooth law and usd by the OSI
 *
 *   y is a container of SiconosVector.
 *
 *
 *
 */
class Interaction : public boost::enable_shared_from_this<Interaction >
{
private:
  /** serialization hooks
  */
  ACCEPT_SERIALIZATION(Interaction);


  /**initialization flag */
  bool _initialized;

  /** name of the Interaction */
  std::string  _id;

  /** number specific to each Interaction */
  int _number;

  /** relative degree of this interaction */
  unsigned int _relativeDegree;

  /** Lowel level for storing output y
   *  y will be initialized from
   *  y[_lowerLevelForOutput] to y[_upperLevelForOutput]
   */
  unsigned int _lowerLevelForOutput;

  /** Upper level for storing output y
    *  y will be initialized from
    *  y[_lowerLevelForOutput] to y[_upperLevelForOutput]
    */
  unsigned int _upperLevelForOutput;

  /** Lowel level for storing input lambda
   *  lambda will be initialized from
   *  lambda[_lowerLevelForIntput] to lambda[_upperLevelForInput]
   */
  unsigned int _lowerLevelForInput;

  /** Upper level for storing input lambda
   *  lambda will be initialized from
   *  lambda[_lowerLevelForIntput] to lambda[_upperLevelForInput]
   */
  unsigned int _upperLevelForInput;

  /** size of the interaction, ie size of y[i] and lambda[i] */
  unsigned int _interactionSize;

  // /** number of relations in the interaction ( equal to
  //     interactionSize / nsLawSize ) */
  // unsigned int _numberOfRelations;

  /** sum of all DS sizes, for DS involved in the interaction */
  unsigned int _sizeOfDS;

  /** sum of all z sizes, for DS involved in the interaction */
  unsigned int _sizeZ;

  /** Absolute position in the "global" vector of constraints (for
      example, the one handled by lsodar) */
  unsigned int _absolutePosition;

  /** Absolute position in the "global" vector of constraints for the proj formulation. */
  unsigned int _absolutePositionProj;


  /** relation between constrained variables and states variables
   * vector of output derivatives
   * y[0] is y, y[1] is yDot and so on
   */
  VectorOfVectors _y;

  /** previous value of Newton iteration for y */
  VectorOfVectors _yOld;

  /** value of the previous time-step */
  VectorOfVectors _y_k;

  /** memory of previous coordinates of the system */
  VectorOfMemories _yMemory;

  /** Size (depth) of the Memory*/
  unsigned int _steps;

  /** result of the computeInput function */
  VectorOfVectors _lambda;

  /** previous step values for lambda */
  VectorOfVectors _lambdaOld;

  /** the Dynamical Systems concerned by this interaction */
  SP::DynamicalSystemsSet _involvedDS;

  /** the Non-smooth Law of the interaction*/
  SP::NonSmoothLaw _nslaw;

  /** the type of Relation of the interaction */
  SP::Relation _relation;

  /** the XML object linked to the Interaction to read XML data */
  SP::InteractionXML _interactionxml;

  /** Work vectors to save pointers to state-related data of the
      dynamical systems involved in the Interaction.*/
  SP::SiconosVector _workX;
  SP::SiconosVector _workXq;
  SP::SiconosVector _workFree;

  SP::SiconosVector _workYp;

  /** Work vector to save pointers to z data of the dynamical systems
      involved in the Interaction.*/
  SP::SiconosVector _workZ;


  // === PRIVATE FUNCTIONS ===



  /** copy constructor => private, no copy nor pass-by-value.
   */
  Interaction(const Interaction& inter);

public:

  /** default constructor */
  Interaction() {};

  /** constructor with XML object of the Interaction
   *  \param InteractionXML* : the XML object corresponding
   *  \param a set of DynamicalSystems
   */
  Interaction(SP::InteractionXML, SP::DynamicalSystemsSet);

  /** constructor with a set of data (only one DS in the Interaction)
      - Note: no id.
   *  \param a SP::DynamicalSystem: the DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int size of the interaction (interactionSize), i.e, the size of the input and output
   *  \param SP::NonSmoothLaw : a pointer to the non smooth law
   *  \param SP::Relation : a pointer to the Relation
   */
  Interaction(SP::DynamicalSystem, int, int, SP::NonSmoothLaw, SP::Relation);
  /** constructor with a set of data (only one DS in the Interaction)
   *  \param string: the id of this Interaction
   *  \param a SP::DynamicalSystem: the DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : size of the interaction (interactionSize), i.e, the size of the input and output
   *  \param SP::NonSmoothLaw : a pointer to the non smooth law
   *  \param SP::Relation : a pointer to the Relation
   */
  Interaction(const std::string&, SP::DynamicalSystem, int, int, SP::NonSmoothLaw, SP::Relation);

  /** constructor with a set of data - Note: no id.
   *  \param a DynamicalSystemsSet: the set of DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : size of the interaction (interactionSize), i.e, the size of the input and output
   *  \param SP::NonSmoothLaw : a pointer to the non smooth law
   *  \param SP::Relation : a pointer to the Relation
   */
  Interaction(DynamicalSystemsSet&, int, int, SP::NonSmoothLaw, SP::Relation);

  /** constructor with a set of data
   *  \param string: the id of this Interaction
   *  \param a DynamicalSystemsSet: the set of DS involved in the Interaction
   *  \param int : the number of this Interaction
   *  \param int : size of the interaction (interactionSize), i.e, the size of the input and output
   *  \param SP::NonSmoothLaw : a pointer to the non smooth law
   *  \param SP::Relation : a pointer to the Relation
   */
  Interaction(const std::string&, DynamicalSystemsSet&, int, int, SP::NonSmoothLaw, SP::Relation);

  /** constructor with no data
   *  \param int : size of the interaction (interactionSize), i.e, the size of the input and output
   *  \param SP::NonSmoothLaw : a pointer to the non smooth law
   *  \param SP::Relation : a pointer to the Relation
   *  \param int : the number of this Interaction (default 0)
   */
  Interaction(int, SP::NonSmoothLaw, SP::Relation, int = 0);

  /** destructor
   */
  ~Interaction();

  /** allocate memory for y[i] and lambda[i] and set them to zero.
   * \param time for initialization.
   * \param the number of required derivatives for y.
   */
  void initialize(double);

  /** check if Interaction is initialized
   */
  bool isInitialized() const
  {
    return _initialized;
  }


  /** build Y and Lambda stl vectors.
  */
  void initializeMemory();

  // === GETTERS/SETTERS ===

  /** get the id of this Interaction
  *  \return the string, id of this Interaction
  */
  inline const std::string  getId() const
  {
    return _id;
  }

  /** set the id of this Interaction
   *  \param the integer to set the id
   */
  inline void setId(const int newId)
  {
    _id = newId;
  }

  /** set the id of this Interaction
   *  \param the string  to set the id
   */
  inline void setId(const std::string& newId)
  {
    _id = newId;
  }

  /** get the value of number
   *  \return the value of number
   */
  inline int number() const
  {
    return _number;
  }

  /** set number
  *  \param int number : the value to set number
  */
  inline void setNumber(const int newNumber)
  {
    _number = newNumber;
  }


  /** get the relative degree
   * \return an unsigned int
   */
  inline unsigned int getRelativeDegree() const
  {
    return _relativeDegree;
  };

  /** set the relative degree
   * \param an unsigned int
   */
  inline void setRelativeDegree(const unsigned int newVal)
  {
    _relativeDegree = newVal;
  };


  /** set the lower level for output y
   * \param an unsigned int
   */
  inline void setLowerLevelForOutput(const unsigned int newVal)
  {
    _lowerLevelForOutput = newVal;
  };

  /** set the upper level for output y
   * \param an unsigned int
   */
  inline void setUpperLevelForOutput(const unsigned int newVal)
  {
    _upperLevelForOutput = newVal;
  };

  /** get the lower level for output y
   */
  inline unsigned int lowerLevelForOutput()
  {
    return _lowerLevelForOutput;
  };

  /** get the upper level for output y
   */
  inline unsigned int  upperLevelForOutput()
  {
    return _upperLevelForOutput;
  };

  /** set the lower level for input Lambda
   * \param an unsigned int
   */
  inline void setLowerLevelForInput(const unsigned int newVal)
  {
    _lowerLevelForInput = newVal;
  };

  /** set the upper level for input Lambda
   * \param an unsigned int
   */
  inline void setUpperLevelForInput(const unsigned int newVal)
  {
    _upperLevelForInput = newVal;
  };


  /** get the lower level for input Lambda
   */
  inline unsigned int lowerLevelForInput()
  {
    return _lowerLevelForInput ;
  };

  /** get the upper level for input Lambda
   */
  inline unsigned int upperLevelForInput()
  {
    return _upperLevelForInput;
  };


  /** get the dimension of the interaction (y and lambda size)
  *  \return an unsigned int
  */
  inline unsigned int getSizeOfY() const
  {
    return _interactionSize;
  }

  /** set the dimension of the Interaction
  *  \param an unsigned int
  */
  inline void setInteractionSize(const unsigned int newVal)
  {
    _interactionSize = newVal;
  }

  //  /** get the number of relations in the interaction
  //  *  \return an unsigned int
  //  */
  // inline unsigned int numberOfRelations() const {return _numberOfRelations;}

  /** get the sum of DS sizes, for DS involved in interaction
   *  \return an unsigned int
   */
  inline unsigned int getSizeOfDS() const
  {
    return _sizeOfDS;
  }

  /** get the sum of z sizes, for DS involved in interaction
   *  \return an unsigned int
   */
  inline unsigned int getSizez() const
  {
    return _sizeZ;
  }

  unsigned int absolutePosition()
  {
    return _absolutePosition;
  };
  void setAbsolutePosition(unsigned int v)
  {
    _absolutePosition = v;
  };
  unsigned int absolutePositionProj()
  {
    return _absolutePositionProj;
  };
  void setAbsolutePositionProj(unsigned int v)
  {
    _absolutePositionProj = v;
  };




  // -- y --

  /** get y[i], derivative number i of output
   *  \return BlockVector
   */
  inline const SiconosVector getCopyOfy(const unsigned int i) const
  {
    return *(_y[i]);
  }

  /** get y[i], derivative number i of output
   *  \return BlockVector
   */
  inline const SiconosVector getCopyOfyOld(const unsigned int i) const
  {
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
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector y(const unsigned int i) const
  {
    return _y[i];
  }

  /** set the output vector y to newVector with copy of the y[i] (ie
      memory allocation)
  *  \param VectorOfVectors
  */
  void setY(const VectorOfVectors&);

  /** set the output vector y to newVector with direct pointer
  *  equality for the y[i]
  * \param VectorOfVectors
  */
  void setYPtr(const VectorOfVectors&);

  /** set y[i] to newValue
  *  \param a SiconosVector and an unsigned int
  */
  void setY(const unsigned int , const SiconosVector&);

  /** set y[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setYPtr(const unsigned int , SP::SiconosVector newPtr);

  // -- yOld --

  /** get vector of output derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getYOld() const
  {
    return _yOld;
  }

  /** get yOld[i], derivative number i of output
  *  \return BlockVector
  */
  inline const SiconosVector getYOld(const unsigned int i) const
  {
    return *(_yOld[i]);
  }

  /** get yOld[i], derivative number i of output
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector yOld(const unsigned int i) const
  {
    return _yOld[i];
  }
  /*get y_k[i]
    \return pointer on a SiconosVector
   */
  inline SP::SiconosVector y_k(const unsigned int i) const
  {
    return _y_k[i];
  }

  /** set the output vector yOld to newVector
  *  \param VectorOfVectors
  */
  void setYOld(const VectorOfVectors&);

  /** set vector yOld to newVector with direct pointer equality for
  *  the yOld[i]
  * \param VectorOfVectors
  */
  void setYOldPtr(const VectorOfVectors&);

  /** set yOld[i] to newValue
  *  \param a SiconosVector and an unsigned int
  */
  void setYOld(const unsigned int , const SiconosVector&);

  /** set yOld[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setYOldPtr(const unsigned int , SP::SiconosVector newPtr);


  /** set yOld[i] to pointer newPtr
   *  \param a SP::SiconosVector  and an unsigned int
   */
  void setSteps(unsigned int newval)
  {
    _steps = newval;
  };


  /** get all the values of the state vector y stored in memory
   *  \return a memory
   */
  inline SP::SiconosMemory yMemory(unsigned int level) const
  {
    return _yMemory[level];
  }

  // -- lambda --

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getLambda() const
  {
    return _lambda;
  }

  /** get lambda[i], derivative number i of input
  *  \return SiconosVector
  */
  inline const SiconosVector getLambda(const unsigned int i) const
  {
    assert(_lambda[i]);
    return *(_lambda[i]);
  }

  /** get lambda[i], derivative number i of input
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector lambda(const unsigned int i) const
  {
    assert(_lambda[i]);
    return _lambda[i];
  }

  /** set the input vector lambda to newVector
  *  \param VectorOfVectors
  */
  void setLambda(const VectorOfVectors&);

  /** set vector lambda to newVector with direct pointer equality for the lambda[i]
  *  \param VectorOfVectors
  */
  void setLambdaPtr(const VectorOfVectors&);

  /** set lambda[i] to newValue
  *  \param a SiconosVector and an unsigned int
  */
  void setLambda(const unsigned int , const SiconosVector&);

  /** set lambda[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setLambdaPtr(const unsigned int , SP::SiconosVector newPtr);

  // -- lambdaOld --

  /** get vector of input derivatives
  *  \return a VectorOfVectors
  */
  inline const VectorOfVectors getLambdaOld() const
  {
    return _lambdaOld;
  }

  /** get lambdaOld[i], derivative number i of input
  *  \return SiconosVector
  */
  inline const SiconosVector getLambdaOld(const unsigned int i) const
  {
    return *(_lambdaOld[i]);
  }

  /** get lambdaOld[i], derivative number i of input
  *  \return pointer on a SiconosVector
  */
  inline SP::SiconosVector lambdaOld(const unsigned int i) const
  {
    return _lambdaOld[i];
  }

  /** set the input vector lambdaOld to newVector
  *  \param VectorOfVectors
  */
  void setLambdaOld(const VectorOfVectors&);

  /** set vector lambdaOld to newVector with direct pointer equality for the lambdaOld[i]
  *  \param VectorOfVectors
  */
  void setLambdaOldPtr(const VectorOfVectors&);

  /** set lambdaOld[i] to newValue
  *  \param a SiconosVector and an unsigned int
  */
  void setLambdaOld(const unsigned int , const SiconosVector&);

  /** set lambdaOld[i] to pointer newPtr
  *  \param a SP::SiconosVector  and an unsigned int
  */
  void setLambdaOldPtr(const unsigned int , SP::SiconosVector newPtr);

  /** insert a Dynamical system
   *  \param a SP::DynamicalSystem
   */
  void insert(SP::DynamicalSystem ds)
  {
    _involvedDS->insert(ds);
  };

  /** gets an iterator to the first element of the involvedDS set.
   *  \return a DSIterator.
   */
  inline DSIterator dynamicalSystemsBegin()
  {
    return _involvedDS->begin();
  };

  /** gets an iterator equal to involvedDS->end().
   *  \return a DSIterator.
   */
  inline DSIterator dynamicalSystemsEnd()
  {
    return _involvedDS->end();
  };

  /** gets a const iterator to the first element of the involvedDS set.
   *  \return a ConstDSIterator.
   */
  inline ConstDSIterator dynamicalSystemsBegin() const
  {
    return _involvedDS->begin();
  };

  /** gets a const iterator equal to _involvedDS->end().
   *  \return a ConstDSIterator.
   */
  inline ConstDSIterator dynamicalSystemsEnd() const
  {
    return _involvedDS->end();
  };

  /** get a pointer to the DynamicalSystems of this Interaction
   *  \return a DynamicalSystemsSet*
   */
  inline SP::DynamicalSystemsSet dynamicalSystems()
  {
    return _involvedDS;
  }

  /** set the _involvedDS
  *  \param a DynamicalSystemsSet
  */
  void setDynamicalSystems(const DynamicalSystemsSet&) ;

  /** get a specific DynamicalSystem
  *  \param the identification number of the wanted DynamicalSystem
  *  \return a pointer on Dynamical System
  */
  SP::DynamicalSystem dynamicalSystem(int);

  /** get the Relation of this Interaction
   *  \return a pointer on this Relation
   */
  inline SP::Relation relation() const
  {
    return _relation;
  }

  /** set the Relation of this Interaction
  *  \param the SP::relation to set
  */
  void setRelationPtr(SP::Relation newRelation) ;

  /** get the NonSmoothLaw of this Interaction
  *  \return a pointer on this NonSmoothLaw
  */
  inline SP::NonSmoothLaw nonSmoothLaw() const
  {
    return _nslaw;
  }
  inline SP::NonSmoothLaw nslaw() const
  {
    return _nslaw;
  }

  /** set the NonSmoothLaw of this Interaction
  *  \param the SP::NonSmoothLaw to set
  */
  void setNonSmoothLawPtr(SP::NonSmoothLaw newNslaw) ;

  /** get the size of the NonSmoothLaw of this Interaction
   * \return the size of the NonSmoothLaw
   */
  inline unsigned int getNonSmoothLawSize() const
  {
    return _nslaw->size();
  }

  /** function used to sort Interaction in SiconosSet<SP::Interaction>
   *  \return a double* (warning: must be const, despite intel compilers warning, because of SiconosSet Cmp function arguments)
   */
  inline double* getSort() const
  {
    return (double*)this;
  }

  // --- OTHER FUNCTIONS ---

  /** compute sum of all interaction-involved DS sizes
  */
  void computeSizeOfDS();

  /**   put values of y into yOld, the same for lambda
  */
  void swapInMemory();

  /** Must be call to fill _y_k. (after convergence of the Newton iterations)
   */
  void swapTimeStepInMemory();
  /** print the data to the screen
  */
  void display() const;

  /** Computes output y; depends on the relation type.
   *  \param double : current time
   *  \param unsigned int: number of the derivative to compute,
   *  optional, default = 0.
   */
  void computeOutput(double, unsigned int = 0);

  /** Compute input r of all Dynamical Systems involved in the present
   *   Interaction.
   *  \param double : current time
   *  \param unsigned int: order of lambda used to compute input.
   */
  void computeInput(double, unsigned int);


  /** to set _workFree content.
   */
  inline void setWorkFree()
  {
    assert(_workFree) ;

    if (_workFree->size() == 0)
    {
      _workFree->resize(_sizeOfDS);
    }

    ConstDSIterator itDS;
    int index;
    for (index = 0, itDS = dynamicalSystemsBegin();
         itDS != dynamicalSystemsEnd();
         ++itDS)
    {
      _workFree->setBlock(index, *((*itDS)->workFree()));
      index += (*itDS)->getDim();
    }
  };

  inline void setWorkX()
  {
    assert(_workX) ;

    if (_workX->size() == 0)
    {
      _workX->resize(_sizeOfDS);
    }

    ConstDSIterator itDS;
    int index;
    for (index = 0, itDS = dynamicalSystemsBegin();
         itDS != dynamicalSystemsEnd();
         ++itDS)
    {
      _workX->setBlock(index, *((*itDS)->x()));
      index += (*itDS)->getDim();
    }
  };

  inline void setWorkXFromVelocity()
  {
    assert(_workX) ;

    if (_workX->size() == 0)
    {
      _workX->resize(_sizeOfDS);
    }

    ConstDSIterator itDS;
    int index;
    for (index = 0, itDS = dynamicalSystemsBegin();
         itDS != dynamicalSystemsEnd();
         ++itDS)
    {
      assert(Type::value(**itDS) == Type::LagrangianDS ||
             Type::value(**itDS) == Type::LagrangianLinearTIDS);
      _workX->setBlock(index, *boost::static_pointer_cast<LagrangianDS>(*itDS)->velocity());
      index += (*itDS)->getDim();
    }
  };

  inline void setWorkXq()
  {
    assert(_workXq) ;

    if (_workXq->size() == 0)
    {
      _workXq->resize(_sizeOfDS);
    }

    ConstDSIterator itDS;
    int index;
    for (index = 0, itDS = dynamicalSystemsBegin();
         itDS != dynamicalSystemsEnd();
         ++itDS)
    {
      assert(Type::value(**itDS) == Type::FirstOrderNonLinearDS ||
             Type::value(**itDS) == Type::FirstOrderLinearDS ||
             Type::value(**itDS) == Type::FirstOrderLinearTIDS);
      _workXq->setBlock(index, *boost::static_pointer_cast<FirstOrderNonLinearDS>(*itDS)->xq());
      index += (*itDS)->getDim();
    }
  };

  inline void setWorkZ()
  {
    assert(_workZ) ;

    if (_workZ->size() == 0)
    {
      _workZ->resize(_sizeZ);
    }

    ConstDSIterator itDS;
    int index;
    for (index = 0, itDS = dynamicalSystemsBegin();
         itDS != dynamicalSystemsEnd();
         ++itDS)
    {
      SiconosVector& tmpz = *((*itDS)->z());
      _workZ->setBlock(index, tmpz);
      index += tmpz.size();
    }
  };

  /** Get a pointer to workX */
  inline SP::SiconosVector workX()
  {
    return _workX;
  };

  inline SP::SiconosVector workXq()
  {
    return _workXq;
  };

  inline SP::SiconosVector workFree()
  {
    return _workFree;
  };

  inline SP::SiconosVector yp()
  {
    return _workYp;
  };


  /** Get a pointer to workZ */
  inline SP::SiconosVector workZ()
  {
    return _workZ;
  };

  /** gets the matrix used in interactionBlock computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting interactionBlock matrix
   */
  void getLeftInteractionBlockForDS(SP::DynamicalSystem, SP::SiconosMatrix) const;

  /** gets the matrix used in interactionBlock computation. Used only for the formulation projecting on the constraints.
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting interactionBlock matrix
   */
  void getLeftInteractionBlockForDSProjectOnConstraints(SP::DynamicalSystem ds, SP::SiconosMatrix InteractionBlock) const;
  /** gets the matrix used in interactionBlock computation, (left * W * rigth), depends on the relation type (ex, LinearTIR, left = C, right = B).
   *         We get only the part corresponding to ds.
   *  \param a pointer to a dynamical system
   *  \param a pointer to SiconosMatrix (in-out parameter): the resulting interactionBlock matrix
   */
  void getRightInteractionBlockForDS(SP::DynamicalSystem, SP::SiconosMatrix) const;

  /** gets extra interactionBlock corresponding to the present Interaction (see the
   *  top of this files for extra interactionBlock meaning)
   * \param[in,out] InteractionBlock SP::SiconosMatrix
   */
  void getExtraInteractionBlock(SP::SiconosMatrix InteractionBlock) const;

  inline double getYRef(unsigned int i) const
  {
    // get the single value used to build indexSets Warning: the
    // relativePosition depends on NsLawSize and/or type.  This means
    // that at the time, for the interactionBlock of y that corresponds to
    // the present relation, the first scalar value is used.  For
    // example, for friction, normal part is in first position, followed
    // by the tangential parts.
    return (*y(i))(0);
  }

  inline double getLambdaRef(unsigned int i) const
  {
    // get the single value used to build indexSets
    return (*lambda(i))(0);
  }

  // --- XML RELATED FUNCTIONS ---

  /** get the InteractionXML* of the Interaction
   *  \return InteractionXML* : the pointer on the InteractionXML
   */
  inline SP::InteractionXML interactionXML() const
  {
    return _interactionxml;
  }

  /** set the InteractionXML* of the Interaction
  *  \param InteractionXML* :  the pointer to set
  */
  inline void setInteractionXMLPtr(SP::InteractionXML interxml)
  {
    _interactionxml = interxml;
  }

  /** copy the data of the Interaction to the XML tree
  *  \exception RuntimeException
  */
  void saveInteractionToXML();

  /** returns the type of the embedded relation.
   */
  inline RELATION::TYPES getRelationType() const
  {
    return _relation->getType();
  };

  /** returns the subtype of the embedded relation.
   */
  inline RELATION::SUBTYPES getRelationSubType() const
  {
    return _relation->getSubType();
  } ;

};

#endif // INTERACTION_H

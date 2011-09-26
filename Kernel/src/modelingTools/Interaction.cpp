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
#include <assert.h>
#include "Interaction.hpp"
#include "InteractionXML.hpp"
#include "NonSmoothLawXML.hpp"
#include "RelationXML.hpp"
#include "RelationTypes.hpp"
#include "ComplementarityConditionNSL.hpp"
#include "RelayNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "DynamicalSystem.hpp"

#include "NewtonEulerR.hpp" // ??
#include "NewtonEulerDS.hpp" // ??

using namespace std;
using namespace RELATION;

// --- CONSTRUCTORS ---

// --- XML constructor ---
Interaction::Interaction(SP::InteractionXML interxml, SP::DynamicalSystemsSet nsdsSet):
  _initialized(false), _id("undefined"), _number(0), _interactionSize(0), _numberOfRelations(0),
  _sizeOfDS(0), _sizeZ(0), _interactionxml(interxml)
{
  assert(_interactionxml && "NULL pointer");

  // id and number
  if (_interactionxml->hasId()) _id = _interactionxml->getId();
  _number = _interactionxml->number();

  // interaction size
  _interactionSize = _interactionxml->getSize();

  // --- Non smooth law ---
  string NslawType = _interactionxml->getNonSmoothLawXML()->getType();
  // ComplementarityConditionNSL
  if (NslawType == COMPLEMENTARITY_CONDITION_NSLAW_TAG)
  {
    _nslaw.reset(new ComplementarityConditionNSL(_interactionxml->getNonSmoothLawXML()));
  }

  // RelayNSL
  else if (NslawType == RELAY_NSLAW_TAG)
  {
    _nslaw.reset(new RelayNSL(_interactionxml->getNonSmoothLawXML()));
  }

  // NewtonImpactNSL
  else if (NslawType == NEWTON_IMPACT_NSLAW_TAG)
  {
    _nslaw.reset(new NewtonImpactNSL(_interactionxml->getNonSmoothLawXML()));
  }
  // Newton impact friction law
  else if (NslawType == NEWTON_IMPACT_FRICTION_NSLAW_TAG)
  {
    _nslaw.reset(new NewtonImpactFrictionNSL(_interactionxml->getNonSmoothLawXML()));
  }
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown NSLAW type");

  // --- Dynamical Systems ---
  _involvedDS.reset(new DynamicalSystemsSet());
  if (!nsdsSet->isEmpty())
  {
    // Get a list of DS concerned from xml

    if (_interactionxml->hasAllDS())
      _involvedDS->insert(nsdsSet->begin(), nsdsSet->end());

    else
    {
      // get numbers of DS involved in the interaction from the xml input file.
      std::vector<int> dsNumbers;
      _interactionxml->getDSNumbers(dsNumbers);
      // get corresponding DS and insert them into the involvedDS set.
      for (vector<int>::iterator it = dsNumbers.begin(); it != dsNumbers.end(); ++it)
        _involvedDS->insert(nsdsSet->getPtr(*it));
    }
  }
  else cout << "Interaction constructor, warning: no dynamical systems linked to the interaction!" << endl;

  // --- Relation ---
  RELATION::TYPES relationType = _interactionxml->getRelationXML()->getType();

  // First Order Non Linear Relation
  if (relationType == FirstOrder)
  {
    RELATION::SUBTYPES relationSubType = _interactionxml->getRelationXML()->getSubType();
    if (relationSubType == Type1R)
      _relation.reset(new FirstOrderType1R(_interactionxml->getRelationXML()));
    // Linear relation
    else if (relationSubType == LinearR)
      _relation.reset(new FirstOrderLinearR(_interactionxml->getRelationXML()));
    // Linear time-invariant coef. relation
    else if (relationSubType == LinearTIR)
      _relation.reset(new FirstOrderLinearTIR(_interactionxml->getRelationXML()));
  }
  // Lagrangian non-linear relation
  else if (relationType == Lagrangian)
  {
    RELATION::SUBTYPES relationSubType = _interactionxml->getRelationXML()->getSubType();
    // \todo create a factory to avoid "if" list for Relation
    // construction according to subType.
    if (relationSubType == ScleronomousR)
      _relation.reset(new LagrangianScleronomousR(_interactionxml->getRelationXML()));
    else if (relationSubType == RheonomousR)
      _relation.reset(new LagrangianRheonomousR(_interactionxml->getRelationXML()));
    else if (relationSubType == CompliantR)
      _relation.reset(new LagrangianCompliantR(_interactionxml->getRelationXML()));
    // Lagrangian linear relation
    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      _relation.reset(new LagrangianLinearTIR(_interactionxml->getRelationXML()));
  }
  else RuntimeException::selfThrow("Interaction::xml constructor, unknown relation type " + _relation->getType());

  // check coherence between interactionSize and nsLawSize
  if ((_interactionSize % _nslaw->size()) != 0)
    RuntimeException::selfThrow("Interaction::xml constructor, inconsistency between interaction size and non smooth law size.");

  // download xml values for y[0] and lambda[0] forbidden since y and
  // lambda are not allocated yet (this is done in initialize
  // function, called by simulation !!)
  if (_interactionxml->hasY() ||  _interactionxml->hasLambda())
    RuntimeException::selfThrow("Interaction::xml constructor, y or lambda download is forbidden.");
}

// --- Constructors from a set of data ---

Interaction::Interaction(SP::DynamicalSystem ds, int newNumber, int nInter,
                         SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false) , _id("none"), _number(newNumber), _interactionSize(nInter), _numberOfRelations(1),
  _sizeOfDS(0), _sizeZ(0), _y(2), _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  _involvedDS->insert(ds); // Warning: insert pointer to DS!!

}
Interaction::Interaction(const string& newId, SP::DynamicalSystem ds,
                         int newNumber, int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false), _id(newId), _number(newNumber), _interactionSize(nInter), _numberOfRelations(1),
  _sizeOfDS(0), _sizeZ(0), _y(2), _nslaw(newNSL),  _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  _involvedDS->insert(ds); // Warning: insert pointer to DS!!
}


Interaction::Interaction(DynamicalSystemsSet& dsConcerned, int newNumber, int nInter,
                         SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false) , _id("none"), _number(newNumber), _interactionSize(nInter), _numberOfRelations(1),
  _sizeOfDS(0), _sizeZ(0), _y(2), _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  DSIterator itDS;
  for (itDS = dsConcerned.begin(); itDS != dsConcerned.end(); ++itDS)
    _involvedDS->insert(*itDS); // Warning: insert pointers to DS!!
}

Interaction::Interaction(const string& newId, DynamicalSystemsSet& dsConcerned, int newNumber,
                         int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false) , _id(newId), _number(newNumber), _interactionSize(nInter), _numberOfRelations(1), _sizeOfDS(0), _sizeZ(0),
  _y(2),  _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  DSIterator itDS;
  for (itDS = dsConcerned.begin(); itDS != dsConcerned.end(); ++itDS)
    _involvedDS->insert(*itDS); // Warning: insert pointers to DS!!
}

/* initialisation with empty set */
Interaction::Interaction(int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel, int newNumber):
  _initialized(false), _number(newNumber), _interactionSize(nInter), _numberOfRelations(1), _sizeOfDS(0), _sizeZ(0),
  _y(2),  _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
}


// --- DESTRUCTOR ---
Interaction::~Interaction()
{
}

void Interaction::initialize(double t0)
{

  if (!_initialized)
  {
    assert(relation() && "Interaction::initialize failed, relation() == NULL");

    assert(nslaw() && "Interaction::initialize failed, non smooth law == NULL");

    computeSizeOfDS();

    DSIterator itDS;
    for (itDS = dynamicalSystemsBegin(); itDS != dynamicalSystemsEnd(); ++itDS)
    {
      assert(_lowerLevelForInput <= _upperLevelForInput);
      for (unsigned int k = _lowerLevelForInput ; k < _upperLevelForInput + 1; k++)
      {
        (*itDS)->initializeNonSmoothInput(k);
      }
    }

    relation()->setInteractionPtr(shared_from_this());

    // compute number of relations.
    _numberOfRelations = _interactionSize / nslaw()->size();

    if (_numberOfRelations > 1)
    {
      RuntimeException::selfThrow("Interaction::initialize() - _numberOfRelations > 1. Obsolete !");
    }

    initializeMemory();
    relation()->initialize(shared_from_this());


    if (_steps > 1) // Multi--step methods
    {
      // Comoyte the old Values of Output with stored values in Memory
      for (unsigned int k = 0; k < _steps - 1; k++)
      {
        /** ComputeOutput to fill the Memory
         * We assume the state x is stored in xMemory except for the  initial
         * condition which has not been swap yet.
         */
        relation()->LinkDataFromMemory(k);
        for (unsigned int i = 0; i < _upperLevelForOutput + 1; ++i)
        {
          computeOutput(t0, i);
          // std::cout <<"_y["<<i<<"] = (link Memory level ["<< k <<"]) " << std::endl;
          // _y[i]->display();
          _yMemory[i]->swap(_y[i]);
        }

      }
      relation()->LinkData();
    }

    // Compute y values for t0
    for (unsigned int i = 0; i < _upperLevelForOutput + 1; ++i)
    {
      computeOutput(t0, i);
      // std::cout <<"_y["<<i<<"] = (current value) " << std::endl;
      // _y[i]->display();
    }


    _workX.reset(new BlockVector());
    _workZ.reset(new BlockVector());
    _workXq.reset(new BlockVector());
    _workFree.reset(new BlockVector());
    _workYp.reset(new SimpleVector(nslaw()->size()));

    _initialized = true;
  }

}

// Initialize and InitializeMemory are separated in two functions
// since we need to know the relative degree to know
// "numberOfDerivatives", while numberOfRelations and the size of the
// non smooth law are required inputs to compute the relative degree.
void Interaction::initializeMemory()
{
  // Warning: this function is called from Simulation initialize,
  // since we need to know :
  // the levels _lowerLevelForOutput and _upperLevelForOutput to size Y
  // and the levels _lowerLevelForInput and _upperLevelForInput to size  Lambda.
  // this depends on many criteria (simulation type, osi type, ds type, nonsmooth type)
  // and they are computed in Simulation::ComputeLevelsForInputAndOutput

  // Memory allocation for y and lambda

  assert(_upperLevelForOutput >= 0);
  assert(_upperLevelForOutput >= _lowerLevelForOutput);
  assert(_upperLevelForInput >= 0);
  assert(_upperLevelForInput >= _lowerLevelForInput);

  // in order to simplify we size from 0 to _upperLevelForXXX
  _y.resize(_upperLevelForOutput + 1) ;
  _yOld.resize(_upperLevelForOutput + 1);
  _y_k.resize(_upperLevelForOutput + 1);
  _lambda.resize(_upperLevelForInput + 1);
  _lambdaOld.resize(_upperLevelForInput + 1);

  _yMemory.resize(_upperLevelForOutput + 1);



  // get the dimension of the non smooth law, ie the size of a unitary blocks (one per relation)
  unsigned int nslawSize = nslaw()->size();
  relation()->initializeMemory();

  for (unsigned int i = _lowerLevelForOutput ;
       i < _upperLevelForOutput + 1 ;
       i++)
  {
    _y[i].reset(new SimpleVector(nslawSize));
    _yOld[i].reset(new SimpleVector(nslawSize));
    _y_k[i].reset(new SimpleVector(nslawSize));
    assert(_steps > 0);
    _yMemory[i].reset(new SiconosMemory(_steps));
    _y[i]->zero();
    _yOld[i]->zero();
    _y_k[i]->zero();
  }
  for (unsigned int i = _lowerLevelForInput ;
       i < _upperLevelForInput + 1 ;
       i++)
  {
    _lambda[i].reset(new SimpleVector(nslawSize));
    _lambdaOld[i].reset(new SimpleVector(nslawSize));
    _lambdaOld[i]->zero();
  }

}

// --- GETTERS/SETTERS ---

void Interaction::setY(const VectorOfVectors& newVector)
{

  unsigned int size = newVector.size();

  _y.clear();
  _y.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _y[i].reset(new SimpleVector(*(newVector[i]))); // -> copy !
}

void Interaction::setYPtr(const VectorOfVectors& newVector)
{
  _y.clear();

  // copy
  _y = newVector; // smart ptr
}

void Interaction::setY(const unsigned int  index, const SiconosVector& newY)
{
  assert(_y.size() > index &&
         "Interaction::setY, index out of range ");

  // set y[index]
  if (! _y[index])
  {
    _y[index].reset(new SimpleVector(newY));
  }
  else
  {
    assert(_y[index]->size() == newY.size() &&
           "Interaction::setY(index,newY), inconsistent sizes between y(index) and newY ");
    *(_y[index]) = newY;
  }
}

void Interaction::setYPtr(const unsigned int  index, SP::SiconosVector newY)
{
  assert(_y.size() > index &&
         "Interaction::setYPtr, index out of range");

  assert(newY->size() == _interactionSize &&
         "Interaction::setYPtr, interactionSize differs from newY vector size");

  assert(!newY->isBlock() &&
         "Interaction::setYPtr(newY), newY is a block vector!");

  _y[index] = boost::static_pointer_cast<SimpleVector>(newY);
}

void Interaction::setYOld(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();
  _yOld.clear();
  _yOld.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _yOld[i].reset(new SimpleVector(*(newVector[i]))); // -> copy !
}

void Interaction::setYOldPtr(const VectorOfVectors& newVector)
{
  // clear _yOld

  _yOld.clear();

  // copy
  _yOld = newVector; // smart ptr
}

void Interaction::setYOld(const unsigned int  index, const SiconosVector& newYOld)
{
  if (_yOld.size() <= index)
    RuntimeException::selfThrow("Interaction::setYOld, index out of range ");

  // set _yOld[index]
  if (! _yOld[index])
  {
    _yOld[index].reset(new SimpleVector(newYOld));
  }
  else
  {
    assert(_yOld[index]->size() == newYOld.size() &&
           "Interaction::setYOld(index,newYOld), inconsistent sizes between yOld(index) and newYOld");
    *(_yOld[index]) = newYOld;
  }
}

void Interaction::setYOldPtr(const unsigned int  index, SP::SiconosVector newYOld)
{
  assert(_yOld.size() > index &&
         "Interaction::setYOldPtr, index out of range");

  assert(newYOld->size() == _interactionSize &&
         "Interaction::setYOldPtr, interactionSize differs from newYOld vector size");

  assert((! newYOld->isBlock()) &&
         "Interaction::setYOldPtr(newY), newY is a block vector!");

  // set _yOld[index]
  _yOld[index] = boost::static_pointer_cast<SimpleVector>(newYOld);
}

void Interaction::setLambda(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();
  _lambda.clear();
  _lambda.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _lambda[i].reset(new SimpleVector(*(newVector[i]))); // -> copy !
}

void Interaction::setLambdaPtr(const VectorOfVectors& newVector)
{
  _lambda.clear();

  _lambda = newVector; // smart ptr
}

void Interaction::setLambda(const unsigned int  index, const SiconosVector& newLambda)
{
  assert(_lambda.size() <= index &&
         "Interaction::setLambda, index out of range");

  // set lambda[index]
  if (! _lambda[index])
  {
    _lambda[index].reset(new SimpleVector(newLambda));
  }
  else
  {
    assert(_lambda[index]->size() == newLambda.size() &&
           "Interaction::setLambda(index,newLambda), inconsistent sizes between lambda(index) and newLambda");
    *(_lambda[index]) = newLambda;
  }
}

void Interaction::setLambdaPtr(const unsigned int  index, SP::SiconosVector newLambda)
{
  assert(_lambda.size() > index &&
         "Interaction::setLambdaPtr, index out of range ");

  assert(newLambda->size() == _interactionSize &&
         "Interaction::setLambdaPtr, interactionSize differs from newLambda vector size ");

  assert(!newLambda->isBlock() &&
         "Interaction::setLambdaPtr(newLambda), newLambda is  a block vector! ");

  // set lambda[index]
  _lambda[index] = boost::static_pointer_cast<SimpleVector>(newLambda);
}

void Interaction::setLambdaOld(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();

  // clear lambdaOld
  _lambdaOld.clear();
  _lambdaOld.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _lambdaOld[i].reset(new SimpleVector(*(newVector[i]))); // -> copy !
}

void Interaction::setLambdaOldPtr(const VectorOfVectors& newVector)
{
  // clear lambdaOld
  _lambdaOld.clear();

  // copy
  _lambdaOld = newVector; // smart ptrs
}

void Interaction::setLambdaOld(const unsigned int  index, const SiconosVector& newLambdaOld)
{
  assert(_lambdaOld.size() > index &&
         "Interaction::setLambdaOld, index out of range ");

  // set lambdaOld[index]
  if (! _lambdaOld[index])
  {
    _lambdaOld[index].reset(new SimpleVector(newLambdaOld));
  }
  else
  {
    if (_lambdaOld[index]->size() != newLambdaOld.size())
      RuntimeException::selfThrow("Interaction::setLambdaOld(index,newLambdaOld), inconsistent sizes between lambdaOld(index) and newLambdaOld ");
    *(_lambdaOld[index]) = newLambdaOld;
  }
}

void Interaction::setLambdaOldPtr(const unsigned int  index, SP::SiconosVector newLambdaOld)
{
  if (_lambdaOld.size() > index)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, index out of range ");
  if (newLambdaOld->size() != _interactionSize)
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr, interactionSize differs from newLambdaOld vector size ");
  if (!newLambdaOld->isBlock())
    RuntimeException::selfThrow("Interaction::setLambdaOldPtr(newLambda), newLambda is  a block vector! ");

  // set lambdaOld[index]
  _lambdaOld[index] = boost::static_pointer_cast<SimpleVector>(newLambdaOld);
}


void Interaction::setDynamicalSystems(const DynamicalSystemsSet& newSet)
{
  ConstDSIterator itDS;
  for (itDS = newSet.begin(); itDS != newSet.end(); ++itDS)
    _involvedDS->insert(*itDS); // smart ptrs

  computeSizeOfDS();
}

SP::DynamicalSystem Interaction::dynamicalSystem(int nb)
{
  assert(_involvedDS->isIn(nb) &&  // if ds number nb is not in the set ...
         "Interaction::dynamicalSystem(nb), DS number nb is not in the set.");
  return _involvedDS->getPtr(nb);
}

void Interaction::setRelationPtr(SP::Relation newRelation)
{
  _relation = newRelation;
}

void Interaction::setNonSmoothLawPtr(SP::NonSmoothLaw newNslaw)
{
  _nslaw = newNslaw;
}

// --- OTHER FUNCTIONS ---

void Interaction::computeSizeOfDS()
{
  _sizeOfDS = 0;
  _sizeZ = 0;
  DSIterator it;
  SP::SiconosVector ZP;
  for (it = _involvedDS->begin(); it != _involvedDS->end(); it++)
  {
    _sizeOfDS += (*it)->getDim();
    ZP = (*it)->z();
    if (ZP) _sizeZ += ZP->size();
  }
}

void Interaction::swapInMemory()
{
  // i corresponds to the derivative number and j the relation number.
  for (unsigned int i = _lowerLevelForOutput; i < _upperLevelForOutput + 1 ; i++)
  {
    for (unsigned int j = 0; j < _numberOfRelations; ++j)
    {
      assert(_y[i]->vector(j));
      assert(_yOld[i]->vector(j));

      *(_yOld[i]->vector(j)) = *(_y[i]->vector(j)) ;
    }


  }

  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1  ; i++)
  {
    for (unsigned int j = 0; j < _numberOfRelations; ++j)
    {

      assert(_lambdaOld[i]->vector(j));
      assert(_lambda[i]->vector(j));

      *(_lambdaOld[i]->vector(j)) = *(_lambda[i]->vector(j));
    }
  }



}

void Interaction::swapTimeStepInMemory()
{
  // i corresponds to the derivative number and j the relation number.
  for (unsigned int i = 0; i < _y.size() ; i++)
  {
    for (unsigned int j = 0; j < _numberOfRelations; ++j)
    {
      *(_y_k[i]->vector(j)) = *(_y[i]->vector(j)) ;
    }
    _yMemory[i]->swap(_y[i]);
  }


}

void Interaction::display() const
{
  cout << "======= Interaction display =======" << endl;

  if (_initialized)
    cout << "The interaction is initialized" << endl;
  else
    cout << "The interaction is not initialized" << endl;
  cout << "| id : " << _id << endl;
  cout << "| number : " << _number << endl;
  cout << "| relativeDegree : " << _relativeDegree << endl;
  cout << "| lowerLevelForOutput : " << _lowerLevelForOutput << endl;
  cout << "| upperLevelForOutput : " << _upperLevelForOutput << endl;
  cout << "| lowerLevelForInput : " << _lowerLevelForInput << endl;
  cout << "| upperLevelForInput : " << _upperLevelForInput << endl;
  cout << "| interactionSize : " << _interactionSize << endl;
  cout << "| numberOfRelations : " << _numberOfRelations << endl;
  cout << "|  _sizeOfDS : " << _sizeOfDS << endl;
  cout << "|  _sizeZ: " << _sizeZ << endl;

  cout << "| involved DS :" << endl;
  cout << _involvedDS << endl;
  _relation->display();
  if (_initialized)
  {
    for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++)
    {
      cout << "| y[" << i  << "] : " << endl;
      if (_y[i]) _y[i]->display();
      else cout << "->NULL" << endl;
    }
    for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++)
    {
      cout << "| yOld[" << i  << "] : " << endl;
      if (_yOld[i]) _yOld[i]->display();
      else cout << "->NULL" << endl;
    }
    for (unsigned int i = 0; i < _upperLevelForInput + 1; i++)
    {
      cout << "| lambda[" << i  << "] : " << endl;
      if (_lambda[i]) _lambda[i]->display();
      else cout << "->NULL" << endl;
    }

  }
  cout << "===================================" << endl;
}

void Interaction::computeOutput(double time, unsigned int level)
{
  relation()->computeOutput(time, level);
}

void Interaction::computeInput(double time, unsigned int level)
{
  relation()->computeInput(time, level);
}




// --- XML RELATED FUNCTIONS ---

void Interaction::saveInteractionToXML()
{
  /*
   * save the data of the Interaction
   */

  if (_interactionxml)
  {
    //  _interactionxml->setDSConcerned( involvedDS );
    _interactionxml->setId(_id);
    _interactionxml->setNumber(_number);
    _interactionxml->setSize(_interactionSize);
    _interactionxml->setY(*(_y[0]));
    _interactionxml->setLambda(*(_lambda[0]));
  }
  else RuntimeException::selfThrow("Interaction::saveInteractionToXML - object InteractionXML does not exist");

  /*
   * save the data of the Relation
   */
  // Main type of the relation: FirstOrder or Lagrangian
  RELATION::TYPES type = relation()->getType();
  // Subtype of the relation
  RELATION::SUBTYPES subType = relation()->getSubType();

  if (type == FirstOrder)
  {
    if (subType == NonLinearR)
      relation()->saveRelationToXML();
    else if (subType == LinearR)
      (boost::static_pointer_cast<FirstOrderLinearR>(relation()))->saveRelationToXML();
    else if (subType == LinearTIR)
      (boost::static_pointer_cast<FirstOrderLinearTIR>(relation()))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation subtype: " + subType);
  }
  else if (type == Lagrangian)
  {
    if (subType == LinearTIR)
      (boost::static_pointer_cast<LagrangianLinearTIR>(relation()))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Not yet implemented for relation subtype " + subType);
  }
  else
    RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation type: " + type);

  /*
   * save the data of the NonSmoothLaw
   */

  nslaw()->saveNonSmoothLawToXML();

}


void Interaction::getLeftUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
{

  unsigned int k = 0;
  unsigned int NumDS = 0;
  DSIterator itDS;
  // itDS = dynamicalSystemsBegin();

  itDS =    _involvedDS->begin();

  // look for ds and its position in G
  while (*itDS != ds && itDS != dynamicalSystemsEnd())
  {
    k += (*itDS)->getDim();
    itDS++;
    NumDS++;
  }

  // check dimension (1)
  if ((*itDS)->getDim() != UnitaryBlock->size(1))
    RuntimeException::selfThrow("Interaction::getLeftUnitaryBlockForDS(DS, UnitaryBlock, ...): inconsistent sizes between UnitaryBlock and DS");

  SP::SiconosMatrix originalMatrix;

  RELATION::TYPES relationType = relation()->getType();

  if (relationType == FirstOrder)
  {
    SP::FirstOrderR r = boost::static_pointer_cast<FirstOrderR> (relation());
    originalMatrix = r->jachx();
  }
  else if (relationType == Lagrangian)
  {
    SP::LagrangianR r = boost::static_pointer_cast<LagrangianR> (relation());
    originalMatrix = r->jachq();
  }
  else if (relationType == NewtonEuler)
  {
    SP::NewtonEulerR r = boost::static_pointer_cast<NewtonEulerR> (relation());
    originalMatrix = r->jachqT();
  }
  else
    RuntimeException::selfThrow("Interaction::getLeftUnitaryBlockForDS, not yet implemented for relations of type " + relationType);

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = 0; //_relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}
void Interaction::getLeftUnitaryBlockForDSProjectOnConstraints(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
{
  unsigned int k = 0;
  unsigned int NumDS = 0;
  DSIterator itDS;

  itDS = _involvedDS->begin();

  Type::Siconos dsType = Type::value(*ds);
  if (dsType != Type::NewtonEulerDS)
    RuntimeException::selfThrow("Interaction::getLeftUnitaryBlockForDSForProject- ds is not from NewtonEulerDS.");

  RELATION::TYPES relationType = relation()->getType();
  if (relationType != NewtonEuler)
    RuntimeException::selfThrow("Interaction::getLeftUnitaryBlockForDSForProject- relation is not from NewtonEulerR.");

  // look for ds and its position in G
  while (*itDS != ds && itDS != dynamicalSystemsEnd())
  {
    k += (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
    itDS++;
    NumDS++;
  }

  // check dimension (1)
  unsigned int sizeDS = (boost::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
  if (sizeDS != UnitaryBlock->size(1))
    RuntimeException::selfThrow("Interaction::getLeftUnitaryBlockForDSForProject(DS, UnitaryBlock, ...): inconsistent sizes between UnitaryBlock and DS");


  SP::SiconosMatrix originalMatrix;

  SP::NewtonEulerR r = boost::static_pointer_cast<NewtonEulerR> (relation());
  //proj_with_q originalMatrix = r->jachqProj();
  originalMatrix = r->jachq();

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = 0;//_relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}
void Interaction::getRightUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
{
  unsigned int k = 0;
  DSIterator itDS;
  itDS = _involvedDS->begin();

  // look for ds and its position in G
  while (*itDS != ds && itDS != dynamicalSystemsEnd())
  {
    k += (*itDS)->getDim();
    itDS++;
  }

  // check dimension (1)
  if ((*itDS)->getDim() != UnitaryBlock->size(0))
    RuntimeException::selfThrow("Interaction::getRightUnitaryBlockForDS(DS, UnitaryBlock, ...): inconsistent sizes between UnitaryBlock and DS");


  SP::SiconosMatrix originalMatrix; // Complete matrix, Relation member.
  RELATION::TYPES relationType = relation()->getType();

  if (relationType == FirstOrder)
  {
    originalMatrix = relation()->jacglambda();
  }
  else if (relationType == Lagrangian || relationType == NewtonEuler)
  {
    RuntimeException::selfThrow("Interaction::getRightUnitaryBlockForDS, call not permit " + relationType);
  }
  else
    RuntimeException::selfThrow("Interaction::getRightUnitaryBlockForDS, not yet implemented for relations of type " + relationType);

  if (! originalMatrix)
    RuntimeException::selfThrow("Interaction::getRightUnitaryBlockForDS(DS, UnitaryBlock, ...): the right unitaryBlock is a NULL pointer (miss matrix B or H or gradients ...in relation ?)");

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = k;
  subPos[1] = 0;//_relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}

void Interaction::getExtraUnitaryBlock(SP::SiconosMatrix UnitaryBlock) const
{
  // !!! Warning: we suppose that D is unitaryBlock diagonal, ie that
  // there is no coupling between Interaction through D !!!  Any
  // coupling between relations through D must be taken into account
  // thanks to the nslaw (by "increasing" its dimension).

  SP::SiconosMatrix D;
  //  if(relation()->getNumberOfJacobiansForH()>1)
  D = relation()->jachlambda();

  if (! D)
  {
    UnitaryBlock->zero();
    return; //ie no extra unitaryBlock
  }

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = 0;//_relativePosition;
  subPos[1] = 0;//_relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(D, UnitaryBlock, subDim, subPos);
}


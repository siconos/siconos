/* Siconos-Kernel, Copyright INRIA 2005-2012.
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
  _initialized(false), _id("undefined"), _number(0), _interactionSize(0),
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
  _initialized(false) , _id("none"), _number(newNumber), _interactionSize(nInter),
  _sizeOfDS(0), _sizeZ(0), _y(2), _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  _involvedDS->insert(ds); // Warning: insert pointer to DS!!

}
Interaction::Interaction(const string& newId, SP::DynamicalSystem ds,
                         int newNumber, int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false), _id(newId), _number(newNumber), _interactionSize(nInter),
  _sizeOfDS(0), _sizeZ(0), _y(2), _nslaw(newNSL),  _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  _involvedDS->insert(ds); // Warning: insert pointer to DS!!
}


Interaction::Interaction(DynamicalSystemsSet& dsConcerned, int newNumber, int nInter,
                         SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false) , _id("none"), _number(newNumber), _interactionSize(nInter),
  _sizeOfDS(0), _sizeZ(0), _y(2), _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  DSIterator itDS;
  for (itDS = dsConcerned.begin(); itDS != dsConcerned.end(); ++itDS)
    _involvedDS->insert(*itDS); // Warning: insert pointers to DS!!
}

Interaction::Interaction(const string& newId, DynamicalSystemsSet& dsConcerned, int newNumber,
                         int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel):
  _initialized(false) , _id(newId), _number(newNumber), _interactionSize(nInter),  _sizeOfDS(0), _sizeZ(0),
  _y(2),  _nslaw(newNSL), _relation(newRel)
{
  _involvedDS.reset(new DynamicalSystemsSet());
  DSIterator itDS;
  for (itDS = dsConcerned.begin(); itDS != dsConcerned.end(); ++itDS)
    _involvedDS->insert(*itDS); // Warning: insert pointers to DS!!
}

/* initialisation with empty set */
Interaction::Interaction(int nInter, SP::NonSmoothLaw newNSL, SP::Relation newRel, int newNumber):
  _initialized(false), _number(newNumber), _interactionSize(nInter), _sizeOfDS(0), _sizeZ(0),
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

    // compute number of relations.

    if (_interactionSize != nslaw()->size())
    {
      RuntimeException::selfThrow("Interaction::initialize() - _interactionSize != nslaw()->size() . Obsolete !");
    }


    initData();
    initializeMemory();
    _relation->initialize(*this);

    if (_steps > 1) // Multi--step methods
    {
      // Comoyte the old Values of Output with stored values in Memory
      for (unsigned int k = 0; k < _steps - 1; k++)
      {
        /** ComputeOutput to fill the Memory
         * We assume the state x is stored in xMemory except for the  initial
         * condition which has not been swap yet.
         */
        //        relation()->LinkDataFromMemory(k);
        for (unsigned int i = 0; i < _upperLevelForOutput + 1; ++i)
        {
          computeOutput(t0, i);
          _yMemory[i]->swap(_y[i]);
        }
      }
      //      relation()->LinkData();
    }

    // Compute y values for t0
    for (unsigned int i = 0; i < _upperLevelForOutput + 1; ++i)
    {
      computeOutput(t0, i);
    }
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
  // and they are computed in Simulation::computeLevelsForInputAndOutput

  // Memory allocation for y and lambda

  //  assert(_upperLevelForOutput >=0);
  assert(_upperLevelForOutput >= _lowerLevelForOutput);
  //  assert(_upperLevelForInput >=0);
  assert(_upperLevelForInput >= _lowerLevelForInput);

  // in order to simplify we size from 0 to _upperLevelForXXX
  _y.resize(_upperLevelForOutput + 1) ;
  _yOld.resize(_upperLevelForOutput + 1);
  _y_k.resize(_upperLevelForOutput + 1);
  _lambda.resize(_upperLevelForInput + 1);
  _lambdaOld.resize(_upperLevelForInput + 1);

  _yMemory.resize(_upperLevelForOutput + 1);



  // get the dimension of the non smooth law, ie the size of an Interaction blocks (one per relation)
  unsigned int nslawSize = nslaw()->size();
  _Residuy.reset(new SiconosVector(nslawSize));
  _h_alpha.reset(new SiconosVector(nslawSize));
  _workYp.reset(new SiconosVector(nslawSize));

  for (unsigned int i = _lowerLevelForOutput ;
       i < _upperLevelForOutput + 1 ;
       i++)
  {
    _y[i].reset(new SiconosVector(nslawSize));
    _yOld[i].reset(new SiconosVector(nslawSize));
    _y_k[i].reset(new SiconosVector(nslawSize));
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
    _lambda[i].reset(new SiconosVector(nslawSize));
    _lambdaOld[i].reset(new SiconosVector(nslawSize));
    _lambdaOld[i]->zero();
  }

}

void Interaction::initData()
{
  RELATION::TYPES relationType = _relation->getType();
  if (relationType == FirstOrder)
    initDataFirstOrder();
  else if (relationType == Lagrangian)
    initDataLagrangian();
  else if (relationType == NewtonEuler)
    initDataNewtonEuler();
  else
    RuntimeException::selfThrow("Interaction::initData unknown initialization procedure for \
        a relation of type: " + relationType);

}

// It could be interesting to make Interaction a pure virtual class and to derive 3
// classes, one for each type of relation
void Interaction::initDataFirstOrder()
{
  // Get the DS concerned by the interaction of this relation
  _data[free].reset(new BlockVector());
  _data[x].reset(new BlockVector()); // displacements
  _data[xq].reset(new BlockVector());
  _data[deltax].reset(new BlockVector()); // displacements
  _data[z].reset(new BlockVector());
  _data[r].reset(new BlockVector());
  _data[residu_r].reset(new BlockVector());
  _data[ds_xp].reset(new BlockVector());
  _data[g_alpha].reset(new BlockVector());

  for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
  {
    // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
    FirstOrderNonLinearDS& ds = static_cast<FirstOrderNonLinearDS&>(**it);
    _data[free]->insertPtr(ds.workFree());
    _data[x]->insertPtr(ds.x());
    _data[xq]->insertPtr(ds.xq());
    _data[deltax]->insertPtr(ds.getWorkVector(DynamicalSystem::local_buffer));
    _data[z]->insertPtr(ds.z());
    _data[r]->insertPtr(ds.r());
    _data[residu_r]->insertPtr(ds.residur());
    _data[ds_xp]->insertPtr(ds.xp());
    _data[g_alpha]->insertPtr(ds.gAlpha());

  }
}

void Interaction::initDataLagrangian()
{
  _data[free].reset(new BlockVector());
  _data[q0].reset(new BlockVector()); // displacement
  _data[q1].reset(new BlockVector()); // velocity
  _data[q2].reset(new BlockVector()); // acceleration
  _data[z].reset(new BlockVector()); // z vector
  _data[p0].reset(new BlockVector());
  _data[p1].reset(new BlockVector());
  _data[p2].reset(new BlockVector());

  SP::LagrangianDS lds;
  for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
  {
    // check dynamical system type
    assert((Type::value(**it) == Type::LagrangianLinearTIDS ||
            Type::value(**it) == Type::LagrangianDS));

    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = std11::static_pointer_cast<LagrangianDS> (*it);

    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    _data[free]->insertPtr(lds->workFree());
    _data[q0]->insertPtr(lds->q());
    _data[q1]->insertPtr(lds->velocity());
    _data[q2]->insertPtr(lds->acceleration());
    _data[z]->insertPtr(lds->z());

    // Put NonsmoothInput _p of each DS into a block. (Pointers links, no copy!!)
    for (unsigned int k = _lowerLevelForInput;
         k < _upperLevelForInput + 1; k++)
    {
      assert(lds->p(k));
      assert(_data[p0 + k]);
      _data[p0 + k]->insertPtr(lds->p(k));
    }
  }
}

void Interaction::initDataNewtonEuler()
{
  DSIterator it;
  _data[free].reset(new BlockVector());
  _data[q0].reset(new BlockVector()); // displacement
  _data[velo].reset(new BlockVector()); // velocity
  _data[deltaq].reset(new BlockVector());
  _data[q1].reset(new BlockVector()); // qdot
  //  data[q2].reset(new BlockVector()); // acceleration
  _data[z].reset(new BlockVector()); // z vector
  _data[p0].reset(new BlockVector());
  _data[p1].reset(new BlockVector());
  _data[p2].reset(new BlockVector());
  SP::NewtonEulerDS lds;
  unsigned int sizeForAllxInDs = 0;
  for (it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
  {
    // check dynamical system type
    assert((Type::value(**it) == Type::NewtonEulerDS) && "NewtonEulerR::initialize failed, not implemented for dynamical system of that type.\n");

    // convert vDS systems into NewtonEulerDS and put them in vLDS
    lds = std11::static_pointer_cast<NewtonEulerDS> (*it);
    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    _data[free]->insertPtr(lds->workFree());
    _data[q0]->insertPtr(lds->q());
    _data[velo]->insertPtr(lds->velocity());
    _data[deltaq]->insertPtr(lds->deltaq());
    _data[q1]->insertPtr(lds->dotq());
    //    data[q2]->insertPtr( lds->acceleration());
    if (lds->p(0))
      _data[p0]->insertPtr(lds->p(0));
    if (lds->p(1))
      _data[p1]->insertPtr(lds->p(1));
    if (lds->p(2))
      _data[p2]->insertPtr(lds->p(2));

    _data[z]->insertPtr(lds->z());
    sizeForAllxInDs += lds->p(1)->size();
  }
}

void Interaction::LinkDataFromMemory(unsigned int memoryLevel)
{
  if (_relation->getType() == Lagrangian)
    LinkDataFromMemoryLagrangian(memoryLevel);
  else
    RuntimeException::selfThrow("Interaction::LinkDataFromMemory: not yet implemented for Relation of type " + _relation->getType());
}

void Interaction::LinkDataFromMemoryLagrangian(unsigned int memoryLevel)
{

  _data[q0].reset(new BlockVector()); // displacement
  _data[q1].reset(new BlockVector()); // velocity
  _data[q2].reset(new BlockVector()); // acceleration
  _data[z].reset(new BlockVector()); // z vector
  _data[p0].reset(new BlockVector());
  _data[p1].reset(new BlockVector());
  _data[p2].reset(new BlockVector());

  SP::LagrangianDS lds;
  for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
  {
    // check dynamical system type
    assert((Type::value(**it) == Type::LagrangianLinearTIDS ||
            Type::value(**it) == Type::LagrangianDS));
    // convert vDS systems into LagrangianDS and put them in vLDS
    lds = std11::static_pointer_cast<LagrangianDS> (*it);

    // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
    _data[q0]->insertPtr(lds->qMemory()->getSiconosVector(memoryLevel));
    _data[q1]->insertPtr(lds->velocityMemory()->getSiconosVector(memoryLevel));


    // Do nothing for the remaining of data since there are no Memory
    // An access to the content ofdata[q2] based on a link on Memory
    //must throw an exeption

  }
}
// --- GETTERS/SETTERS ---

void Interaction::setY(const VectorOfVectors& newVector)
{

  unsigned int size = newVector.size();

  _y.clear();
  _y.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _y[i].reset(new SiconosVector(*(newVector[i]))); // -> copy !
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
    _y[index].reset(new SiconosVector(newY));
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

  _y[index] = newY;
}

void Interaction::setYOld(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();
  _yOld.clear();
  _yOld.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _yOld[i].reset(new SiconosVector(*(newVector[i]))); // -> copy !
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
    _yOld[index].reset(new SiconosVector(newYOld));
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

  _yOld[index] = newYOld;
}

void Interaction::setLambda(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();
  _lambda.clear();
  _lambda.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _lambda[i].reset(new SiconosVector(*(newVector[i]))); // -> copy !
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
    _lambda[index].reset(new SiconosVector(newLambda));
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

  _lambda[index] = newLambda;
}

void Interaction::setLambdaOld(const VectorOfVectors& newVector)
{
  unsigned int size = newVector.size();

  // clear lambdaOld
  _lambdaOld.clear();
  _lambdaOld.resize(size);

  for (unsigned int i = 0; i < size; i++)
    _lambdaOld[i].reset(new SiconosVector(*(newVector[i]))); // -> copy !
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
    _lambdaOld[index].reset(new SiconosVector(newLambdaOld));
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

  _lambdaOld[index] = newLambdaOld;
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
    assert(_y[i]);
    assert(_yOld[i]);

    *(_yOld[i]) = *(_y[i]) ;
  }

  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1  ; i++)
  {
    assert(_lambdaOld[i]);
    assert(_lambda[i]);
    *(_lambdaOld[i]) = *(_lambda[i]);
  }
}

void Interaction::swapTimeStepInMemory()
{
  // i corresponds to the derivative number and j the relation number.
  for (unsigned int i = 0; i < _y.size() ; i++)
  {

    *(_y_k[i]) = *(_y[i]) ;

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
  relation()->computeOutput(time, *this, level);
}

void Interaction::computeInput(double time, unsigned int level)
{
  relation()->computeInput(time, *this, level);
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
      (std11::static_pointer_cast<FirstOrderLinearR>(relation()))->saveRelationToXML();
    else if (subType == LinearTIR)
      (std11::static_pointer_cast<FirstOrderLinearTIR>(relation()))->saveRelationToXML();
    else
      RuntimeException::selfThrow("Interaction::saveInteractionToXML - Unknown relation subtype: " + subType);
  }
  else if (type == Lagrangian)
  {
    if (subType == LinearTIR)
      (std11::static_pointer_cast<LagrangianLinearTIR>(relation()))->saveRelationToXML();
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


void Interaction::getLeftInteractionBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix InteractionBlock) const
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
  if ((*itDS)->getDim() != InteractionBlock->size(1))
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDS(DS, InteractionBlock, ...): inconsistent sizes between InteractionBlock and DS");

  SP::SiconosMatrix originalMatrix;

  RELATION::TYPES relationType = relation()->getType();

  if (relationType == FirstOrder)
  {
    SP::FirstOrderR r = std11::static_pointer_cast<FirstOrderR> (relation());
    originalMatrix = r->jachx();
  }
  else if (relationType == Lagrangian)
  {
    SP::LagrangianR r = std11::static_pointer_cast<LagrangianR> (relation());
    originalMatrix = r->jachq();
  }
  else if (relationType == NewtonEuler)
  {
    SP::NewtonEulerR r = std11::static_pointer_cast<NewtonEulerR> (relation());
    originalMatrix = r->jachqT();
  }
  else
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDS, not yet implemented for relations of type " + relationType);

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  Index subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  Index subPos(4);
  subPos[0] = 0; //_relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, InteractionBlock, subDim, subPos);
}
void Interaction::getLeftInteractionBlockForDSProjectOnConstraints(SP::DynamicalSystem ds, SP::SiconosMatrix InteractionBlock) const
{
  unsigned int k = 0;
  unsigned int NumDS = 0;
  DSIterator itDS;

  itDS = _involvedDS->begin();

  Type::Siconos dsType = Type::value(*ds);
  if (dsType != Type::NewtonEulerDS)
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDSForProject- ds is not from NewtonEulerDS.");

  RELATION::TYPES relationType = relation()->getType();
  if (relationType != NewtonEuler)
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDSForProject- relation is not from NewtonEulerR.");

  // look for ds and its position in G
  while (*itDS != ds && itDS != dynamicalSystemsEnd())
  {
    k += (std11::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
    itDS++;
    NumDS++;
  }

  // check dimension (1)
  unsigned int sizeDS = (std11::static_pointer_cast<NewtonEulerDS>(ds))->getqDim();
  if (sizeDS != InteractionBlock->size(1))
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDSForProject(DS, InteractionBlock, ...): inconsistent sizes between InteractionBlock and DS");


  SP::SiconosMatrix originalMatrix;

  SP::NewtonEulerR r = std11::static_pointer_cast<NewtonEulerR> (relation());
  //proj_with_q originalMatrix = r->jachqProj();
  originalMatrix = r->jachq();

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  Index subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  Index subPos(4);
  subPos[0] = 0;//_relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, InteractionBlock, subDim, subPos);
}
void Interaction::getRightInteractionBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix InteractionBlock) const
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
  if ((*itDS)->getDim() != InteractionBlock->size(0))
    RuntimeException::selfThrow("Interaction::getRightInteractionBlockForDS(DS, InteractionBlock, ...): inconsistent sizes between InteractionBlock and DS");


  SP::SiconosMatrix originalMatrix; // Complete matrix, Relation member.
  RELATION::TYPES relationType = relation()->getType();

  if (relationType == FirstOrder)
  {
    originalMatrix = std11::static_pointer_cast<FirstOrderR>(relation())->jacglambda();
  }
  else if (relationType == Lagrangian || relationType == NewtonEuler)
  {
    RuntimeException::selfThrow("Interaction::getRightInteractionBlockForDS, call not permit " + relationType);
  }
  else
    RuntimeException::selfThrow("Interaction::getRightInteractionBlockForDS, not yet implemented for relations of type " + relationType);

  if (! originalMatrix)
    RuntimeException::selfThrow("Interaction::getRightInteractionBlockForDS(DS, InteractionBlock, ...): the right interactionBlock is a NULL pointer (miss matrix B or H or gradients ...in relation ?)");

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  Index subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  Index subPos(4);
  subPos[0] = k;
  subPos[1] = 0;//_relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, InteractionBlock, subDim, subPos);
}

void Interaction::getExtraInteractionBlock(SP::SiconosMatrix InteractionBlock) const
{
  // !!! Warning: we suppose that D is interactionBlock diagonal, ie that
  // there is no coupling between Interaction through D !!!  Any
  // coupling between relations through D must be taken into account
  // thanks to the nslaw (by "increasing" its dimension).

  SP::SiconosMatrix D = relation()->jachlambda();
  //  if(relation()->getNumberOfJacobiansForH()>1)

  if (!D)
  {
    InteractionBlock->zero();
    return; //ie no extra interactionBlock
  }

  // copy sub-interactionBlock of originalMatrix into InteractionBlock
  // dim of the sub-interactionBlock
  Index subDim(2);
  subDim[0] = InteractionBlock->size(0);
  subDim[1] = InteractionBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in InteractionBlock
  Index subPos(4);
  subPos[0] = 0;//_relativePosition;
  subPos[1] = 0;//_relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(D, InteractionBlock, subDim, subPos);
}

void Interaction::computeResiduY(const double time)
{
  //Residu_y = y_alpha_k+1 - H_alpha;
  *_Residuy = *_h_alpha;
  scal(-1, *_Residuy, *_Residuy);

  (*_Residuy) += *(y(0));

}

void Interaction::computeResiduR(const double time)
{
  //Residu_r = r_alpha_k+1 - g_alpha;
  *_data[residu_r] = *_data[r];
  *_data[residu_r] -= *_data[g_alpha];

  // std::cout<< "Interaction::computeResiduR(const double time)" << std::endl;
  // std::cout<< "_data[r] = " << std::endl ;
  // _data[r]->display();
  // std::cout<< "_data[g_alpha] = " << std::endl ;
  // _data[g_alpha]->display();
  // std::cout<< "_data[residu_r] = " << std::endl ;
  // _data[residu_r]->display();

  //RuntimeException::selfThrow("Interaction::computeResiduR do not use this function");
}

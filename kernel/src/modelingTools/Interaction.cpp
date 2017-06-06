/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2016 INRIA.
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
#include <assert.h>
#include <iostream>
//#define DEBUG_BEGIN_END_ONLY
// #define DEBUG_STDOUT
// #define DEBUG_NOCOLOR
// #define DEBUG_MESSAGES
#include "debug.h"

#include "Interaction.hpp"
#include "RelationTypes.hpp"
#include "ComplementarityConditionNSL.hpp"
#include "RelayNSL.hpp"
#include "NewtonImpactNSL.hpp"
#include "NewtonImpactFrictionNSL.hpp"
#include "DynamicalSystem.hpp"

#include "LagrangianDS.hpp"

#include "FirstOrderR.hpp"
#include "LagrangianR.hpp"
#include "NewtonEulerR.hpp" // ??
#include "NewtonEulerDS.hpp" // ??

#include "BlockVector.hpp"
#include "FirstOrderNonLinearDS.hpp"

#include "SimulationGraphs.hpp"

using namespace std;
using namespace RELATION;


unsigned int Interaction::__count = 0;

// --- CONSTRUCTORS ---
struct Interaction::_setLevels : public SiconosVisitor
{
  /* we set the _lowerLevelForOutput, _upperLevelForOutput,
     _lowerLevelForOutput, _upperLevelForOutput
     w.r.t to the choice of the nslaw and the relation
   */
  using SiconosVisitor::visit;

  Interaction* _interaction;

  _setLevels(Interaction * inter) :
    _interaction(inter) {};

  void visit(const ComplementarityConditionNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    RELATION::SUBTYPES subType = _interaction->relation()->getSubType();

    if (relationType == FirstOrder)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(0);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(0);
    }
    else if (relationType == Lagrangian && subType == CompliantLinearTIR )
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(0);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(0);
    }
    else
    {
      RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    };
  }

  void visit(const RelayNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == FirstOrder)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(0);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(0);
    }
    else if (relationType == Lagrangian || relationType == NewtonEuler)
    {
      // For friction
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(1);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(1);
    }
    else
    {
      RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    };
  }

 void visit(const NormalConeNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == FirstOrder)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(0);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(0);
    }
    else
    {
      RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    };
  }

  void visit(const MixedComplementarityConditionNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == FirstOrder)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(0);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(0);
    }
    else
    {
      RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    };
  }

  void visit(const EqualityConditionNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == Lagrangian || relationType == NewtonEuler)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(1);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(1);

    }
    else
    {
	RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    }
    ;
  }
  void visit(const NewtonImpactNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == Lagrangian || relationType == NewtonEuler)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(1);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(1);

    }
    else
    {
	RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    }
  }

  void visit(const NewtonImpactFrictionNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == Lagrangian || relationType == NewtonEuler)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(1);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(1);

    }
    else
    {
	RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    }
  }
  void visit(const MultipleImpactNSL& nslaw)
  {
    RELATION::TYPES relationType = _interaction->relation()->getType();
    if (relationType == Lagrangian || relationType == NewtonEuler)
    {
      _interaction->setLowerLevelForOutput(0);
      _interaction->setUpperLevelForOutput(1);

      _interaction->setLowerLevelForInput(0);
      _interaction->setUpperLevelForInput(1);

    }
    else
    {
	RuntimeException::selfThrow("Interaction::_setLevels::visit - unknown relation type for the nslaw ");
    }
  }
};


void Interaction::reset()
{
  // Check levels values and
  // resize all containers-like attributes according to these levels.

  // This function must be called at the first instanciation of
  // an interaction (in __init) and may be called by simulation and/or
  // OSI if levels are updated.

  //  assert(_upperLevelForOutput >=0);
  assert(_upperLevelForOutput >= _lowerLevelForOutput);
  //  assert(_upperLevelForInput >=0);
  assert(_upperLevelForInput >= _lowerLevelForInput);

  // --  Memory allocation for y and lambda --
   // in order to simplify we size from 0 to _upperLevelForXXX
  _y.resize(_upperLevelForOutput + 1) ;
  _lambda.resize(_upperLevelForInput + 1);

  // -- Memory allocation for buffers (OSI related ! Must be moved to the graph)
  _yOld.resize(_upperLevelForOutput + 1);
  _y_k.resize(_upperLevelForOutput + 1);
  _lambdaOld.resize(_upperLevelForInput + 1);

  // get the dimension of the non smooth law, ie the size of an Interaction blocks (one per relation)
  unsigned int nslawSize = _nslaw->size();

  for (unsigned int i = _lowerLevelForOutput ;
       i < _upperLevelForOutput + 1 ;
       i++)
  {
    _y[i].reset(new SiconosVector(nslawSize));
    _yOld[i].reset(new SiconosVector(nslawSize));
    _y_k[i].reset(new SiconosVector(nslawSize));

    _y[i]->zero();
    _yOld[i]->zero();
    _y_k[i]->zero();
  }

  
  for (unsigned int i = _lowerLevelForInput ;
       i < _upperLevelForInput + 1 ;
       i++)
  {
    DEBUG_PRINTF("Interaction::initializeMemory(). _lambda[%i].reset()\n",i)
      _lambda[i].reset(new SiconosVector(nslawSize));
    _lambdaOld[i].reset(new SiconosVector(nslawSize));
    _lambdaOld[i]->zero();
  }
}


void Interaction::__init()
{
  // -- Delagated constructor --
  // i.e. what should be done when (and only there) the interaction
  // is instanciated.
  // Other operations (like levels review and y, lambda resizing)
  // occur in reset function, potentially called during
  // simulation phase (in OSI indeed).

  assert(_relation && "Interaction::__init failed, relation() == NULL");
  assert(_nslaw && "Interaction::__inits, non smooth law == NULL");

  // -- Set upper/lower levels, according to the nslaw --
  std11::shared_ptr<_setLevels> setLevels;
  setLevels.reset(new _setLevels(this));
  _nslaw->accept(*(setLevels.get()));

  // Ensure consistency between interaction and nslaw sizes
  if (_interactionSize != _nslaw->size())
    RuntimeException::selfThrow("Interaction::__init - Nonsmooth law and relation are not consistent.");

  // Check levels and resize attributes (y, lambda ...) if needed.
  reset();
}

Interaction::Interaction(SP::NonSmoothLaw NSL,
                         SP::Relation rel):
  _number(__count++), _interactionSize(NSL->size()),
  _sizeOfDS(0), _has2Bodies(false), _y(2),  _nslaw(NSL), _relation(rel)
{
  __init();
}




// void Interaction::setDSLinkAndWorkspace(InteractionProperties& interProp,
// 					DynamicalSystem& ds1, VectorOfVectors& workV1,
// 					DynamicalSystem& ds2, VectorOfVectors& workV2)
// {
//   DEBUG_BEGIN("Interaction::setDSLinkAndWorkspace(...)\n");

//   VectorOfBlockVectors& DSlink = *interProp.DSlink;
//   VectorOfVectors& workVInter = *interProp.workVectors;
//   VectorOfSMatrices& workMInter = *interProp.workMatrices;

//   initData(DSlink);
//   // Initialize interaction work vectors, depending on Dynamical systems
//   // linked to the interaction.

//   initDSData(ds1, DSlink);

//   if(&ds1 != &ds2)
//     {
//       DEBUG_PRINT("ds1 != ds2\n");
//       DEBUG_PRINTF("ds1 number %i", ds1.number())
//       DEBUG_PRINTF("ds2 number %i", ds2.number())
//       initDSData(ds2, DSlink);
//     }

//   bool computeResidu = _relation->requireResidu();

//   // Relation initializes the work vectors and matrices
//   _relation->initialize(*this, DSlink, workVInter, workMInter);

//   if (computeResidu)
//     {
//       RELATION::TYPES relationType = _relation->getType();
//       if (relationType == FirstOrder)
// 	{
// 	  if (!workVInter[FirstOrderR::g_alpha])
// 	    workVInter[FirstOrderR::g_alpha].reset(new SiconosVector(_sizeOfDS));
// 	  if (!workVInter[FirstOrderR::vec_residuR])
// 	    workVInter[FirstOrderR::vec_residuR].reset(new SiconosVector(_sizeOfDS));
// 	}
//       else if (relationType == Lagrangian)
//         RuntimeException::selfThrow("Interaction::initialize() - computeResiduR for LagrangianR is not implemented");
//       else if (relationType == NewtonEuler)
//         RuntimeException::selfThrow("Interaction::initialize() - computeResiduR for NewtonEulerR is not implemented");
//     }

//   DEBUG_END(" Interaction::setDSLinkAndWorkspace(...)\n");
// }



void Interaction::initialize_ds_links(InteractionProperties& interaction_properties, DynamicalSystem& ds1,
				      DynamicalSystem& ds2)
{
  // Initialize DSlink property

  interaction_properties.DSlink.reset(new VectorOfBlockVectors);
  // Get (from graph) DSLink property.
  // This container of vectors is supposed to handle
  // pointer links to dynamical system(s) attributes
  // that may be used to compute input and output.
  // The list of potential keys depends on the relation type
  // and is defined in an enum, in XXR.hpp, XX being the relation type
  // (Lagrangian, NewtonEuler or FirstOrder)
  VectorOfBlockVectors& DSlink = *interaction_properties.DSlink;
  
  // The dynamical systems linked to the interaction (2 at most, ds2 may be equal to ds1).
  RELATION::TYPES relationType = _relation->getType();

  if (relationType == FirstOrder)
    __initDataFirstOrder(DSlink, ds1, ds2);
  
  else if (relationType == Lagrangian)
    __initDataLagrangian(DSlink, ds1, ds2);

  else if (relationType == NewtonEuler)
    __initDataNewtonEuler(DSlink, ds1, ds2);
  else
    RuntimeException::selfThrow("Interaction::initData unknown initialization procedure for \
        a relation of type: " + relationType);

  // -- Stage 2 : create buffers (in the graph) that will be used for relation/interaction internal operations --
  // Relation initializes the work vectors and matrices
  //
  interaction_properties.workVectors.reset(new VectorOfVectors);
  interaction_properties.workMatrices.reset(new VectorOfSMatrices);
  VectorOfVectors& workVInter = *interaction_properties.workVectors;
  VectorOfSMatrices& workMInter = *interaction_properties.workMatrices;
  _relation->initialize(*this, DSlink, workVInter, workMInter);
}



// Initialize and InitializeMemory are separated in two functions
// since we need to know the relative degree to know
// "numberOfDerivatives", while numberOfRelations and the size of the
// non smooth law are required inputs to compute the relative degree.
void Interaction::initializeMemory(bool computeResidu, unsigned int steps)
{

  DEBUG_PRINT("Interaction::initializeMemory() \n");
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

  _yMemory.resize(_upperLevelForOutput + 1);
  _lambdaMemory.resize(_upperLevelForInput + 1);
  unsigned int nslawSize = _nslaw->size();

  for (unsigned int i = _lowerLevelForOutput ; i < _upperLevelForOutput + 1 ; i++)
    _yMemory[i].reset(new SiconosMemory(steps, nslawSize));

  for (unsigned int i = _lowerLevelForInput ; i < _upperLevelForInput + 1 ; i++)
  {
    DEBUG_PRINTF("Interaction::initializeMemory(). _lambda[%i].reset()\n",i)
    _lambdaMemory[i].reset(new SiconosMemory(steps, nslawSize));
  }


}

void Interaction::resetAllLambda()
{
   for (unsigned int i = _lowerLevelForInput ;
       i < _upperLevelForInput + 1 ;
       i++)
  {
    if (_lambda[i])
      _lambda[i]->zero();
  }

}


void Interaction::resetLambda(unsigned int level)
{
  if (_lambda[level])
    _lambda[level]->zero();
}


// It could be interesting to make Interaction a pure virtual class and to derive 3
// classes, one for each type of relation
 void Interaction::__initDataFirstOrder(VectorOfBlockVectors& DSlink, DynamicalSystem& ds1, DynamicalSystem& ds2)
{

  DSlink.resize(FirstOrderR::DSlinkSize);
  DSlink[FirstOrderR::x].reset(new BlockVector());
  DSlink[FirstOrderR::r].reset(new BlockVector());
  DSlink[FirstOrderR::z].reset(new BlockVector());
  __initDSDataFirstOrder(ds1, DSlink);
  if(&ds1 != &ds2)
    __initDSDataFirstOrder(ds2, DSlink);
 
}

void Interaction::__initDSDataFirstOrder(DynamicalSystem& ds, VectorOfBlockVectors& DSlink)
{
  // Put x/r ... of each DS into a block. (Pointers links, no copy!!)
  FirstOrderNonLinearDS& lds = static_cast<FirstOrderNonLinearDS&>(ds);
  DSlink[FirstOrderR::x]->insertPtr(lds.x());
  DSlink[FirstOrderR::r]->insertPtr(lds.r());
  DSlink[FirstOrderR::z]->insertPtr(lds.z());
}

void Interaction::__initDataLagrangian(VectorOfBlockVectors& DSlink, DynamicalSystem& ds1, DynamicalSystem& ds2)
{

  DEBUG_PRINT("Interaction::initDataLagrangian()\n");
  DSlink.resize(LagrangianR::DSlinkSize);
  DSlink[LagrangianR::q0].reset(new BlockVector()); // displacement
  DSlink[LagrangianR::q1].reset(new BlockVector()); // velocity
  DSlink[LagrangianR::q2].reset(new BlockVector()); // acceleration
  DSlink[LagrangianR::p0].reset(new BlockVector());
  DSlink[LagrangianR::p1].reset(new BlockVector());
  DSlink[LagrangianR::p2].reset(new BlockVector());
  DSlink[LagrangianR::z].reset(new BlockVector());
  __initDSDataLagrangian(ds1, DSlink);
  if(&ds1 != &ds2)
    __initDSDataLagrangian(ds2, DSlink);

}

void Interaction::__initDSDataLagrangian(DynamicalSystem& ds, VectorOfBlockVectors& DSlink)
{
  // check dynamical system type
  assert((Type::value(ds) == Type::LagrangianLinearTIDS ||
          Type::value(ds) == Type::LagrangianDS ||
          Type::value(ds) == Type::LagrangianLinearDiagonalDS));

  // convert vDS systems into LagrangianDS and put them in vLDS
  LagrangianDS& lds = static_cast<LagrangianDS&> (ds);

  // Put q, velocity and acceleration of each DS into a block. (Pointers links, no copy!!)
  DSlink[LagrangianR::q0]->insertPtr(lds.q());

  DSlink[LagrangianR::q1]->insertPtr(lds.velocity());
  if(lds.acceleration())
    DSlink[LagrangianR::q2]->insertPtr(lds.acceleration());
  DSlink[LagrangianR::z]->insertPtr(lds.z());

  for (unsigned int k = 0; k < 3; k++)
  {
    if(lds.p(k) && DSlink[LagrangianR::p0 + k])
      DSlink[LagrangianR::p0 + k]->insertPtr(lds.p(k));
  }
}

void Interaction::__initDataNewtonEuler(VectorOfBlockVectors& DSlink, DynamicalSystem& ds1, DynamicalSystem& ds2)
{
  DEBUG_BEGIN("Interaction::initDataNewtonEuler(VectorOfBlockVectors& DSlink)\n");
  DSlink.resize(NewtonEulerR::DSlinkSize);
  //DSlink[NewtonEulerR::xfree].reset(new BlockVector());
  DSlink[NewtonEulerR::q0].reset(new BlockVector()); // displacement
  DSlink[NewtonEulerR::velocity].reset(new BlockVector()); // velocity
//  DSlink[NewtonEulerR::deltaq].reset(new BlockVector());
  DSlink[NewtonEulerR::dotq].reset(new BlockVector()); // qdot
  //  data[NewtonEulerR::q2].reset(new BlockVector()); // acceleration
  DSlink[NewtonEulerR::z].reset(new BlockVector()); // z vector
  DSlink[NewtonEulerR::p0].reset(new BlockVector());
  DSlink[NewtonEulerR::p1].reset(new BlockVector());
  DSlink[NewtonEulerR::p2].reset(new BlockVector());
  DEBUG_END("Interaction::initDataNewtonEuler(VectorOfBlockVectors& DSlink)\n");
  __initDSDataNewtonEuler(ds1, DSlink);
  if(&ds1 != &ds2)
    __initDSDataNewtonEuler(ds2, DSlink);

}

void Interaction::__initDSDataNewtonEuler(DynamicalSystem& ds, VectorOfBlockVectors& DSlink)
{
  DEBUG_BEGIN("Interaction::initDSDataNewtonEuler(DynamicalSystem& ds, VectorOfBlockVectors& DSlink)\n");
  // check dynamical system type
  assert((Type::value(ds) == Type::NewtonEulerDS) && "Interaction initDSData failed, not implemented for dynamical system of that type.\n");

  // convert vDS systems into NewtonEulerDS and put them in vLDS
  NewtonEulerDS& neds = static_cast<NewtonEulerDS&>(ds);
  // Put q/velocity/acceleration of each DS into a block. (Pointers links, no copy!!)
  DSlink[NewtonEulerR::q0]->insertPtr(neds.q());
  DSlink[NewtonEulerR::velocity]->insertPtr(neds.twist());
  //  DSlink[NewtonEulerR::deltaq]->insertPtr(neds.deltaq());
  DSlink[NewtonEulerR::dotq]->insertPtr(neds.dotq());
  //    data[NewtonEulerR::q2]->insertPtr( neds.acceleration());
  if (neds.p(0))
      DSlink[NewtonEulerR::p0]->insertPtr(neds.p(0));
  if (neds.p(1))
    DSlink[NewtonEulerR::p1]->insertPtr(neds.p(1));
  if (neds.p(2))
    DSlink[NewtonEulerR::p2]->insertPtr(neds.p(2));

  DSlink[NewtonEulerR::z]->insertPtr(neds.z());
  DEBUG_END("Interaction::initDSDataNewtonEuler(DynamicalSystem& ds, VectorOfBlockVectors& DSlink)\n");

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


// --- OTHER FUNCTIONS ---

void Interaction::swapInOldVariables()
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

void Interaction::swapInMemory()
{
  // i corresponds to the derivative number and j the relation number.
  for (unsigned int  i = _lowerLevelForOutput; i < _upperLevelForOutput + 1 ; i++)
  {
    *(_y_k[i]) = *(_y[i]) ;
    _yMemory[i]->swap(*_y[i]);
  }

  for (unsigned int i = _lowerLevelForInput; i < _upperLevelForInput + 1  ; i++)
  {
    _lambdaMemory[i]->swap(*_lambda[i]);
  }

}

void Interaction::display() const
{
  std::cout << "======= Interaction display number " << _number <<" =======" <<std::endl;

  cout << "| lowerLevelForOutput : " << _lowerLevelForOutput << endl;
  cout << "| upperLevelForOutput : " << _upperLevelForOutput << endl;
  cout << "| lowerLevelForInput : " << _lowerLevelForInput << endl;
  cout << "| upperLevelForInput : " << _upperLevelForInput << endl;
  cout << "| interactionSize : " << _interactionSize << endl;
  cout << "| _sizeOfDS : " << _sizeOfDS << endl;

  cout << "| "  ; _relation->display();
  for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++)
    {

      std::cout << "| y[" << i  << "] : ";
      if (_y[i])
	{
	  if (_y[i]->size() >= 5) std::cout <<std::endl;
	  _y[i]->display();
	}
      else std::cout << "->NULL" <<std::endl;
    }
  for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++)
    {
      std::cout << "| yOld[" << i  << "] : ";
      if (_yOld[i])
	{
	  if (_yOld[i]->size() >= 5) std::cout <<std::endl;
	  _yOld[i]->display();
	}
      else std::cout << "->NULL" <<std::endl;
    }
  for (unsigned int i = 0; i < _upperLevelForOutput + 1; i++)
    {
      std::cout << "| y_k[" << i  << "] : ";
      if (_y_k[i])
	{
	  if (_y_k[i]->size() >= 5) std::cout <<std::endl;
	  _y_k[i]->display();
	}
      else std::cout << "->NULL" <<std::endl;
    }
  for (unsigned int i = 0; i < _upperLevelForInput + 1; i++)
    {
      std::cout << "| lambda[" << i  << "] : ";
      if (_lambda[i])
	{
	  if (_lambda[i]->size() >= 5) std::cout <<std::endl;
	  _lambda[i]->display();
	}
      else std::cout << "->NULL" <<std::endl;
    }


  std::cout << "===================================" <<std::endl;
}

void Interaction::computeOutput(double time, InteractionProperties& interProp, unsigned int derivativeNumber)
{

  DEBUG_BEGIN("Interaction::computeOutput(...)\n");
  DEBUG_PRINTF("time= %f\t",time);
  DEBUG_PRINTF("derivativeNumber= %i\n",derivativeNumber);
  relation()->computeOutput(time, *this, interProp, derivativeNumber);
  DEBUG_END("Interaction::computeOutput(...)\n");

}

void Interaction::computeInput(double time, InteractionProperties& interProp, unsigned int level)
{
  DEBUG_BEGIN("Interaction::computeInput(...)\n");
  DEBUG_PRINTF("time= %f\t",time);
  DEBUG_PRINTF("level= %i\n",level);
  relation()->computeInput(time, *this, interProp, level);
  DEBUG_END("Interaction::computeInput(...)\n");
}




void Interaction::getLeftInteractionBlockForDS(unsigned int pos, SP::SiconosMatrix InteractionBlock, VectorOfSMatrices& workM) const
{
  SP::SiconosMatrix originalMatrix;
  RELATION::TYPES relationType = relation()->getType();
  RELATION::SUBTYPES relationSubType = relation()->getSubType();

  if (relationType == FirstOrder)
  {
    SP::SiconosMatrix CMat = std11::static_pointer_cast<FirstOrderR> (relation())->C();
    if (CMat)
      originalMatrix = CMat;
    else if (relationSubType != LinearTIR)
      originalMatrix = workM[FirstOrderR::mat_C];
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
  subPos[1] = pos;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, InteractionBlock, subDim, subPos);
}

SiconosMatrix& Interaction::getLeftInteractionBlock(VectorOfSMatrices& workM) const
{
  RELATION::TYPES relationType = relation()->getType();
  RELATION::SUBTYPES relationSubType = relation()->getSubType();

  if (relationType == FirstOrder)
  {
    SP::SiconosMatrix CMat = std11::static_pointer_cast<FirstOrderR> (relation())->C();
    if (CMat)
      return *CMat;
    else if (relationSubType != LinearTIR)
      return *workM[FirstOrderR::mat_C];
  }
  else if (relationType == Lagrangian)
  {
    SP::LagrangianR r = std11::static_pointer_cast<LagrangianR> (relation());
    return *r->jachq();
  }
  else if (relationType == NewtonEuler)
  {
    SP::NewtonEulerR r = std11::static_pointer_cast<NewtonEulerR> (relation());
    return *r->jachqT();
  }
  else
  {
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDS, not yet implemented for relations of type " + relationType);
  }
  // stupid compiler check
  return *workM[FirstOrderR::mat_C];

}
void Interaction::getLeftInteractionBlockForDSProjectOnConstraints(unsigned int pos, SP::SiconosMatrix InteractionBlock) const
{
  DEBUG_PRINT("Interaction::getLeftInteractionBlockForDSProjectOnConstraints(unsigned int pos, SP::SiconosMatrix InteractionBlock) \n");
  DEBUG_PRINTF("pos = %i\n", pos);

  if (pos==6)
    pos = pos + 1 ;


  //Type::Siconos dsType = Type::value(*ds);
  //if (dsType != Type::NewtonEulerDS)
  //  RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDSForProject- ds is not from NewtonEulerDS.");

  RELATION::TYPES relationType = relation()->getType();
  if (relationType != NewtonEuler)
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDSForProject- relation is not from NewtonEulerR.");

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
  subPos[1] = pos;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, InteractionBlock, subDim, subPos);
}

void Interaction::getRightInteractionBlockForDS(unsigned int pos, SP::SiconosMatrix InteractionBlock, VectorOfSMatrices& workM) const
{
  SP::SiconosMatrix originalMatrix; // Complete matrix, Relation member.
  RELATION::TYPES relationType = relation()->getType();
  RELATION::SUBTYPES relationSubType = relation()->getSubType();

  if (relationType == FirstOrder)
  {
    SP::SiconosMatrix BMat = std11::static_pointer_cast<FirstOrderR> (relation())->B();
    if (BMat)
      originalMatrix = BMat;
    else if (relationSubType != LinearTIR)
      originalMatrix = workM[FirstOrderR::mat_B];
    else
       RuntimeException::selfThrow("Interaction::getRightInteractionBlockForDS, FirstOrderLinearTIR relation but no B matrix found!");
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
  subPos[0] = pos;
  subPos[1] = 0;//_relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, InteractionBlock, subDim, subPos);
}

void Interaction::getExtraInteractionBlock(SP::SiconosMatrix InteractionBlock, VectorOfSMatrices& workM) const
{
  // !!! Warning: we suppose that D is interactionBlock diagonal, ie that
  // there is no coupling between Interaction through D !!!  Any
  // coupling between relations through D must be taken into account
  // thanks to the nslaw (by "increasing" its dimension).

  RELATION::TYPES relationType = relation()->getType();
  RELATION::SUBTYPES relationSubType = relation()->getSubType();
  SP::SiconosMatrix D;

  if (relationType == FirstOrder)
  {
    SP::SiconosMatrix DMat = std11::static_pointer_cast<FirstOrderR> (relation())->D();
    if (DMat)
      D = DMat;
    else if (relationSubType != LinearTIR)
      D = workM[FirstOrderR::mat_D];
  }
  else if (relationType == Lagrangian)
  {
    D = std11::static_pointer_cast<LagrangianR> (relation())->jachlambda();
  }
  else if (relationType == NewtonEuler)
  {
    D = std11::static_pointer_cast<NewtonEulerR> (relation())->jachlambda();
  }
  else
    RuntimeException::selfThrow("Interaction::getLeftInteractionBlockForDS, not yet implemented for relations of type " + relationType);

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

void Interaction::computeKhat(SiconosMatrix& m, VectorOfSMatrices& workM, double h) const
{
  RELATION::TYPES relationType = relation()->getType();

  if ((relationType == FirstOrder) && (workM[FirstOrderR::mat_Khat]))
  {
    SP::SiconosMatrix K = std11::static_pointer_cast<FirstOrderR>(_relation)->K();
    if (!K) K = workM[FirstOrderR::mat_K];
    prod(*K, m, *workM[FirstOrderR::mat_Khat], true);
    *workM[FirstOrderR::mat_Khat] *= h;
  }
}

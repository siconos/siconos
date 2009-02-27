/* Siconos-Kernel version 3.0.0, Copyright INRIA 2005-2008.
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
#include "LagrangianDS.h"
#include "UnitaryRelation.h"
#include "RelationTypes.hpp"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"
#include "RuntimeException.h"

using namespace std;
using namespace RELATION;

// --- CONSTRUCTORS ---

// Data constructor
UnitaryRelation::UnitaryRelation(SP::Interaction inter, unsigned int pos, unsigned int num): mainInteraction(inter), relativePosition(pos), number(num)
{}

// --- DESTRUCTOR ---
UnitaryRelation::~UnitaryRelation()
{
}


const VectorOfVectors UnitaryRelation::getY() const
{
  // A new object of type VectorOfVectors is created but it handles
  // pointers to BlockVectors, thus there is no copy of the "basic"
  // SimpleVectors.

  VectorOfVectors tmp;
  VectorOfVectors interactionUnitaryBlocks = mainInteraction->getY();

  for (unsigned int i = 0; i < interactionUnitaryBlocks.size(); ++i)
    tmp[i] = interactionUnitaryBlocks[i]->getVectorPtr(number);

  return tmp;
}

SP::SiconosVector UnitaryRelation::getYPtr(unsigned int i) const
{
  // i is the derivative number.
  return ((mainInteraction->getYPtr(i))->getVectorPtr(number));
}

SP::SiconosVector UnitaryRelation::getYOldPtr(unsigned int i) const
{
  // i is the derivative number.
  return ((mainInteraction->getYOldPtr(i))->getVectorPtr(number));
}

const VectorOfVectors UnitaryRelation::getLambda() const
{
  // A new object of type VectorOfVectors is created but it handles
  // pointers to BlockVectors, thus there is no copy of the "basic"
  // SimpleVectors.
  VectorOfVectors tmp;
  VectorOfVectors interactionUnitaryBlocks = mainInteraction->getLambda();

  for (unsigned int i = 0; i < interactionUnitaryBlocks.size(); ++i)
    tmp[i] = interactionUnitaryBlocks[i]->getVectorPtr(number);

  return tmp;
}

SP::SiconosVector UnitaryRelation::getLambdaPtr(unsigned int i) const
{
  // i is the derivative number.
  return ((mainInteraction->getLambdaPtr(i))->getVectorPtr(number));
}

const double UnitaryRelation::getYRef(unsigned int i) const
{
  // get the single value used to build indexSets Warning: the
  // relativePosition depends on NsLawSize and/or type.  This means
  // that at the time, for the unitaryBlock of y that corresponds to
  // the present relation, the first scalar value is used.  For
  // example, for friction, normal part is in first position, followed
  // by the tangential parts.
  return (*getYPtr(i))(0);
}

const double UnitaryRelation::getLambdaRef(unsigned int i) const
{
  // get the single value used to build indexSets
  return (*getLambdaPtr(i))(0);
}

const unsigned int UnitaryRelation::getNonSmoothLawSize() const
{
  return mainInteraction->getNonSmoothLawPtr()->getNsLawSize();
}

const string UnitaryRelation::getNonSmoothLawType() const
{
  return mainInteraction->getNonSmoothLawPtr()->getType();
}

const RELATION::TYPES UnitaryRelation::getRelationType() const
{
  return mainInteraction->getRelationPtr()->getType();
}

const RELATION::SUBTYPES UnitaryRelation::getRelationSubType() const
{
  return mainInteraction->getRelationPtr()->getSubType();
}

SP::DynamicalSystemsSet UnitaryRelation::getDynamicalSystemsPtr()
{
  return mainInteraction->getDynamicalSystemsPtr();
}

void UnitaryRelation::initialize(const std::string& simulationType)
{
  if (!mainInteraction)
    RuntimeException::selfThrow("UnitaryRelation::initialize() failed: the linked interaction is NULL.");

  workX.reset(new BlockVector());
  workZ.reset(new BlockVector());

  for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
    workZ->insertPtr((*it)->getZPtr());

  //   if(simulationType == "TimeStepping")
  //     {
  //     }
  if (simulationType == "EventDriven")
  {
    RELATION::TYPES pbType = getRelationType();
    if (pbType == FirstOrder)
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
        workX->insertPtr((*it)->getXPtr());
    }
    else // Lagrangian
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
        workX->insertPtr((boost::static_pointer_cast<LagrangianDS>(*it))->getVelocityPtr());
    }
  }
  //   else
  //     RuntimeException::selfThrow("UnitaryRelation::initialize(simulationType) failed: unknown simulation type.");
}

void UnitaryRelation::getLeftUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock, unsigned index) const
{
  unsigned int k = 0;
  DSIterator itDS;
  itDS = mainInteraction->dynamicalSystemsBegin();

  // look for ds and its position in G
  while (*itDS != ds && itDS != mainInteraction->dynamicalSystemsEnd())
  {
    k += (*itDS)->getDim();
    itDS++;
  }

  // check dimension (1)
  if ((*itDS)->getDim() != UnitaryBlock->size(1))
    RuntimeException::selfThrow("UnitaryRelation::getLeftUnitaryBlockForDS(DS, UnitaryBlock, ...): inconsistent sizes between UnitaryBlock and DS");

  SP::SiconosMatrix originalMatrix;

  RELATION::TYPES relationType = getRelationType();
  RELATION::SUBTYPES relationSubType = getRelationSubType();

  if (relationType == FirstOrder)
  {
    originalMatrix = mainInteraction->getRelationPtr()->getJacHPtr(0);
  }
  else if (relationType == Lagrangian)
  {
    originalMatrix = mainInteraction->getRelationPtr()->getJacHPtr(index);
  }
  else
    RuntimeException::selfThrow("UnitaryRelation::getLeftUnitaryBlockForDS, not yet implemented for relations of type " + relationType);

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}

void UnitaryRelation::getRightUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock, unsigned index) const
{
  unsigned int k = 0;
  DSIterator itDS;
  itDS = mainInteraction->dynamicalSystemsBegin();

  // look for ds and its position in G
  while (*itDS != ds && itDS != mainInteraction->dynamicalSystemsEnd())
  {
    k += (*itDS)->getDim();
    itDS++;
  }

  // check dimension (1)
  if ((*itDS)->getDim() != UnitaryBlock->size(0))
    RuntimeException::selfThrow("UnitaryRelation::getRightUnitaryBlockForDS(DS, UnitaryBlock, ...): inconsistent sizes between UnitaryBlock and DS");


  SP::SiconosMatrix originalMatrix; // Complete matrix, Relation member.
  RELATION::TYPES relationType = getRelationType();
  RELATION::SUBTYPES relationSubType = getRelationSubType();

  if (relationType == FirstOrder)
  {
    originalMatrix = mainInteraction->getRelationPtr()->getJacGPtr(0);
  }
  else if (relationType == Lagrangian) // Note: the transpose will be done in LCP or MLCP ...
  {
    originalMatrix = mainInteraction->getRelationPtr()->getJacHPtr(index);
  }
  else
    RuntimeException::selfThrow("UnitaryRelation::getRightUnitaryBlockForDS, not yet implemented for relations of type " + relationType);

  if (! originalMatrix)
    RuntimeException::selfThrow("UnitaryRelation::getRightUnitaryBlockForDS(DS, UnitaryBlock, ...): the right unitaryBlock is a NULL pointer (miss matrix B or H or gradients ...in relation ?)");

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix
  // and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = k;
  subPos[1] = relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}

void UnitaryRelation::getExtraUnitaryBlock(SP::SiconosMatrix UnitaryBlock) const
{
  // !!! Warning: we suppose that D is unitaryBlock diagonal, ie that
  // there is no coupling between UnitaryRelation through D !!!  Any
  // coupling between relations through D must be taken into account
  // thanks to the nslaw (by "increasing" its dimension).

  RELATION::TYPES relationType = getRelationType();
  RELATION::SUBTYPES relationSubType = getRelationSubType();

  SP::SiconosMatrix D;
  if (mainInteraction->getRelationPtr()->getNumberOfJacobiansForH() > 1)
    D = mainInteraction->getRelationPtr()->getJacHPtr(1);

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
  subPos[0] = relativePosition;
  subPos[1] = relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(D, UnitaryBlock, subDim, subPos);
}


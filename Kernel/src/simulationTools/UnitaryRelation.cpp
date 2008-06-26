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
#include "RelationTypes.h"
#include "NewtonImpactNSL.h"
#include "NewtonImpactFrictionNSL.h"
#include "RuntimeException.h"

using namespace std;

// --- CONSTRUCTORS ---

// Data constructor
UnitaryRelation::UnitaryRelation(Interaction* inter, unsigned int pos, unsigned int num): mainInteraction(inter), relativePosition(pos), number(num),
  workX(NULL), workZ(NULL)
{}

// --- DESTRUCTOR ---
UnitaryRelation::~UnitaryRelation()
{
  mainInteraction = NULL;
  delete workX;
  delete workZ;
}

const VectorOfVectors UnitaryRelation::getY() const
{
  // A new object of type VectorOfVectors is created but it handles pointers to BlockVectors,
  // thus there is no copy of the "basic" SimpleVectors.

  VectorOfVectors tmp;
  VectorOfVectors interactionUnitaryBlocks = mainInteraction->getY();

  for (unsigned int i = 0; i < interactionUnitaryBlocks.size(); ++i)
    tmp[i] = interactionUnitaryBlocks[i]->getVectorPtr(number);

  return tmp;
}

SiconosVector* UnitaryRelation::getYPtr(unsigned int i) const
{
  // i is the derivative number.
  return ((mainInteraction->getYPtr(i))->getVectorPtr(number));
}

SiconosVector* UnitaryRelation::getYOldPtr(unsigned int i) const
{
  // i is the derivative number.
  return ((mainInteraction->getYOldPtr(i))->getVectorPtr(number));
}

const VectorOfVectors UnitaryRelation::getLambda() const
{
  // A new object of type VectorOfVectors is created but it handles pointers to BlockVectors,
  // thus there is no copy of the "basic" SimpleVectors.
  VectorOfVectors tmp;
  VectorOfVectors interactionUnitaryBlocks = mainInteraction->getLambda();

  for (unsigned int i = 0; i < interactionUnitaryBlocks.size(); ++i)
    tmp[i] = interactionUnitaryBlocks[i]->getVectorPtr(number);

  return tmp;
}

SiconosVector* UnitaryRelation::getLambdaPtr(unsigned int i) const
{
  // i is the derivative number.
  return ((mainInteraction->getLambdaPtr(i))->getVectorPtr(number));
}

const double UnitaryRelation::getYRef(unsigned int i) const
{
  // get the single value used to build indexSets
  // Warning: the relativePosition depends on NsLawSize and/or type.
  // This means that at the time, for the unitaryBlock of y that corresponds to the present relation, the first scalar value is used.
  // For example, for friction, normal part is in first position, followed by the tangential parts.
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

const RELATIONTYPES UnitaryRelation::getRelationType() const
{
  return mainInteraction->getRelationPtr()->getType();
}

const RELATIONSUBTYPES UnitaryRelation::getRelationSubType() const
{
  return mainInteraction->getRelationPtr()->getSubType();
}

DynamicalSystemsSet * UnitaryRelation::getDynamicalSystemsPtr()
{
  return mainInteraction->getDynamicalSystemsPtr();
}

void UnitaryRelation::initialize(const std::string& simulationType)
{
  if (mainInteraction == NULL)
    RuntimeException::selfThrow("UnitaryRelation::initialize() failed: the linked interaction is NULL.");

  workX = new BlockVector();
  workZ = new BlockVector();

  for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
    workZ->insertPtr((*it)->getZPtr());

  //   if(simulationType == "TimeStepping")
  //     {
  //     }
  if (simulationType == "EventDriven")
  {
    RELATIONTYPES pbType = getRelationType();
    if (pbType == FirstOrder)
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
        workX->insertPtr((*it)->getXPtr());
    }
    else // Lagrangian
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
        workX->insertPtr((static_cast<LagrangianDS*>(*it))->getVelocityPtr());
    }
  }
  //   else
  //     RuntimeException::selfThrow("UnitaryRelation::initialize(simulationType) failed: unknown simulation type.");
}

void UnitaryRelation::getLeftUnitaryBlockForDS(DynamicalSystem * ds, SiconosMatrix* UnitaryBlock, unsigned index) const
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

  SiconosMatrix * originalMatrix = NULL; // Complete matrix, Relation member.

  RELATIONTYPES relationType = getRelationType();
  RELATIONSUBTYPES relationSubType = getRelationSubType();

  if (relationType == FirstOrder)
  {
    if (relationSubType == Type1R)//|| relationType =="FirstOrderType2R" || relationType =="FirstOrderType3R")
      originalMatrix = (static_cast<FirstOrderR*>(mainInteraction->getRelationPtr()))->getJacobianHPtr(0);

    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      originalMatrix = (static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr()))->getCPtr();
  }
  else if (relationType == Lagrangian)
  {
    if (relationSubType == ScleronomousR || relationSubType == RheonomousR || relationSubType == CompliantR)
      originalMatrix = (static_cast<LagrangianR*>(mainInteraction->getRelationPtr()))->getGPtr(index);
    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      originalMatrix = (static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr()))->getHPtr();
  }
  else
    RuntimeException::selfThrow("UnitaryRelation::getLeftUnitaryBlockForDS, not yet implemented for relations of type " + relationType);

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}

void UnitaryRelation::getRightUnitaryBlockForDS(DynamicalSystem * ds, SiconosMatrix* UnitaryBlock, unsigned index) const
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

  SiconosMatrix * originalMatrix = NULL; // Complete matrix, Relation member.
  RELATIONTYPES relationType = getRelationType();
  RELATIONSUBTYPES relationSubType = getRelationSubType();

  if (relationType == FirstOrder)
  {
    if (relationSubType == Type1R)//|| relationType =="FirstOrderType2R" || relationType =="FirstOrderType3R")
      originalMatrix = (static_cast<FirstOrderR*>(mainInteraction->getRelationPtr()))->getJacobianGPtr(0);

    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      originalMatrix = (static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr()))->getBPtr();
  }
  else if (relationType == Lagrangian)
  {
    if (relationSubType == ScleronomousR || relationSubType == RheonomousR || relationSubType == CompliantR)
      originalMatrix = (static_cast<LagrangianR*>(mainInteraction->getRelationPtr()))->getGPtr(index);
    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      originalMatrix = (static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr()))->getHPtr();
  }
  else
    RuntimeException::selfThrow("UnitaryRelation::getRightUnitaryBlockForDS, not yet implemented for relations of type " + relationType);

  if (originalMatrix == NULL)
    RuntimeException::selfThrow("UnitaryRelation::getRightUnitaryBlockForDS(DS, UnitaryBlock, ...): the right unitaryBlock is a NULL pointer (miss matrix B or H or gradients ...in relation ?)");

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = k;
  subPos[1] = relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}

void UnitaryRelation::getExtraUnitaryBlock(SiconosMatrix* UnitaryBlock) const
{
  // !!! Warning: we suppose that D is unitaryBlock diagonal, ie that there is no coupling between UnitaryRelation through D !!!
  // Any coupling between relations through D must be taken into account thanks to the nslaw (by "increasing" its dimension).

  RELATIONTYPES relationType = getRelationType();
  RELATIONSUBTYPES relationSubType = getRelationSubType();
  SiconosMatrix * D = NULL;

  if (relationType == FirstOrder)
  {
    if (relationSubType == Type1R)//|| relationType =="FirstOrderType2R" || relationType =="FirstOrderType3R")
    {
      // nothing, D = NULL
    }

    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      D = static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr())->getDPtr();
  }
  else if (relationType == Lagrangian)
  {
    if (relationSubType == ScleronomousR || relationSubType == RheonomousR)
    {
      // nothing, D = NULL
    }
    else if (relationSubType == CompliantR)
      D = static_cast<LagrangianCompliantR*>(mainInteraction->getRelationPtr())->getGPtr(1);
    else if (relationSubType == LinearR || relationSubType == LinearTIR)
      D = static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr())->getDPtr();
  }
  else
    RuntimeException::selfThrow("UnitaryRelation::getExtraUnitaryBlock, not yet implemented for relations of type " + relationType);

  if (D == NULL)
  {
    UnitaryBlock->zero();
    return; //ie no extra unitaryBlock
  }

  // copy sub-unitaryBlock of originalMatrix into UnitaryBlock
  // dim of the sub-unitaryBlock
  Index subDim(2);
  subDim[0] = UnitaryBlock->size(0);
  subDim[1] = UnitaryBlock->size(1);
  // Position (row,col) of first element to be read in originalMatrix and of first element to be set in UnitaryBlock
  Index subPos(4);
  subPos[0] = relativePosition;
  subPos[1] = relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(D, UnitaryBlock, subDim, subPos);
}




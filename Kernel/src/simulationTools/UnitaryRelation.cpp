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
#include "FirstOrderNonLinearDS.h"
#include "FirstOrderR.h"

using namespace std;
using namespace RELATION;

// --- CONSTRUCTORS ---

// Data constructor
UnitaryRelation::UnitaryRelation(SP::Interaction inter,
                                 unsigned int pos,
                                 unsigned int num): _mainInteraction(inter),
  _relativePosition(pos),
  _number(num)
{}

// --- DESTRUCTOR ---
UnitaryRelation::~UnitaryRelation()
{
}

SP::SiconosVector UnitaryRelation::y(unsigned int i) const
{
  // i is the derivative number.
  return (interaction()->y(i)->vector(_number));
}

SP::SiconosVector UnitaryRelation::yOld(unsigned int i) const
{
  // i is the derivative number.
  return (interaction()->yOld(i)->vector(_number));
}

const VectorOfVectors UnitaryRelation::getLambda() const
{
  // A new object of type VectorOfVectors is created but it handles
  // pointers to BlockVectors, thus there is no copy of the "basic"
  // SimpleVectors.
  VectorOfVectors tmp;
  VectorOfVectors interactionUnitaryBlocks = interaction()->getLambda();

  for (unsigned int i = 0; i < interactionUnitaryBlocks.size(); ++i)
    tmp[i] = interactionUnitaryBlocks[i]->vector(_number);

  return tmp;
}

SP::SiconosVector UnitaryRelation::lambda(unsigned int i) const
{
  // i is the derivative number.
  return ((interaction()->lambda(i))->vector(_number));
}

const double UnitaryRelation::getYRef(unsigned int i) const
{
  // get the single value used to build indexSets Warning: the
  // relativePosition depends on NsLawSize and/or type.  This means
  // that at the time, for the unitaryBlock of y that corresponds to
  // the present relation, the first scalar value is used.  For
  // example, for friction, normal part is in first position, followed
  // by the tangential parts.
  return (*y(i))(0);
}

const double UnitaryRelation::getLambdaRef(unsigned int i) const
{
  // get the single value used to build indexSets
  return (*lambda(i))(0);
}

const unsigned int UnitaryRelation::getNonSmoothLawSize() const
{
  return interaction()->nonSmoothLaw()->size();
}

const RELATION::TYPES UnitaryRelation::getRelationType() const
{
  return interaction()->relation()->getType();
}

const RELATION::SUBTYPES UnitaryRelation::getRelationSubType() const
{
  return interaction()->relation()->getSubType();
}

SP::DynamicalSystemsSet UnitaryRelation::dynamicalSystems()
{
  return interaction()->dynamicalSystems();
}

void UnitaryRelation::initialize(const std::string& simulationType)
{
  if (!interaction())
    RuntimeException::selfThrow("UnitaryRelation::initialize() failed: the linked interaction is NULL.");

  _workX.reset(new BlockVector());
  _workZ.reset(new BlockVector());
  mWorkXq.reset(new BlockVector());
  _workFree.reset(new BlockVector());

  for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
    _workZ->insertPtr((*it)->z());

  if (simulationType == "TimeStepping")
  {
    RELATION::TYPES pbType = getRelationType();
    if (pbType == FirstOrder)
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
      {
        SP::FirstOrderNonLinearDS fds = boost::static_pointer_cast<FirstOrderNonLinearDS>(*it);
        _workX->insertPtr(fds->x());
        _workFree->insertPtr(fds->workFree());
        mWorkXq->insertPtr(fds->xq());
      }
    }
  }

  if (simulationType == "EventDriven")
  {
    RELATION::TYPES pbType = getRelationType();
    if (pbType == FirstOrder)
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
      {
        _workX->insertPtr((*it)->x());
        //      mWorkXq->insertPtr((*it)->xq());
      }
    }
    else // Lagrangian
    {
      for (DSIterator it = dynamicalSystemsBegin(); it != dynamicalSystemsEnd(); ++it)
        _workX->insertPtr((boost::static_pointer_cast<LagrangianDS>(*it))->velocity());
    }
  }
  //   else
  //     RuntimeException::selfThrow("UnitaryRelation::initialize(simulationType) failed: unknown simulation type.");
}

void UnitaryRelation::getLeftUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
{
  unsigned int k = 0;
  DSIterator itDS;
  itDS = interaction()->dynamicalSystemsBegin();

  // look for ds and its position in G
  while (*itDS != ds && itDS != interaction()->dynamicalSystemsEnd())
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
    SP::FirstOrderR r = boost::static_pointer_cast<FirstOrderR> (interaction()->relation());
    originalMatrix = r->jacXH();
  }
  else if (relationType == Lagrangian)
  {
    SP::LagrangianR r = boost::static_pointer_cast<LagrangianR> (interaction()->relation());
    originalMatrix = r->jacQH();
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
  subPos[0] = _relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, UnitaryBlock, subDim, subPos);
}

void UnitaryRelation::getRightUnitaryBlockForDS(SP::DynamicalSystem ds, SP::SiconosMatrix UnitaryBlock) const
{
  unsigned int k = 0;
  DSIterator itDS;
  itDS = interaction()->dynamicalSystemsBegin();

  // look for ds and its position in G
  while (*itDS != ds && itDS != interaction()->dynamicalSystemsEnd())
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
    originalMatrix = interaction()->relation()->jacLG();
  }
  else if (relationType == Lagrangian) // Note: the transpose will be done in LCP or MLCP ...
  {
    originalMatrix = interaction()->relation()->jacLH();
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
  subPos[1] = _relativePosition;
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
  //  if(interaction()->relation()->getNumberOfJacobiansForH()>1)
  D = interaction()->relation()->jacLH();

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
  subPos[0] = _relativePosition;
  subPos[1] = _relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(D, UnitaryBlock, subDim, subPos);
}


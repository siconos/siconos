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
  VectorOfVectors interactionBlocks = mainInteraction->getY();

  for (unsigned int i = 0; i < interactionBlocks.size(); ++i)
    tmp[i] = interactionBlocks[i]->getVectorPtr(number);

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
  VectorOfVectors interactionBlocks = mainInteraction->getLambda();

  for (unsigned int i = 0; i < interactionBlocks.size(); ++i)
    tmp[i] = interactionBlocks[i]->getVectorPtr(number);

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
  // This means that at the time, for the block of y that corresponds to the present relation, the first scalar value is used.
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

const string UnitaryRelation::getRelationType() const
{
  return mainInteraction->getRelationPtr()->getType();
}

const string UnitaryRelation::getRelationSubType() const
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
    string pbType = getRelationType();
    if (pbType == "FirstOrder")
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

void UnitaryRelation::getLeftBlockForDS(DynamicalSystem * ds, SiconosMatrix* Block, unsigned index) const
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
  if ((*itDS)->getDim() != Block->size(1))
    RuntimeException::selfThrow("UnitaryRelation::getLeftBlockForDS(DS, Block, ...): inconsistent sizes between Block and DS");

  SiconosMatrix * originalMatrix = NULL; // Complete matrix, Relation member.

  string relationType = getRelationType() + getRelationSubType();

  if (relationType == "FirstOrderType1R" || relationType == "FirstOrderType2R" || relationType == "FirstOrderType3R")
    originalMatrix = (static_cast<FirstOrderR*>(mainInteraction->getRelationPtr()))->getJacobianHPtr(0);

  else if (relationType == "FirstOrderLinearR" || relationType == "FirstOrderLinearTIR")
    originalMatrix = (static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr()))->getCPtr();

  else if (relationType == "LagrangianScleronomousR" || relationType == "LagrangianRheonomousR" || relationType == "LagrangianCompliantR")
    originalMatrix = (static_cast<LagrangianR*>(mainInteraction->getRelationPtr()))->getGPtr(index);

  else if (relationType == "LagrangianLinearR")
    originalMatrix = (static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr()))->getHPtr();

  else
    RuntimeException::selfThrow("UnitaryRelation::getLeftBlockForDS, not yet implemented for relations of type " + relationType);

  // copy sub-block of originalMatrix into Block
  // dim of the sub-block
  Index subDim(2);
  subDim[0] = Block->size(0);
  subDim[1] = Block->size(1);
  // Position (row,col) of first element to be read in originalMatrix and of first element to be set in Block
  Index subPos(4);
  subPos[0] = relativePosition;
  subPos[1] = k;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, Block, subDim, subPos);
}

void UnitaryRelation::getRightBlockForDS(DynamicalSystem * ds, SiconosMatrix* Block, unsigned index) const
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
  if ((*itDS)->getDim() != Block->size(0))
    RuntimeException::selfThrow("UnitaryRelation::getRightBlockForDS(DS, Block, ...): inconsistent sizes between Block and DS");

  SiconosMatrix * originalMatrix = NULL; // Complete matrix, Relation member.
  string relationType = getRelationType() + getRelationSubType();

  if (relationType == "FirstOrderType1R" || relationType == "FirstOrderType2R" || relationType == "FirstOrderType3R")
    originalMatrix = (static_cast<FirstOrderR*>(mainInteraction->getRelationPtr()))->getJacobianGPtr(0);

  else if (relationType == "FirstOrderLinearR" || relationType == "FirstOrderLinearTIR")
    originalMatrix = (static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr()))->getBPtr();

  else if (relationType == "LagrangianScleronomousR" || relationType == "LagrangianRheonomousR" || relationType == "LagrangianCompliantR")
    originalMatrix = (static_cast<LagrangianR*>(mainInteraction->getRelationPtr()))->getGPtr(index);

  else if (relationType == "LagrangianLinearR")
    originalMatrix = (static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr()))->getHPtr();

  else RuntimeException::selfThrow("UnitaryRelation::getRightBlockForDS, not yet implemented for relation of type " + relationType);

  if (originalMatrix == NULL)
    RuntimeException::selfThrow("UnitaryRelation::getRightBlockForDS(DS, Block, ...): the right block is a NULL pointer (miss matrix B or H or gradients ...in relation ?)");

  // copy sub-block of originalMatrix into Block
  // dim of the sub-block
  Index subDim(2);
  subDim[0] = Block->size(0);
  subDim[1] = Block->size(1);
  // Position (row,col) of first element to be read in originalMatrix and of first element to be set in Block
  Index subPos(4);
  subPos[0] = k;
  subPos[1] = relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(originalMatrix, Block, subDim, subPos);
}

void UnitaryRelation::getExtraBlock(SiconosMatrix* Block) const
{
  // !!! Warning: we suppose that D is block diagonal, ie that there is no coupling between UnitaryRelation through D !!!
  // Any coupling between relations through D must be taken into account thanks to the nslaw (by "increasing" its dimension).

  string relationType = getRelationType() + getRelationSubType();

  SiconosMatrix * D = NULL;
  if (relationType == "FirstOrderType1R" || relationType == "LagrangianScleronomousR" || relationType == "LagrangianRheonomousR")
  {
    // nothing, D = NULL
  }
  else if (relationType == "FirstOrderType2R" || relationType == "FirstOrderType3R")
    D = static_cast<FirstOrderR*>(mainInteraction->getRelationPtr())->getJacobianHPtr(1);
  else if (relationType == "FirstOrderLinearTIR")
    D = static_cast<FirstOrderLinearTIR*>(mainInteraction->getRelationPtr())->getDPtr();
  else if (relationType == "FirstOrderLinearR")
    D = static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr())->getDPtr();
  else if (relationType == "LagrangianCompliantR")
    D = static_cast<LagrangianCompliantR*>(mainInteraction->getRelationPtr())->getGPtr(1);
  else if (relationType == "LagrangianLinearR")
    D = static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr())->getDPtr();
  else
    RuntimeException::selfThrow("UnitaryRelation::getExtraBlock, not yet implemented for first order relations of subtype " + relationType);

  if (D == NULL)
  {
    Block->zero();
    return; //ie no extra block
  }

  // copy sub-block of originalMatrix into Block
  // dim of the sub-block
  Index subDim(2);
  subDim[0] = Block->size(0);
  subDim[1] = Block->size(1);
  // Position (row,col) of first element to be read in originalMatrix and of first element to be set in Block
  Index subPos(4);
  subPos[0] = relativePosition;
  subPos[1] = relativePosition;
  subPos[2] = 0;
  subPos[3] = 0;
  setBlock(D, Block, subDim, subPos);
}

void UnitaryRelation::computeEquivalentY(double time, unsigned int level, const string& simulationType, SiconosVector* yOut, unsigned int pos)
{

  // Get relation and non smooth law types
  string relationType = getRelationType() + getRelationSubType();
  string nslawType = getNonSmoothLawType();
  // Warning: first version with OneStepNSProblem as an input argument. But this means "inclusion" of simulationTools class into a
  // modelingTools class => no!
  //   string simulationType = osns->getSimulationPtr()->getType();
  // unsigned int level = osns->getLevelMin(); // this corresponds to the derivative order (for y) used to compute yOut.



  //   if(simulationType == "TimeStepping")
  //     {
  unsigned int sizeY = getNonSmoothLawSize();
  std::vector<unsigned int> coord(8);
  coord[0] = relativePosition;
  coord[1] = relativePosition + sizeY;
  coord[2] = 0;
  coord[4] = 0;
  coord[6] = pos;
  coord[7] = pos + sizeY;

  //mainInteraction->computeFreeOutput(time,level);
  SiconosMatrix * H = NULL;
  if (relationType == "FirstOrderType1R" || relationType == "FirstOrderType2R" || relationType == "FirstOrderType3R")
  {
    H = static_cast<FirstOrderR*>(mainInteraction->getRelationPtr())->getJacobianHPtr(0);
    if (H != NULL)
    {
      coord[3] = H->size(1);
      coord[5] = H->size(1);
      subprod(*H, *workX, *yOut, coord, true);
    }
  }

  else if (relationType == "FirstOrderLinearTIR" || relationType == "FirstOrderLinearR")
  {
    // yOut = HXfree + e + Fz
    H = static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr())->getCPtr();
    if (H != NULL)
    {
      coord[3] = H->size(1);
      coord[5] = H->size(1);
      subprod(*H, *workX, (*yOut), coord, true);
    }
    SiconosVector * e = static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr())->getEPtr();
    if (e != NULL)
      static_cast<SimpleVector*>(yOut)->addBlock(pos, *e);

    H = static_cast<FirstOrderLinearR*>(mainInteraction->getRelationPtr())->getFPtr();
    if (H != NULL)
    {
      coord[3] = H->size(1);
      coord[5] = H->size(1);
      subprod(*H, *workZ, *yOut, coord, false);
    }
  }
  else if (relationType == "LagrangianCompliantR" || relationType == "LagrangianScleronomousR" || relationType == "LagrangianRheonomousR")
  {
    // yOut = jacobian_q h().v_free
    H = static_cast<LagrangianR*>(mainInteraction->getRelationPtr())->getGPtr(0);
    if (H != NULL)
    {
      coord[3] = H->size(1);
      coord[5] = H->size(1);
      subprod(*H, *workX, *yOut, coord, true);
    }
  }

  else if (relationType == "LagrangianLinearR")
  {
    // yOut = H.v_free
    H = static_cast<LagrangianLinearR*>(mainInteraction->getRelationPtr())->getHPtr();
    if (H != NULL)
    {
      coord[3] = H->size(1);
      coord[5] = H->size(1);
      subprod(*H, *workX, *yOut, coord, true);
    }
  }
  else
    RuntimeException::selfThrow("UnitaryRelation::getExtraBlock, not yet implemented for first order relations of subtype " + relationType);

  //     }
  //   else  if(simulationType == "EventDriven")
  //     {
  //       //       mainInteraction->computeOutput(time,level);
  //       //       (*yOut) = *(getYPtr(level));
  //     }
  // Add "non-smooth law effect" on yOut
  if (getRelationType() == "Lagrangian")
  {
    double e;
    if (nslawType == NEWTONIMPACTNSLAW)
    {
      e = (static_cast<NewtonImpactNSL*>(mainInteraction->getNonSmoothLawPtr()))->getE();
      std::vector<unsigned int> subCoord(4);
      if (simulationType == "TimeStepping")
      {
        subCoord[0] = 0;
        subCoord[1] = getNonSmoothLawSize();
        subCoord[2] = pos;
        subCoord[3] = pos + subCoord[1];
        subscal(e, *getYOldPtr(level), *yOut, subCoord, false);
      }
      else if (simulationType == "EventDriven")
      {
        subCoord[0] = pos;
        subCoord[1] = pos + getNonSmoothLawSize();
        subCoord[2] = pos;
        subCoord[3] = subCoord[1];
        subscal(e, *yOut, *yOut, subCoord, false); // yOut = yOut + e * yOut
      }
      else
        RuntimeException::selfThrow("UnitaryRelation::computeEquivalentY not yet implemented for relation of type " + relationType + " and non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);
    }
    else if (nslawType == NEWTONIMPACTFRICTIONNSLAW)
    {
      e = (static_cast<NewtonImpactFrictionNSL*>(mainInteraction->getNonSmoothLawPtr()))->getEn();
      // Only the normal part is multiplied by e
      if (simulationType == "TimeStepping")
        (*yOut)(pos) +=  e * (*getYOldPtr(level))(0);

      else RuntimeException::selfThrow("UnitaryRelation::computeEquivalentY not yet implemented for relation of type " + relationType + " and non smooth law of type " + nslawType + " for a simulaton of type " + simulationType);

    }
    else
      RuntimeException::selfThrow("UnitaryRelation::computeEquivalentY not yet implemented for relation of type " + relationType + " and non smooth law of type " + nslawType);
  }
}



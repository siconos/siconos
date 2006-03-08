/* Siconos-Kernel version 1.1.2, Copyright INRIA 2005-2006.
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
#include "FrictionContact.h"

// includes to be deleted thanks to factories
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactFrictionNSL.h"
#include "LinearTIR.h"

using namespace std;

// Default constructor
FrictionContact::FrictionContact(): OneStepNSProblem(), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{}

// xml constructor
FrictionContact::FrictionContact(OneStepNSProblemXML* osNsPbXml, Strategy* newStrat):
  OneStepNSProblem(osNsPbXml, newStrat), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  if (osNsPbXml != NULL)
  {
    FrictionContactXML * xmllcp = (static_cast<FrictionContactXML*>(osNsPbXml));

    // If both q and M are given in xml file, check if sizes are consistent
    if (xmllcp->hasQ() && xmllcp->hasM() && ((xmllcp->getM()).size(0) != (xmllcp->getQ()).size()))
      RuntimeException::selfThrow("FrictionContact: xml constructor, inconsistent sizes between given q and M");

    // first nlcp is given by M matrix size in xml file
    if (xmllcp->hasM())
    {
      dim = (xmllcp->getM()).size(0);
      M = new SiconosMatrix(xmllcp->getM());
      isMAllocatedIn = true;
    }

    if (xmllcp->hasQ())
    {
      // get dim if necessary
      if (M == NULL)
        dim = (xmllcp->getQ()).size();

      q = new SimpleVector(xmllcp->getQ());
      isQAllocatedIn = true;
    }

  }
  else RuntimeException::selfThrow("FrictionContact: xml constructor, xml file=NULL");
}

// Constructor from a set of data
FrictionContact::FrictionContact(Strategy* newStrat):
  OneStepNSProblem(newStrat), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{}

// Constructor from a set of data (2)
FrictionContact::FrictionContact(Strategy* newStrat, Solver*  newSolver):
  OneStepNSProblem(newStrat, newSolver), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{}

// destructor
FrictionContact::~FrictionContact()
{
  if (isWAllocatedIn)
  {
    delete w;
    w = NULL;
  }
  if (isZAllocatedIn)
  {
    delete z;
    z = NULL;
  }
  if (isMAllocatedIn)
  {
    delete M;
    M = NULL;
  }
  if (isQAllocatedIn)
  {
    delete q;
    q = NULL;
  }
}

// Setters

void FrictionContact::setW(const SimpleVector& newValue)
{
  if (dim != newValue.size())
    RuntimeException::selfThrow("FrictionContact: setW, inconsistent size between given w size and problem size. You should set dim before");

  if (w == NULL)
  {
    w = new SimpleVector(dim);
    isWAllocatedIn = true;
  }
  else if (w->size() != dim)
    RuntimeException::selfThrow("FrictionContact: setW, w size differs from dim");

  *w = newValue;
}

void FrictionContact::setWPtr(SimpleVector* newPtr)
{
  if (dim != newPtr->size())
    RuntimeException::selfThrow("FrictionContact: setWPtr, inconsistent size between given w size and problem size. You should set dim before");

  if (isWAllocatedIn) delete w;
  w = newPtr;
  isWAllocatedIn = false;
}


void FrictionContact::setZ(const SimpleVector& newValue)
{
  if (dim != newValue.size())
    RuntimeException::selfThrow("FrictionContact: setZ, inconsistent size between given z size and problem size. You should set dim before");

  if (z == NULL)
  {
    z = new SimpleVector(dim);
    isZAllocatedIn = true;
  }

  *z = newValue;
}

void FrictionContact::setZPtr(SimpleVector* newPtr)
{
  if (dim != newPtr->size())
    RuntimeException::selfThrow("FrictionContact: setZPtr, inconsistent size between given z size and problem size. You should set dim before");

  if (isZAllocatedIn) delete z;
  z = newPtr;
  isZAllocatedIn = false;
}

void FrictionContact::setM(const SiconosMatrix& newValue)
{
  if (dim != newValue.size(0) || dim != newValue.size(1))
    RuntimeException::selfThrow("FrictionContact: setM, inconsistent size between given M size and problem size. You should set dim before");

  if (M == NULL)
  {
    M = new SiconosMatrix(dim, dim);
    isMAllocatedIn = true;
  }

  *M = newValue;
}

void FrictionContact::setMPtr(SiconosMatrix* newPtr)
{
  if (dim != newPtr->size(0) || dim != newPtr->size(1))
    RuntimeException::selfThrow("FrictionContact: setMPtr, inconsistent size between given M size and problem size. You should set dim before");

  if (isMAllocatedIn) delete M;
  M = newPtr;
  isMAllocatedIn = false;
}


void FrictionContact::setQ(const SimpleVector& newValue)
{
  if (dim != newValue.size())
    RuntimeException::selfThrow("FrictionContact: setQ, inconsistent size between given q size and problem size. You should set dim before");

  if (q == NULL)
  {
    q = new SimpleVector(dim);
    isQAllocatedIn = true;
  }

  *q = newValue;
}

void FrictionContact::setQPtr(SimpleVector* newPtr)
{
  if (dim != newPtr->size())
    RuntimeException::selfThrow("FrictionContact: setQPtr, inconsistent size between given q size and problem size. You should set dim before");

  if (isQAllocatedIn) delete q;
  q = newPtr;
  isQAllocatedIn = false;
}

void FrictionContact::initialize()
{

  // This function performs all steps that are not time dependant
  // General initialize for OneStepNSProblem
  OneStepNSProblem::initialize();

  // compute all block-matrices and save them in OneStepNSProblem
  computeAllBlocks();

  // get topology
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // if all relative degrees are equal to 0 or 1
  if (topology->isTimeInvariant())
  {
    dim = topology->getEffectiveSizeOutput();
    assembleM();
  }
}

void FrictionContact::computeAllBlocks()
{
  // --- get topology ---
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // --- get linked-Interactions map ---
  map< Interaction*, vector<InteractionLink*> > linkedInteractionMap = topology->getLinkedInteractionMap();
  map< Interaction*, vector<InteractionLink*> >::const_iterator itMapInter;

  // --- get interactions list ---
  vector<Interaction*> listInteractions = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
  vector<Interaction*>::iterator iter;

  // --- get time step ---
  double h = strategy->getTimeDiscretisationPtr()->getH(); // time step

  string relationType, dsType;
  unsigned int globalDSSize, sizeInteraction, linkedInteractionSize;

  unsigned int sizeBlock, sizeLinkedBlock; // For each interaction, this size depends on relative degree(r), indexMin and the number of relations (j).
  // sizeBlock = (r - indexMin)*j  ou = j si r=0.
  vector<unsigned int> r; // relative degrees
  unsigned int rMax;
  vector<unsigned int> indexMin;
  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // the current interaction and its size
    Interaction *currentInteraction = *iter;
    sizeInteraction = currentInteraction->getNInteraction();
    // relative degrees
    r = topology->getRelativeDegrees(currentInteraction);
    rMax = r[0]; // !!! we suppose the interaction is homogeneous !!!
    if (rMax == 0) rMax = 1 ; // warning: review r=0 case
    indexMin = topology->getIndexMin(currentInteraction);
    sizeBlock = (rMax - indexMin[0]) * sizeInteraction;

    // --- DIAGONAL BLOCKS MANAGEMENT ---

    // matrix that corresponds to diagonal block for the current interaction
    diagonalBlocksMap[ currentInteraction ] = new SiconosMatrix(sizeBlock, sizeBlock);
    SiconosMatrix * currentMatrixBlock = diagonalBlocksMap[ currentInteraction ];

    // get DS list of the current interaction
    vector<DynamicalSystem*> vDS = currentInteraction ->getDynamicalSystems();
    unsigned int nbDS = vDS.size(); // number of DS

    // sum of all DS sizes
    globalDSSize = 0;
    vector<DynamicalSystem*>::iterator itDS;
    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
    {
      dsType = (*itDS)->getType();
      if (dsType == LNLDS || dsType == LTIDS)
        globalDSSize += (*itDS)->getN() / 2;
      else
        globalDSSize += (*itDS)->getN();
    }

    // Get Wi matrix of each DS concerned by the interaction and assemble global matrix W
    OneStepIntegrator * Osi;
    map<DynamicalSystem*, SiconosMatrix*> W;

    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
    {
      Osi = strategy->getIntegratorOfDSPtr(*itDS); // get OneStepIntegrator of current dynamical system
      if (Osi->getType() == MOREAU_INTEGRATOR)
        W[*itDS] = (static_cast<Moreau*>(Osi))->getWPtr();  // get its W matrix
      else
        RuntimeException::selfThrow("FrictionContact::computeAllBlocks not yet implemented for Integrator of type " + Osi->getType());
    }

    // get the relation of the current interaction and its type
    Relation * RCurrent = currentInteraction -> getRelationPtr();
    relationType = RCurrent->getType();

    if (relationType == LAGRANGIANLINEARRELATION || relationType == LAGRANGIANRELATION)
      computeDiagonalBlocksLagrangianR(RCurrent, sizeInteraction, vDS, W, h, currentMatrixBlock);
    else RuntimeException::selfThrow("FrictionContact::computeAllBlocks not yet implemented for relation of type " + relationType);

    // --- EXTRA-DIAGONAL BLOCKS MANAGEMENT ---

    // check if there are linked interactions with current one, and if so get them.
    itMapInter = linkedInteractionMap.find(currentInteraction);
    if (itMapInter != linkedInteractionMap.end())
    {
      vector<InteractionLink*> ILV = linkedInteractionMap[currentInteraction];
      vector<InteractionLink*>::iterator itLink;

      //map<Interaction*, SiconosMatrix*> tmpMap;
      SiconosMatrix * coupledInteractionsBlock, *coupledInteractionsBlockSym;

      // loop over LinkedInteractions
      for (itLink = ILV.begin(); itLink != ILV.end(); itLink++)
      {
        Interaction * linkedInteraction = (*itLink)->getLinkedInteractionPtr();
        linkedInteractionSize = linkedInteraction->getNInteraction();
        // relative degrees
        r = topology->getRelativeDegrees(linkedInteraction);
        rMax = r[0]; // !!! we suppose the interaction is homogeneous !!!
        if (rMax == 0) rMax = 1 ; // warning: review r=0 case
        indexMin = topology->getIndexMin(linkedInteraction);
        sizeLinkedBlock = (rMax - indexMin[0]) * linkedInteractionSize;

        (extraDiagonalBlocksMap[currentInteraction])[linkedInteraction] = new SiconosMatrix(sizeBlock, sizeLinkedBlock);
        (extraDiagonalBlocksMap[linkedInteraction])[currentInteraction] = new SiconosMatrix(sizeLinkedBlock, sizeBlock);
        coupledInteractionsBlock = extraDiagonalBlocksMap[currentInteraction][linkedInteraction];
        coupledInteractionsBlockSym = extraDiagonalBlocksMap[linkedInteraction][currentInteraction];

        // get the list of common DS and their number
        vector<DynamicalSystem*> commonDS = (*itLink)->getCommonDS();
        nbDS = commonDS.size();

        // get the relation of the current interaction
        Relation * RLinked = linkedInteraction -> getRelationPtr();
        string RlinkedType = RLinked->getType();

        if ((relationType == LAGRANGIANRELATION || relationType == LAGRANGIANLINEARRELATION) && (RlinkedType == LAGRANGIANRELATION ||  RlinkedType == LAGRANGIANLINEARRELATION))
        {
          computeExtraDiagonalBlocksLagrangianR(RCurrent, RLinked, sizeInteraction, linkedInteractionSize, commonDS, W, h, coupledInteractionsBlock);
          computeExtraDiagonalBlocksLagrangianR(RLinked, RCurrent, linkedInteractionSize, sizeInteraction, commonDS, W, h, coupledInteractionsBlockSym);
        }
        else RuntimeException::selfThrow("FrictionContact::computeAllBlocks not yet implemented for relation of type " + relationType);
      } // end of loop over linked interactions
    }
  } // end of loop over interactions -> increment current interaction

}

void FrictionContact::computeDiagonalBlocksLagrangianR(Relation * R, const unsigned int& sizeInteraction, vector<DynamicalSystem*> vDS,
    map<DynamicalSystem*, SiconosMatrix*> W, const double& h, SiconosMatrix* currentMatrixBlock)
{

  // Warning: we suppose that at this point, G and/or h have been computed through plug-in mechanism.

  // convert to lagrangianR
  LagrangianR *LR = static_cast<LagrangianR*>(R);

  // For 1 relation, we need:
  // - a block-row of G, that corresponds to a specific DS
  // - a block of tG, that corresponds to a specific DS
  SiconosMatrix *G;
  vector<DynamicalSystem*>::iterator itDS;
  unsigned int sizeDS;

  bool isHomogeneous = false; // \todo to be defined with relative degrees
  // isHomogeneous = true means that all relative degrees of the interaction are equal.
  if (!isHomogeneous)
  {
    // the resulting row-block
    SimpleVector * currentLine = new SimpleVector(sizeInteraction);
    // a row of G corresponding to a specific DS (its size depends on the DS size)
    SimpleVector * Grow ;

    // compute currentLine = sum(j) Grow,j Wj  tGj, with j the list of DS of the
    // current interaction.
    // loop over the relations of the current interaction
    for (unsigned int i = 0; i < sizeInteraction; i++)
    {
      currentLine->zero();
      // loop over the DS of the current interaction
      for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
      {
        sizeDS = (*itDS)->getN() / 2; // divided by 2 to get nDof
        // get blocks corresponding to the current DS
        G = new SiconosMatrix(sizeInteraction, sizeDS);
        LR->getGBlockDS(*itDS, *G);
        // get row i of G
        Grow = new SimpleVector(sizeDS);
        G->getRow(i, *Grow);

        // compute currentLine
        *currentLine +=  *Grow * (W[*itDS])->multTranspose(*G);
        delete G;
        delete Grow;
      }
      // save the result in currentMatrixBlock (and so in diagonalBlocksMap)
      currentMatrixBlock->setRow(i, *currentLine);
    }
    delete currentLine;
  }
  else
  {
    // the resulting block
    currentMatrixBlock->zero();
    // loop over the DS of the current interaction
    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
    {

      sizeDS = (*itDS)->getN() / 2; // divided by 2 to get nDof
      // get blocks corresponding to the current DS
      G = new SiconosMatrix(sizeInteraction, sizeDS);
      LR->getGBlockDS(*itDS, *G);
      *currentMatrixBlock +=  *G * (W[*itDS])->multTranspose(*G);
      delete G;
    }
  }
}

void FrictionContact::computeExtraDiagonalBlocksLagrangianR(Relation * RCurrent, Relation* RLinked, const unsigned int& sizeInteraction,
    const unsigned int& linkedInteractionSize, vector<DynamicalSystem*> commonDS,
    map<DynamicalSystem*, SiconosMatrix*> W, const double& h, SiconosMatrix* coupledInteractionsBlock)
{
  // Warning: we suppose that at this point, G and/or h have been computed through plug-in mechanism.

  // convert to LagrangianLinear relation
  LagrangianR *LR1 = static_cast<LagrangianR*>(RCurrent);
  LagrangianR *LR2 = static_cast<LagrangianR*>(RLinked);
  SiconosMatrix *Gcurrent, *Glinked;
  coupledInteractionsBlock->zero();
  vector<DynamicalSystem*>::iterator itDS;
  unsigned int sizeDS;
  bool isHomogeneous = false; // \todo to be defined with relative degrees

  if (!isHomogeneous)
  {
    // the resulting row-block
    SimpleVector * currentLine = new SimpleVector(linkedInteractionSize);
    // a row of G corresponding to a specific DS (its size depends on the DS size)
    SimpleVector * Grow ;
    currentLine->zero();
    // compute currentLine =  sum(j) Grow,j Wj  tGj, with j the list of common DS
    // Grow belongs to current Interaction and tG to linked interaction
    // loop over the relations of the current interaction
    for (unsigned int i = 0; i < sizeInteraction; i++)
    {
      // loop over common DS
      for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
      {
        sizeDS = (*itDS)->getN() / 2; // divided by 2 to get nDof
        // get blocks corresponding to the current DS
        Gcurrent = new SiconosMatrix(sizeInteraction, sizeDS);
        Glinked = new SiconosMatrix(linkedInteractionSize, sizeDS);
        LR1->getGBlockDS(*itDS, *Gcurrent);
        LR2->getGBlockDS(*itDS, *Glinked);
        // get row i of G
        Grow = new SimpleVector(sizeDS);
        Gcurrent->getRow(i, *Grow);
        // compute currentLine
        *currentLine +=  *Grow * (W[*itDS])->multTranspose(*Glinked);
        delete Grow;
        delete Gcurrent;
        delete Glinked;
      }
      // save the result in currentMatrixBlock (and so in diagonalBlocksMap)
      coupledInteractionsBlock->setRow(i, *currentLine);
    }
    delete currentLine;
  }
  else
  {
    // the resulting block
    // loop over the DS of the current interaction
    for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
    {
      sizeDS = (*itDS)->getN() / 2; // divided by 2 to get nDof
      // get blocks corresponding to the current DS
      Gcurrent = new SiconosMatrix(sizeInteraction, sizeDS);
      Glinked = new SiconosMatrix(linkedInteractionSize, sizeDS);
      LR1->getGBlockDS(*itDS, *Gcurrent);
      LR2->getGBlockDS(*itDS, *Glinked);
      *coupledInteractionsBlock += *Gcurrent * (W[*itDS])->multTranspose(*Glinked);
      delete Gcurrent;
      delete Glinked;
    }
  }
}

void FrictionContact::updateBlocks()
{
  // --- get topology ---
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // --- get linked-Interactions map ---
  map< Interaction*, vector<InteractionLink*> > linkedInteractionMap = topology->getLinkedInteractionMap();
  map< Interaction*, vector<InteractionLink*> >::const_iterator itMapInter;

  // --- get interactions list ---
  vector<Interaction*> listInteractions = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
  vector<Interaction*>::iterator iter;

  // --- get time step ---
  double h = strategy->getTimeDiscretisationPtr()->getH(); // time step

  string relationType, dsType;
  unsigned int globalDSSize, sizeInteraction, linkedInteractionSize;

  unsigned int sizeBlock, sizeLinkedBlock; // For each interaction, this size depends on relative degree(r), indexMin and the number of relations (j).
  // sizeBlock = (r - indexMin)*j  ou = j si r=0.
  vector<unsigned int> r; // relative degrees
  unsigned int rMax;
  vector<unsigned int> indexMin;
  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // the current interaction and its size
    Interaction *currentInteraction = *iter;
    sizeInteraction = currentInteraction->getNInteraction();
    // relative degrees
    r = topology->getRelativeDegrees(currentInteraction);
    rMax = r[0]; // !!! we suppose the interaction is homogeneous !!!
    if (rMax == 0) rMax = 1 ; // warning: review r=0 case
    indexMin = topology->getIndexMin(currentInteraction);
    sizeBlock = (rMax - indexMin[0]) * sizeInteraction;

    // --- DIAGONAL BLOCKS MANAGEMENT ---

    // matrix that corresponds to diagonal block for the current interaction
    SiconosMatrix * currentMatrixBlock = diagonalBlocksMap[ currentInteraction ];

    // get DS list of the current interaction
    vector<DynamicalSystem*> vDS = currentInteraction ->getDynamicalSystems();
    unsigned int nbDS = vDS.size(); // number of DS

    // sum of all DS sizes
    globalDSSize = 0;
    vector<DynamicalSystem*>::iterator itDS;
    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
    {
      dsType = (*itDS)->getType();
      if (dsType == LNLDS || dsType == LTIDS)
        globalDSSize += (*itDS)->getN() / 2;
      else
        globalDSSize += (*itDS)->getN();
    }

    // Get Wi matrix of each DS concerned by the interaction and assemble global matrix W
    OneStepIntegrator * Osi;
    map<DynamicalSystem*, SiconosMatrix*> W;

    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
    {
      Osi = strategy->getIntegratorOfDSPtr(*itDS); // get OneStepIntegrator of current dynamical system
      if (Osi->getType() == MOREAU_INTEGRATOR)
        W[*itDS] = (static_cast<Moreau*>(Osi))->getWPtr();  // get its W matrix
      else
        RuntimeException::selfThrow("FrictionContact::computeAllBlocks not yet implemented for Integrator of type " + Osi->getType());
    }

    // get the relation of the current interaction and its type
    Relation * RCurrent = currentInteraction -> getRelationPtr();
    relationType = RCurrent->getType();

    if (relationType == LAGRANGIANLINEARRELATION)
    {}// Nothing to be done, blocks allready computed
    else if (relationType == LAGRANGIANRELATION)
      computeDiagonalBlocksLagrangianR(RCurrent, sizeInteraction, vDS, W, h, currentMatrixBlock);
    else RuntimeException::selfThrow("FrictionContact::updateBlocks not yet implemented for relation of type " + relationType);

    // --- EXTRA-DIAGONAL BLOCKS MANAGEMENT ---

    // check if there are linked interactions with current one, and if so get them.
    itMapInter = linkedInteractionMap.find(currentInteraction);
    if (itMapInter != linkedInteractionMap.end())
    {
      vector<InteractionLink*> ILV = linkedInteractionMap[currentInteraction];
      vector<InteractionLink*>::iterator itLink;

      //map<Interaction*, SiconosMatrix*> tmpMap;
      SiconosMatrix * coupledInteractionsBlock, *coupledInteractionsBlockSym;

      // loop over LinkedInteractions
      for (itLink = ILV.begin(); itLink != ILV.end(); itLink++)
      {
        Interaction * linkedInteraction = (*itLink)->getLinkedInteractionPtr();
        linkedInteractionSize = linkedInteraction->getNInteraction();
        // relative degrees
        r = topology->getRelativeDegrees(linkedInteraction);
        rMax = r[0]; // !!! we suppose the interaction is homogeneous !!!
        if (rMax == 0) rMax = 1 ; // warning: review r=0 case
        indexMin = topology->getIndexMin(linkedInteraction);
        sizeLinkedBlock = (rMax - indexMin[0]) * linkedInteractionSize;

        coupledInteractionsBlock = extraDiagonalBlocksMap[currentInteraction][linkedInteraction];
        coupledInteractionsBlockSym = extraDiagonalBlocksMap[linkedInteraction][currentInteraction];

        // get the list of common DS and their number
        vector<DynamicalSystem*> commonDS = (*itLink)->getCommonDS();
        nbDS = commonDS.size();

        // get the relation of the current interaction
        Relation * RLinked = linkedInteraction -> getRelationPtr();
        string RlinkedType = RLinked->getType();
        if (relationType == LAGRANGIANLINEARRELATION && RlinkedType == LAGRANGIANLINEARRELATION)
        {}// Nothing to be done, blocks allready computed
        else if (relationType == LAGRANGIANRELATION || RlinkedType == LAGRANGIANRELATION)
        {
          computeExtraDiagonalBlocksLagrangianR(RCurrent, RLinked, sizeInteraction, linkedInteractionSize, commonDS, W, h, coupledInteractionsBlock);
          computeExtraDiagonalBlocksLagrangianR(RLinked, RCurrent, linkedInteractionSize, sizeInteraction, commonDS, W, h, coupledInteractionsBlockSym);
        }
        else RuntimeException::selfThrow("FrictionContact::computeAllBlocks not yet implemented for relation of type " + relationType);
      } // end of loop over linked interactions
    }
  } // end of loop over interactions -> increment current interaction

}

void FrictionContact::preFrictionContact(const double& time)
{
  IN("FrictionContact::preFrictionContact()\n");
  // compute M and q operators for FrictionContact problem

  // get topology
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // if relative degree is not equal to 0 or 1, check effective output
  if (! topology->isTimeInvariant())
  {
    updateBlocks(); // compute blocks for NonLinear Relations, with G that has been computed during update of previous time step
    computeEffectiveOutput();
    dim = topology->getEffectiveSizeOutput();
    if (dim != 0)
      assembleM();
  }

  if (dim != 0)
  {
    computeQ(time);
    // check z and w sizes and reset if necessary
    if (z == NULL)
    {
      z = new SimpleVector(dim);
      isZAllocatedIn = true;
    }
    else if (z->size() != dim)
    {
      // reset z if it has a wrong size
      if (isZAllocatedIn) delete z;
      z = new SimpleVector(dim);
      isZAllocatedIn = true;
    }

    if (w == NULL)
    {
      w = new SimpleVector(dim);
      isWAllocatedIn = true;
    }
    else if (w->size() != dim)
    {
      // reset w if it has a wrong size
      if (isWAllocatedIn) delete w;
      w = new SimpleVector(dim);
      isWAllocatedIn = true;
    }
    w->zero();
    z->zero();
  }
  OUT("FrictionContact::preFrictionContact()\n");
}

void FrictionContact::assembleM() //
{

  if (M == NULL)
  {
    M = new SiconosMatrix(dim, dim);
    isMAllocatedIn = true;
  }
  else if (M->size(0) != dim || M->size(1) != dim)
  {
    // reset M matrix if it has a wrong size
    if (isMAllocatedIn) delete M;
    M = new SiconosMatrix(dim, dim);
    isMAllocatedIn = true;
  }
  M->zero();

  // get topology
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  map< Interaction* , SiconosMatrix*>::iterator itDiago;
  //map< Interaction*, unsigned int>  interactionPositionMap =  topology->getInteractionPositionMap();

  unsigned int pos, col;
  unsigned int size, sizeLinked;

  map< Interaction*, unsigned int>  interactionEffectivePositionMap =  topology->getInteractionEffectivePositionMap();
  map< Interaction* , map<Interaction *, SiconosMatrix*> >::iterator itCoupled  ;
  map<Interaction *, SiconosMatrix*>::iterator it;
  SiconosMatrix *blockDiago , *coupledBlock;
  // --- loop over all the interactions ---
  for (itDiago = diagonalBlocksMap.begin(); itDiago != diagonalBlocksMap.end(); itDiago++)
  {
    // the current interaction
    Interaction * currentInteraction = itDiago->first;
    // and its corresponding diagonal block
    blockDiago = itDiago->second;
    pos = interactionEffectivePositionMap[ currentInteraction ];

    // get, if they exist, the linked interactions
    itCoupled = extraDiagonalBlocksMap.find(currentInteraction);

    if (topology->isTimeInvariant())
    {
      // diagonal blocks
      M->blockMatrixCopy(*blockDiago, pos, pos); // \todo avoid copy
      // extra-diagonal blocks
      if (itCoupled != extraDiagonalBlocksMap.end())
      {
        for (it = extraDiagonalBlocksMap[currentInteraction].begin(); it != extraDiagonalBlocksMap[currentInteraction].end(); it++)
        {
          coupledBlock = (*it).second;
          col = interactionEffectivePositionMap[(*it).first ];
          M->blockMatrixCopy(*coupledBlock, pos, col); // \todo avoid copy
        }
      }
    }
    else
    {
      // get the "reduced" diagonal block for the current interaction
      size = blockIndexesMap[currentInteraction].size();
      SiconosMatrix * reducedBlock = new SiconosMatrix(size, size);
      blockDiago->getBlock(blockIndexesMap[currentInteraction], blockIndexesMap[currentInteraction], *reducedBlock);

      // diagonal blocks
      M->blockMatrixCopy(*reducedBlock, pos, pos); // \todo avoid copy

      // extra-diagonal blocks
      SiconosMatrix * reducedCoupledBlock;

      // if currentInteraction is coupled with another interaction ...
      if (itCoupled != extraDiagonalBlocksMap.end())
      {
        // loop over linked interactions
        for (it = extraDiagonalBlocksMap[currentInteraction].begin(); it != extraDiagonalBlocksMap[currentInteraction].end(); it++)
        {
          // get list of effective relations for linked interaction
          Interaction * linkedInteraction = (*it).first;
          sizeLinked = blockIndexesMap[linkedInteraction].size();
          // get the corresponding "reduced" block
          coupledBlock = (*it).second;
          reducedCoupledBlock = new SiconosMatrix(size, sizeLinked);
          coupledBlock->getBlock(blockIndexesMap[currentInteraction], blockIndexesMap[linkedInteraction], *reducedCoupledBlock);
          col = interactionEffectivePositionMap[ linkedInteraction ];
          M->blockMatrixCopy(*reducedCoupledBlock, pos, col); // \todo avoid copy
          delete reducedCoupledBlock;
        }
      }
      delete reducedBlock;
    }
  }// --- end of interactions loop ---
}

void FrictionContact::compute(const double& time)
{
  IN("FrictionContact::compute(void)\n");

  // --- Prepare data for FrictionContact2D computing ---
  preFrictionContact(time);

  // --- Call Numerics solver ---
  if (dim != 0)
  {

    int info;
    int Dim = (int)dim;
    // get solving method and friction coefficient value.
    method solvingMethod = *(solver->getSolvingMethodPtr());
    Interaction * currentInteraction = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0);
    // call Numerics method for 2D or 3D problem:
    if (nspbType == "FrictionContact2D")
    {
      solvingMethod.pfc_2D.mu = static_cast<NewtonImpactFrictionNSL*>(currentInteraction->getNonSmoothLawPtr())->getMu();
      info = pfc_2D_solver(M->getArray(), q->getArray(), &Dim, &solvingMethod  , z->getArray(), w->getArray());
    }
    else if (nspbType == "FrictionContact3D")
    {
      Dim = Dim / 3; // in pfc_3D, Dim is the number of contact points.
      solvingMethod.pfc_3D.mu = static_cast<NewtonImpactFrictionNSL*>(currentInteraction->getNonSmoothLawPtr())->getMu();
      info = pfc_3D_solver(M->getArray(), q->getArray(), &Dim, &solvingMethod  , z->getArray(), w->getArray());
    }
    else
      RuntimeException::selfThrow("FrictionContact::compute, unknown or unconsistent non smooth problem type: " + nspbType);
    check_solver(info);
    // --- Recover the desired variables from FrictionContact2D output ---
    postFrictionContact(*w, *z);
  }

  OUT("FrictionContact::compute(void)\n");
}

void FrictionContact::postFrictionContact(const SimpleVector& w, const SimpleVector &z)
{
  // --- get topology ---
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // --- get interactions list ---
  vector<Interaction*> listInteractions = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
  vector<Interaction*>::iterator iter;
  // get Interaction position map
  map< Interaction* , SiconosMatrix*>::iterator itDiago;
  map< Interaction*, unsigned int>  interactionEffectivePositionMap =  topology->getInteractionEffectivePositionMap();

  vector<unsigned int> effectiveIndexes, indexMin;
  unsigned int effectivePosition;
  unsigned int effectiveSize, numberOfRelations ;
  unsigned int i; // index of derivation for Y
  unsigned int j; // number of relation in a specific interaction
  unsigned int k;
  unsigned int rMax; // maximum value for relative degrees
  SimpleVector * tmpY, *tmpLambda;
  vector< SimpleVector* >  Y, Lambda;
  vector<unsigned int> r; // relative degrees

  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // the current interaction and its size
    Interaction *currentInteraction = *iter;
    numberOfRelations = currentInteraction->getNInteraction();

    // Y vector of the interactions
    Y = currentInteraction -> getY();
    // lambda vector
    Lambda = currentInteraction ->getLambda();

    // relative degrees
    r = topology->getRelativeDegrees(currentInteraction);
    rMax = r[0]; // we suppose the interaction is homogeneous
    if (rMax == 0) rMax = 1 ; // warning: review r=0 case

    // get the list of relations that are constrained and the position vector
    effectiveSize = topology->computeEffectiveSizeOutput(currentInteraction) ; // Improvement: save this value in Topology?
    effectiveIndexes = topology->getEffectiveIndexes(currentInteraction);
    indexMin = topology->getIndexMin(currentInteraction);
    // 'offset' due to indexMin
    /*      for(j=0;j<effectiveSize;j++)
    effectiveIndexes[j] += indexMin[effectiveIndexes[j]]*numberOfRelations;
    */
    // Get the 'position' of vector corresponding to the current interaction, in the FrictionContact output:
    effectivePosition = topology->getInteractionEffectivePosition(currentInteraction);

    // we consider that the interaction is homogeneous, ie all degrees are equals
    tmpY = new SimpleVector(effectiveSize); // warning: review r=0 case
    tmpLambda = new SimpleVector(effectiveSize);

    unsigned int pos;
    vector<unsigned int>::iterator itPos;
    // First we get in w results corresponding to the current interaction

    for (j = 0; j < effectiveSize ; j++)
    {
      (*tmpY)(j) =  w(j + effectivePosition);
      (*tmpLambda)(j) = z(j + effectivePosition);
    }

    // then we save these results in Y and Lambda, only for effective relations
    for (i = 0; i < rMax ; i++)
    {
      for (j = 0; j < numberOfRelations; j++)
      {
        pos = j * rMax + i;
        for (k = 0; k < effectiveSize; k++)
        {
          // itPos = find(effectiveIndexes.begin(),effectiveIndexes.end(),pos);
          // -> how can we get k/ effectiveIndex(k)=itPos ????
          //if (itPos!=effectiveIndexes.end())
          if (effectiveIndexes[k] == pos)
          {
            (*(Y[i]))(j) = (*tmpY)(k);
            (*(Lambda[i]))(j) = (*tmpLambda)(k);
          }
        }
      }
    }
    delete tmpY;
    delete tmpLambda;
  }
}


void FrictionContact::display() const
{
  cout << "======= FrictionContact display ======" << endl;
  cout << "| dim : " << dim << endl;
  cout << "| FrictionContact Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| FrictionContact vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;

}

void FrictionContact::saveNSProblemToXML()
{
  IN("FrictionContact::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<FrictionContactXML*>(onestepnspbxml))->setM(*M);
    (static_cast<FrictionContactXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("FrictionContact::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("FrictionContact::saveNSProblemToXML\n");
}

void FrictionContact::saveMToXML()
{
  IN("FrictionContact::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<FrictionContactXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("FrictionContact::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("FrictionContact::saveMToXML\n");
}

void FrictionContact::saveQToXML()
{
  IN("FrictionContact::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<FrictionContactXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("FrictionContact::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("FrictionContact::saveQToXML\n");
}

FrictionContact* FrictionContact::convert(OneStepNSProblem* osnsp)
{
  FrictionContact* fc2d = dynamic_cast<FrictionContact*>(osnsp);
  return fc2d;
}



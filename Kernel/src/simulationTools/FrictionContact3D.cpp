/* Siconos version 1.0, Copyright INRIA 2005.
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
#include "FrictionContact3D.h"

// includes to be deleted thanks to factories
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactFrictionNSL.h"
#include "LinearTIR.h"

using namespace std;

// Default constructor
FrictionContact3D::FrictionContact3D(): OneStepNSProblem(), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = FrictionContact3D_OSNSP;
}

// xml constructor
FrictionContact3D::FrictionContact3D(OneStepNSProblemXML* osNsPbXml, Strategy* newStrat):
  OneStepNSProblem(osNsPbXml, newStrat), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = FrictionContact3D_OSNSP;
  if (osNsPbXml != NULL)
  {
    FrictionContact3DXML * xmllcp = (static_cast<FrictionContact3DXML*>(osNsPbXml));

    // If both q and M are given in xml file, check if sizes are consistent
    if (xmllcp->hasQ() && xmllcp->hasM() && ((xmllcp->getM()).size(0) != (xmllcp->getQ()).size()))
      RuntimeException::selfThrow("FrictionContact3D: xml constructor, inconsistent sizes between given q and M");

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
  else RuntimeException::selfThrow("FrictionContact3D: xml constructor, xml file=NULL");
}

// Constructor from a set of data
FrictionContact3D::FrictionContact3D(Strategy* newStrat, const string& newSolver, const string& newSolvingMethod,
                                     const int& MaxIter, const double & Tolerance, const string & NormType,
                                     const double & SearchDirection):
  OneStepNSProblem(newStrat, newSolver, newSolvingMethod, MaxIter, Tolerance, NormType, SearchDirection),
  w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = FrictionContact3D_OSNSP;
}

// Constructor from a set of data
FrictionContact3D::FrictionContact3D(Strategy* newStrat, Solver*  newSolver):
  OneStepNSProblem(newStrat, newSolver), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = FrictionContact3D_OSNSP;
}

// destructor
FrictionContact3D::~FrictionContact3D()
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

void FrictionContact3D::setW(const SimpleVector& newValue)
{
  if (dim != newValue.size())
    RuntimeException::selfThrow("FrictionContact3D: setW, inconsistent size between given w size and problem size. You should set dim before");

  if (w == NULL)
  {
    w = new SimpleVector(dim);
    isWAllocatedIn = true;
  }
  else if (w->size() != dim)
    RuntimeException::selfThrow("FrictionContact3D: setW, w size differs from dim");

  *w = newValue;
}

void FrictionContact3D::setWPtr(SimpleVector* newPtr)
{
  if (dim != newPtr->size())
    RuntimeException::selfThrow("FrictionContact3D: setWPtr, inconsistent size between given w size and problem size. You should set dim before");

  if (isWAllocatedIn) delete w;
  w = newPtr;
  isWAllocatedIn = false;
}


void FrictionContact3D::setZ(const SimpleVector& newValue)
{
  if (dim != newValue.size())
    RuntimeException::selfThrow("FrictionContact3D: setZ, inconsistent size between given z size and problem size. You should set dim before");

  if (z == NULL)
  {
    z = new SimpleVector(dim);
    isZAllocatedIn = true;
  }

  *z = newValue;
}

void FrictionContact3D::setZPtr(SimpleVector* newPtr)
{
  if (dim != newPtr->size())
    RuntimeException::selfThrow("FrictionContact3D: setZPtr, inconsistent size between given z size and problem size. You should set dim before");

  if (isZAllocatedIn) delete z;
  z = newPtr;
  isZAllocatedIn = false;
}

void FrictionContact3D::setM(const SiconosMatrix& newValue)
{
  if (dim != newValue.size(0) || dim != newValue.size(1))
    RuntimeException::selfThrow("FrictionContact3D: setM, inconsistent size between given M size and problem size. You should set dim before");

  if (M == NULL)
  {
    M = new SiconosMatrix(dim, dim);
    isMAllocatedIn = true;
  }

  *M = newValue;
}

void FrictionContact3D::setMPtr(SiconosMatrix* newPtr)
{
  if (dim != newPtr->size(0) || dim != newPtr->size(1))
    RuntimeException::selfThrow("FrictionContact3D: setMPtr, inconsistent size between given M size and problem size. You should set dim before");

  if (isMAllocatedIn) delete M;
  M = newPtr;
  isMAllocatedIn = false;
}


void FrictionContact3D::setQ(const SimpleVector& newValue)
{
  if (dim != newValue.size())
    RuntimeException::selfThrow("FrictionContact3D: setQ, inconsistent size between given q size and problem size. You should set dim before");

  if (q == NULL)
  {
    q = new SimpleVector(dim);
    isQAllocatedIn = true;
  }

  *q = newValue;
}

void FrictionContact3D::setQPtr(SimpleVector* newPtr)
{
  if (dim != newPtr->size())
    RuntimeException::selfThrow("FrictionContact3D: setQPtr, inconsistent size between given q size and problem size. You should set dim before");

  if (isQAllocatedIn) delete q;
  q = newPtr;
  isQAllocatedIn = false;
}

void FrictionContact3D::initialize()
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

void FrictionContact3D::computeAllBlocks()
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
        RuntimeException::selfThrow("FrictionContact3D::computeAllBlocks not yet implemented for Integrator of type " + Osi->getType());
    }

    // get the relation of the current interaction and its type
    Relation * RCurrent = currentInteraction -> getRelationPtr();
    relationType = RCurrent->getType();

    if (relationType == LINEARTIRELATION)
      computeDiagonalBlocksLinearTIR(RCurrent, sizeInteraction, vDS, W, h, currentMatrixBlock);
    else if (relationType == LAGRANGIANLINEARRELATION || relationType == LAGRANGIANRELATION)
      computeDiagonalBlocksLagrangianR(RCurrent, sizeInteraction, vDS, W, h, currentMatrixBlock);
    //else if (relationType == LAGRANGIANRELATION)
    //computeDiagonalBlocksLagrangianR(RCurrent, sizeInteraction,vDS,W,h,currentMatrixBlock);
    else RuntimeException::selfThrow("FrictionContact3D::computeAllBlocks not yet implemented for relation of type " + relationType);

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
        if (relationType == LINEARTIRELATION && RlinkedType == LINEARTIRELATION)
        {
          computeExtraDiagonalBlocksLinearTIR(RCurrent, RLinked, sizeInteraction, linkedInteractionSize, commonDS, W, h, coupledInteractionsBlock);
          computeExtraDiagonalBlocksLinearTIR(RLinked, RCurrent, linkedInteractionSize, sizeInteraction, commonDS, W, h, coupledInteractionsBlockSym);
        }
        //else if (relationType == LAGRANGIANLINEARRELATION && RlinkedType== LAGRANGIANLINEARRELATION)
        //computeExtraDiagonalBlocksLagrangianLinearR(RCurrent, RLinked, sizeInteraction, linkedInteractionSize, commonDS,W,h,coupledInteractionsBlock);
        else if ((relationType == LAGRANGIANRELATION || relationType == LAGRANGIANLINEARRELATION) && (RlinkedType == LAGRANGIANRELATION ||  RlinkedType == LAGRANGIANLINEARRELATION))
        {
          computeExtraDiagonalBlocksLagrangianR(RCurrent, RLinked, sizeInteraction, linkedInteractionSize, commonDS, W, h, coupledInteractionsBlock);
          computeExtraDiagonalBlocksLagrangianR(RLinked, RCurrent, linkedInteractionSize, sizeInteraction, commonDS, W, h, coupledInteractionsBlockSym);
        }
        else RuntimeException::selfThrow("FrictionContact3D::computeAllBlocks not yet implemented for relation of type " + relationType);
      } // end of loop over linked interactions
    }
  } // end of loop over interactions -> increment current interaction

}

void FrictionContact3D::computeDiagonalBlocksLinearTIR(Relation * R, const unsigned int& sizeInteraction, vector<DynamicalSystem*> vDS,
    map<DynamicalSystem*, SiconosMatrix*> W, const double& h, SiconosMatrix* currentMatrixBlock)
{
  // convert to linearTIR
  LinearTIR *LTIR = static_cast<LinearTIR*>(R);
  // For 1 relation, we need:
  // - the corresponding row of D
  // - a block-row of C, that corresponds to a specific DS
  // - a block of B, that corresponds to a specific DS
  SiconosMatrix* D = LTIR->getDPtr();
  SiconosMatrix * C, *B;
  vector<DynamicalSystem*>::iterator itDS;
  unsigned int sizeDS;

  bool isHomogeneous = false; // \todo to be defined with relative degrees
  // isHomogeneous = true means that all relative degrees of the interaction are equal.

  if (!isHomogeneous)
  {
    // the resulting row-block
    SimpleVector * currentLine = new SimpleVector(sizeInteraction);
    // a row of D corresponding to the a relation of the current interaction
    SimpleVector * Drow = new SimpleVector(sizeInteraction);
    // a row of C corresponding to a specific DS (its size depends on the DS size)
    SimpleVector * Crow ;

    // compute currentLine = Drow + h sum(j) Crow,j Wj  Bj, with j the list of DS of the
    // current interaction.
    // loop over the relations of the current interaction
    for (unsigned int i = 0; i < sizeInteraction; i++)
    {
      // get row i of D
      D->getRow(i, *Drow);
      *currentLine = *Drow;
      // loop over the DS of the current interaction
      for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
      {
        sizeDS = (*itDS)->getN();
        // get blocks corresponding to the current DS
        C = new SiconosMatrix(sizeInteraction, sizeDS);
        B = new SiconosMatrix(sizeDS, sizeInteraction);
        LTIR->getCBlockDSPtr(*itDS, *C);
        LTIR->getBBlockDSPtr(*itDS, *B);
        // get row i of C
        Crow = new SimpleVector(sizeDS);
        C->getRow(i, *Crow);

        // compute currentLine
        *currentLine +=  h* *Crow * (*W[*itDS]* *B);
        delete C;
        delete B;
        delete Crow;
      }
      // save the result in currentMatrixBlock (and so in diagonalBlocksMap)
      currentMatrixBlock->setRow(i, *currentLine);
    }
    delete currentLine;
    delete Drow;
  }
  else
  {
    // the resulting block
    *currentMatrixBlock = *D;
    // loop over the DS of the current interaction
    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
    {
      sizeDS = (*itDS)->getN();
      // get blocks corresponding to the current DS
      C = new SiconosMatrix(sizeInteraction, sizeDS);
      B = new SiconosMatrix(sizeDS, sizeInteraction);
      LTIR->getCBlockDSPtr(*itDS, *C);
      LTIR->getBBlockDSPtr(*itDS, *B);
      *currentMatrixBlock += h * *C * (*W[*itDS] * *B);
      delete C;
      delete B;
    }
  }
}

void FrictionContact3D::computeExtraDiagonalBlocksLinearTIR(Relation * RCurrent, Relation* RLinked, const unsigned int& sizeInteraction,
    const unsigned int& linkedInteractionSize, vector<DynamicalSystem*> commonDS,
    map<DynamicalSystem*, SiconosMatrix*> W, const double& h, SiconosMatrix* coupledInteractionsBlock)
{
  // convert to linearTIR
  LinearTIR *LTIR1 = static_cast<LinearTIR*>(RCurrent);
  LinearTIR *LTIR2 = static_cast<LinearTIR*>(RLinked);
  SiconosMatrix *C, *B;
  coupledInteractionsBlock->zero();
  vector<DynamicalSystem*>::iterator itDS;
  unsigned int sizeDS;
  bool isHomogeneous = false; // \todo to be defined with relative degrees
  if (!isHomogeneous)
  {
    // the resulting row-block
    SimpleVector * currentLine = new SimpleVector(sizeInteraction);
    // a row of C corresponding to a specific DS (its size depends on the DS size)
    SimpleVector * Crow ;
    currentLine->zero();

    // compute currentLine =  h sum(j) Crow,j Wj  Bj, with j the list of common DS
    // C belongs to current Interaction and B to linked interaction
    // loop over the relations of the current interaction
    for (unsigned int i = 0; i < sizeInteraction; i++)
    {
      // loop over common DS
      for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
      {
        sizeDS = (*itDS)->getN();
        // get blocks corresponding to the current DS
        C = new SiconosMatrix(sizeInteraction, sizeDS);
        B = new SiconosMatrix(sizeDS, linkedInteractionSize);
        LTIR1->getCBlockDSPtr(*itDS, *C);
        LTIR2->getBBlockDSPtr(*itDS, *B);
        // get row i of C
        Crow = new SimpleVector(sizeDS);
        C->getRow(i, *Crow);
        // compute currentLine
        *currentLine +=   *Crow * (*W[*itDS]* *B);
        delete Crow;
        delete C;
        delete B;
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
      sizeDS = (*itDS)->getN();
      // get blocks corresponding to the current DS
      C = new SiconosMatrix(sizeInteraction, sizeDS);
      B = new SiconosMatrix(sizeDS, linkedInteractionSize);
      LTIR1->getCBlockDSPtr(*itDS, *C);
      LTIR2->getBBlockDSPtr(*itDS, *B);
      *coupledInteractionsBlock +=  *C * (*W[*itDS] * *B);
      delete C;
      delete B;
    }
  }
}

void FrictionContact3D::computeDiagonalBlocksLagrangianR(Relation * R, const unsigned int& sizeInteraction, vector<DynamicalSystem*> vDS,
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

void FrictionContact3D::computeExtraDiagonalBlocksLagrangianR(Relation * RCurrent, Relation* RLinked, const unsigned int& sizeInteraction,
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

void FrictionContact3D::updateBlocks()
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
        RuntimeException::selfThrow("FrictionContact3D::computeAllBlocks not yet implemented for Integrator of type " + Osi->getType());
    }

    // get the relation of the current interaction and its type
    Relation * RCurrent = currentInteraction -> getRelationPtr();
    relationType = RCurrent->getType();

    if (relationType == LAGRANGIANLINEARRELATION)
    {}// Nothing to be done, blocks allready computed
    else if (relationType == LAGRANGIANRELATION)
      computeDiagonalBlocksLagrangianR(RCurrent, sizeInteraction, vDS, W, h, currentMatrixBlock);
    else RuntimeException::selfThrow("FrictionContact3D::updateBlocks not yet implemented for relation of type " + relationType);

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
        else RuntimeException::selfThrow("FrictionContact3D::computeAllBlocks not yet implemented for relation of type " + relationType);
      } // end of loop over linked interactions
    }
  } // end of loop over interactions -> increment current interaction

}

void FrictionContact3D::preFrictionContact3D(const double& time)
{
  IN("FrictionContact3D::preFrictionContact3D()\n");
  // compute M and q operators for FrictionContact3D problem

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
  OUT("FrictionContact3D::preFrictionContact3D()\n");
}

void FrictionContact3D::assembleM() //
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

void FrictionContact3D::computeQ(const double& time)
{
  IN("FrictionContact3D::computeQ(void)\n");
  if (q == NULL)
  {
    q = new SimpleVector(dim);
    isQAllocatedIn = true;
  }
  else if (q->size() != dim)
  {
    // reset q if it has a wrong size
    if (isQAllocatedIn) delete q;
    q = new SimpleVector(dim);
    isQAllocatedIn = true;
  }

  q->zero();

  // --- get topology ---
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // --- get interactions list ---
  vector<Interaction*> listInteractions = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
  vector<Interaction*>::iterator iter;
  // get Interaction position map
  map< Interaction* , SiconosMatrix*>::iterator itDiago;
  map< Interaction*, unsigned int>  interactionEffectivePositionMap =  topology->getInteractionEffectivePositionMap();
  unsigned int pos;
  Relation *R;
  NonSmoothLaw *nslaw;
  vector<unsigned int> index ;
  unsigned int effectiveSize ;
  Interaction *currentInteraction ;
  SimpleVector * yFree;
  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // get current interaction, its relation and its nslaw
    currentInteraction = *iter;
    unsigned int numberOfRelations = currentInteraction->getNInteraction();
    R = currentInteraction->getRelationPtr();
    string relationType = R->getType();
    nslaw = currentInteraction->getNonSmoothLawPtr();
    // position of current yFree in global (including all interactions) y vector
    pos = interactionEffectivePositionMap[currentInteraction];
    index = topology->getEffectiveIndexes(currentInteraction);
    if (relationType == LINEARTIRELATION)
    {
      LinearTIR *LTIR = static_cast<LinearTIR*>(R);
      string nslawType = nslaw->getType() ;
      if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
      {
        LTIR->computeFreeOutput(); // free output is saved in y
        if (topology->isTimeInvariant())
        {
          yFree = new SimpleVector(currentInteraction->getY(0));   // copy, no pointer equality
          for (unsigned int i = 0; i < numberOfRelations; i++)
            (*q)(i + pos) = (*yFree)(i);
        }
        else
        {
          effectiveSize = index.size();
          yFree = new SimpleVector(effectiveSize); // we get only the "effective" relations
          (currentInteraction->getY(0)).getBlock(index, *yFree); // copy, no pointer equality
          for (unsigned int i = 0; i < effectiveSize; i++)
            (*q)(i + pos) = (*yFree)(i);
        }
      }
      else
        RuntimeException::selfThrow("FrictionContact3D::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else if (relationType == LAGRANGIANLINEARRELATION || relationType == LAGRANGIANRELATION)
    {
      LagrangianLinearR *LLR = static_cast<LagrangianLinearR*>(R);
      if (nslaw->getType() == NEWTONIMPACTFRICTIONNSLAW)
      {
        vector<unsigned int> indexMax = topology->getIndexMax(currentInteraction);

        NewtonImpactFrictionNSL * newton = static_cast<NewtonImpactFrictionNSL*>(nslaw);
        double e = newton->getEn();
        //double mu= newton->getMu();
        LLR->computeFreeOutput(time);
        if (topology->isTimeInvariant())
        {
          yFree = new SimpleVector(currentInteraction -> getY(1)); // copy, no pointer equality
          // Only normal part has to be multiplied by e.
          // even indexes corresponds to normal components of y.

          SimpleVector * yDotOld = currentInteraction->getYOldPtr(1);
          for (unsigned int i = 0; i < numberOfRelations; i = i + 3)
            (*yFree)(i) += e * (*yDotOld)(i);

          for (unsigned int i = 0; i < numberOfRelations; i++)
            (*q)(i + pos) = (*yFree)(i);
        }
        else
        {
          // compute list of effective relations
          unsigned int k = 0;
          for (unsigned int j = 0; j < numberOfRelations; j++)
          {
            for (unsigned int i = 1; i < indexMax[j] + 1; i++)
            {
              index[k] = j;
              k++;
            }
          }

          effectiveSize = index.size();
          yFree = new SimpleVector(effectiveSize); // we get only the "effective" relations
          (currentInteraction->getY(1)).getBlock(index, *yFree); // copy, no pointer equality
          SimpleVector * tmp =  new SimpleVector(effectiveSize);
          (currentInteraction -> getYOld(1)).getBlock(index, *tmp);

          for (unsigned int i = 0; i < numberOfRelations; i = i + 3)
            (*yFree)(i) += e * (*tmp)(i);

          delete tmp;
          for (unsigned int i = 0; i < effectiveSize; i++)
            (*q)(i + pos) = (*yFree)(i);
        }
      }
      else
        RuntimeException::selfThrow("FrictionContact3D::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else
      RuntimeException::selfThrow("FrictionContact3D::computeQ not yet implemented for relation of type " + R->getType());

    delete yFree;
  }
  OUT("FrictionContact3D::computeQ(void)\n");
}

void FrictionContact3D::compute(const double& time)
{
  IN("FrictionContact3D::compute(void)\n");

  // --- Prepare data for FrictionContact3D computing ---
  preFrictionContact3D(time);

  // --- Call Numerics solver ---
  if (dim != 0)
  {

    int info;
    int Dim = (int)dim;
    method solvingMethod = *(solver->getSolvingMethodPtr());
    Interaction * currentInteraction = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractionPtr(0);
    solvingMethod.pfc_3D.mu = static_cast<NewtonImpactFrictionNSL*>(currentInteraction->getNonSmoothLawPtr())->getMu();
    info = pfc_3D_solver(M->getArray(), q->getArray(), &Dim, &solvingMethod  , z->getArray(), w->getArray());
    check_solver(info);
    // --- Recover the desired variables from FrictionContact3D output ---
    postFrictionContact3D(*w, *z);
  }

  OUT("FrictionContact3D::compute(void)\n");
}

void FrictionContact3D::postFrictionContact3D(const SimpleVector& w, const SimpleVector &z)
{
  // --- get topology ---
  cout << " IN POST " << endl;
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
    // Get the 'position' of vector corresponding to the current interaction, in the FrictionContact3D output:
    effectivePosition = topology->getInteractionEffectivePosition(currentInteraction);
    cout << " IN POST 3" << endl;

    // we consider that the interaction is homogeneous, ie all degrees are equals
    cout << "EF " << effectiveSize << endl;
    tmpY = new SimpleVector(effectiveSize); // warning: review r=0 case
    cout << " IN POST 33" << endl;
    tmpLambda = new SimpleVector(effectiveSize);

    unsigned int pos;
    vector<unsigned int>::iterator itPos;
    // First we get in w results corresponding to the current interaction
    cout << " IN POST 4" << endl;

    for (j = 0; j < effectiveSize ; j++)
    {
      (*tmpY)(j) =  w(j + effectivePosition);
      (*tmpLambda)(j) = z(j + effectivePosition);
    }
    cout << " IN POST 5" << endl;

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
  cout << " IN POST FIN" << endl;
}


void FrictionContact3D::display() const
{
  cout << "======= FrictionContact3D display ======" << endl;
  cout << "| dim : " << dim << endl;
  cout << "| FrictionContact3D Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| FrictionContact3D vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;

}

void FrictionContact3D::saveNSProblemToXML()
{
  IN("FrictionContact3D::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<FrictionContact3DXML*>(onestepnspbxml))->setM(*M);
    (static_cast<FrictionContact3DXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("FrictionContact3D::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("FrictionContact3D::saveNSProblemToXML\n");
}

void FrictionContact3D::saveMToXML()
{
  IN("FrictionContact3D::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<FrictionContact3DXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("FrictionContact3D::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("FrictionContact3D::saveMToXML\n");
}

void FrictionContact3D::saveQToXML()
{
  IN("FrictionContact3D::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<FrictionContact3DXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("FrictionContact3D::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("FrictionContact3D::saveQToXML\n");
}

FrictionContact3D* FrictionContact3D::convert(OneStepNSProblem* osnsp)
{
  cout << "FrictionContact3D::convert (OneStepNSProblem* osnsp)" << endl;
  FrictionContact3D* fc3d = dynamic_cast<FrictionContact3D*>(osnsp);
  return fc3d;
}



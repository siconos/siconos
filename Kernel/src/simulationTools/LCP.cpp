/* Siconos-Kernel version 1.1.3, Copyright INRIA 2005-2006.
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
#include "LCP.h"

// includes to be deleted thanks to factories
#include "Moreau.h"
#include "LagrangianLinearR.h"
#include "NewtonImpactLawNSL.h"
#include "LinearTIR.h"

using namespace std;

// Default constructor
LCP::LCP(): OneStepNSProblem(), nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = "LCP";
}

// xml constructor
LCP::LCP(OneStepNSProblemXML* onestepnspbxml, Strategy* newStrat):
  OneStepNSProblem(onestepnspbxml, newStrat), nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = "LCP";
  // === read solver related data ===
  if (onestepnspbxml->hasSolver())
    solver = new Solver(onestepnspbxml->getSolverXMLPtr(), nspbType);
  else // solver = default one
  {
    solver = new Solver(nspbType, DEFAULT_SOLVER, DEFAULT_ITER, DEFAULT_TOL, DEFAULT_NORMTYPE, DEFAULT_SEARCHDIR);
    cout << " Warning, no solver defined in non-smooth problem - Default one is used" << endl;
    solver->display();
  }
  isSolverAllocatedIn = true;

  if (onestepnspbxml != NULL)
  {
    LCPXML * xmllcp = (static_cast<LCPXML*>(onestepnspbxml));

    // If both q and M are given in xml file, check if sizes are consistent
    if (xmllcp->hasQ() && xmllcp->hasM() && ((xmllcp->getM()).size(0) != (xmllcp->getQ()).size()))
      RuntimeException::selfThrow("LCP: xml constructor, inconsistent sizes between given q and M");

    // first nlcp is given by M matrix size in xml file
    if (xmllcp->hasM())
    {
      nLcp = (xmllcp->getM()).size(0);
      M = new SimpleMatrix(xmllcp->getM());
      isMAllocatedIn = true;
    }

    if (xmllcp->hasQ())
    {
      // get nLcp if necessary
      if (M == NULL)
        nLcp = (xmllcp->getQ()).size();

      q = new SimpleVector(xmllcp->getQ());
      isQAllocatedIn = true;
    }

  }
  else RuntimeException::selfThrow("LCP: xml constructor, xml file=NULL");
}

// Constructor from a set of data
LCP::LCP(Strategy* newStrat, const string& newSolver, const unsigned int& MaxIter,
         const double & Tolerance, const string & NormType, const double & SearchDirection):
  OneStepNSProblem(newStrat),
  nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = "LCP";
  // set solver:
  solver = new Solver(nspbType, newSolver, MaxIter, Tolerance, NormType, SearchDirection);
  isSolverAllocatedIn = true;
}

// destructor
LCP::~LCP()
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

void LCP::setW(const SimpleVector& newValue)
{
  if (nLcp != newValue.size())
    RuntimeException::selfThrow("LCP: setW, inconsistent size between given w size and problem size. You should set nLcp before");

  if (w == NULL)
  {
    w = new SimpleVector(nLcp);
    isWAllocatedIn = true;
  }
  else if (w->size() != nLcp)
    RuntimeException::selfThrow("LCP: setW, w size differs from nLcp");

  *w = newValue;
}

void LCP::setWPtr(SimpleVector* newPtr)
{
  if (nLcp != newPtr->size())
    RuntimeException::selfThrow("LCP: setWPtr, inconsistent size between given w size and problem size. You should set nLcp before");

  if (isWAllocatedIn) delete w;
  w = newPtr;
  isWAllocatedIn = false;
}


void LCP::setZ(const SimpleVector& newValue)
{
  if (nLcp != newValue.size())
    RuntimeException::selfThrow("LCP: setZ, inconsistent size between given z size and problem size. You should set nLcp before");

  if (z == NULL)
  {
    z = new SimpleVector(nLcp);
    isZAllocatedIn = true;
  }

  *z = newValue;
}

void LCP::setZPtr(SimpleVector* newPtr)
{
  if (nLcp != newPtr->size())
    RuntimeException::selfThrow("LCP: setZPtr, inconsistent size between given z size and problem size. You should set nLcp before");

  if (isZAllocatedIn) delete z;
  z = newPtr;
  isZAllocatedIn = false;
}

void LCP::setM(const SiconosMatrix& newValue)
{
  if (nLcp != newValue.size(0) || nLcp != newValue.size(1))
    RuntimeException::selfThrow("LCP: setM, inconsistent size between given M size and problem size. You should set nLcp before");

  if (M == NULL)
  {
    M = new SimpleMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }

  *M = newValue;
}

void LCP::setMPtr(SiconosMatrix* newPtr)
{
  if (nLcp != newPtr->size(0) || nLcp != newPtr->size(1))
    RuntimeException::selfThrow("LCP: setMPtr, inconsistent size between given M size and problem size. You should set nLcp before");

  if (isMAllocatedIn) delete M;
  M = newPtr;
  isMAllocatedIn = false;
}


void LCP::setQ(const SimpleVector& newValue)
{
  if (nLcp != newValue.size())
    RuntimeException::selfThrow("LCP: setQ, inconsistent size between given q size and problem size. You should set nLcp before");

  if (q == NULL)
  {
    q = new SimpleVector(nLcp);
    isQAllocatedIn = true;
  }

  *q = newValue;
}

void LCP::setQPtr(SimpleVector* newPtr)
{
  if (nLcp != newPtr->size())
    RuntimeException::selfThrow("LCP: setQPtr, inconsistent size between given q size and problem size. You should set nLcp before");

  if (isQAllocatedIn) delete q;
  q = newPtr;
  isQAllocatedIn = false;
}

void LCP::initialize()
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
    nLcp = topology->getEffectiveSizeOutput();
    assembleM();
  }
}

void LCP::computeAllBlocks()
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
    diagonalBlocksMap[ currentInteraction ] = new SimpleMatrix(sizeBlock, sizeBlock);
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
        RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for Integrator of type " + Osi->getType());
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
    else RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for relation of type " + relationType);

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

        (extraDiagonalBlocksMap[currentInteraction])[linkedInteraction] = new SimpleMatrix(sizeBlock, sizeLinkedBlock);
        (extraDiagonalBlocksMap[linkedInteraction])[currentInteraction] = new SimpleMatrix(sizeLinkedBlock, sizeBlock);
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
        else RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for relation of type " + relationType);
      } // end of loop over linked interactions
    }
  } // end of loop over interactions -> increment current interaction

}

void LCP::computeDiagonalBlocksLinearTIR(Relation * R, const unsigned int& sizeInteraction, vector<DynamicalSystem*> vDS,
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
        C = new SimpleMatrix(sizeInteraction, sizeDS);
        B = new SimpleMatrix(sizeDS, sizeInteraction);
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
      C = new SimpleMatrix(sizeInteraction, sizeDS);
      B = new SimpleMatrix(sizeDS, sizeInteraction);
      LTIR->getCBlockDSPtr(*itDS, *C);
      LTIR->getBBlockDSPtr(*itDS, *B);
      *currentMatrixBlock += h * *C * (*W[*itDS] * *B);
      delete C;
      delete B;
    }
  }
}

void LCP::computeExtraDiagonalBlocksLinearTIR(Relation * RCurrent, Relation* RLinked, const unsigned int& sizeInteraction,
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
        C = new SimpleMatrix(sizeInteraction, sizeDS);
        B = new SimpleMatrix(sizeDS, linkedInteractionSize);
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
      C = new SimpleMatrix(sizeInteraction, sizeDS);
      B = new SimpleMatrix(sizeDS, linkedInteractionSize);
      LTIR1->getCBlockDSPtr(*itDS, *C);
      LTIR2->getBBlockDSPtr(*itDS, *B);
      *coupledInteractionsBlock +=  *C * (*W[*itDS] * *B);
      delete C;
      delete B;
    }
  }
}

void LCP::computeDiagonalBlocksLagrangianR(Relation * R, const unsigned int& sizeInteraction, vector<DynamicalSystem*> vDS,
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
        G = new SimpleMatrix(sizeInteraction, sizeDS);
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
      G = new SimpleMatrix(sizeInteraction, sizeDS);
      LR->getGBlockDS(*itDS, *G);
      *currentMatrixBlock +=  *G * (W[*itDS])->multTranspose(*G);
      delete G;
    }
  }
}

void LCP::computeExtraDiagonalBlocksLagrangianR(Relation * RCurrent, Relation* RLinked, const unsigned int& sizeInteraction,
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
        Gcurrent = new SimpleMatrix(sizeInteraction, sizeDS);
        Glinked = new SimpleMatrix(linkedInteractionSize, sizeDS);
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
      Gcurrent = new SimpleMatrix(sizeInteraction, sizeDS);
      Glinked = new SimpleMatrix(linkedInteractionSize, sizeDS);
      LR1->getGBlockDS(*itDS, *Gcurrent);
      LR2->getGBlockDS(*itDS, *Glinked);
      *coupledInteractionsBlock += *Gcurrent * (W[*itDS])->multTranspose(*Glinked);
      delete Gcurrent;
      delete Glinked;
    }
  }
}

void LCP::updateBlocks()
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
        RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for Integrator of type " + Osi->getType());
    }

    // get the relation of the current interaction and its type
    Relation * RCurrent = currentInteraction -> getRelationPtr();
    relationType = RCurrent->getType();

    if (relationType == LAGRANGIANLINEARRELATION)
    {}// Nothing to be done, blocks allready computed
    else if (relationType == LAGRANGIANRELATION)
      computeDiagonalBlocksLagrangianR(RCurrent, sizeInteraction, vDS, W, h, currentMatrixBlock);
    else RuntimeException::selfThrow("LCP::updateBlocks not yet implemented for relation of type " + relationType);

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
        else RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for relation of type " + relationType);
      } // end of loop over linked interactions
    }
  } // end of loop over interactions -> increment current interaction

}

void LCP::preLCP(const double& time)
{
  IN("LCP::preLCP()\n");
  // compute M and q operators for LCP problem

  // get topology
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // if relative degree is not equal to 0 or 1, check effective output
  if (! topology->isTimeInvariant())
  {
    updateBlocks(); // compute blocks for NonLinear Relations, with G that has been computed during update of previous time step
    computeEffectiveOutput();
    nLcp = topology->getEffectiveSizeOutput();
    if (nLcp != 0)
      assembleM();
  }

  if (nLcp != 0)
  {
    computeQ(time);
    // check z and w sizes and reset if necessary
    if (z == NULL)
    {
      z = new SimpleVector(nLcp);
      isZAllocatedIn = true;
    }
    else if (z->size() != nLcp)
    {
      // reset z if it has a wrong size
      if (isZAllocatedIn) delete z;
      z = new SimpleVector(nLcp);
      isZAllocatedIn = true;
    }

    if (w == NULL)
    {
      w = new SimpleVector(nLcp);
      isWAllocatedIn = true;
    }
    else if (w->size() != nLcp)
    {
      // reset w if it has a wrong size
      if (isWAllocatedIn) delete w;
      w = new SimpleVector(nLcp);
      isWAllocatedIn = true;
    }
    w->zero();
    z->zero();
  }
  OUT("LCP::preLCP()\n");
}

void LCP::assembleM() //
{

  if (M == NULL)
  {
    M = new SimpleMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }
  else if (M->size(0) != nLcp || M->size(1) != nLcp)
  {
    // reset M matrix if it has a wrong size
    if (isMAllocatedIn) delete M;
    M = new SimpleMatrix(nLcp, nLcp);
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
      SiconosMatrix * reducedBlock = new SimpleMatrix(size, size);
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
          reducedCoupledBlock = new SimpleMatrix(size, sizeLinked);
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

void LCP::computeQ(const double& time)
{
  IN("LCP::computeQ(void)\n");
  if (q == NULL)
  {
    q = new SimpleVector(nLcp);
    isQAllocatedIn = true;
  }
  else if (q->size() != nLcp)
  {
    // reset q if it has a wrong size
    if (isQAllocatedIn) delete q;
    q = new SimpleVector(nLcp);
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
        RuntimeException::selfThrow("LCP::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else if (relationType == LAGRANGIANLINEARRELATION || relationType == LAGRANGIANRELATION)
    {
      LagrangianLinearR *LLR = static_cast<LagrangianLinearR*>(R);
      if (nslaw->getType() == NEWTONIMPACTNSLAW)
      {
        vector<unsigned int> indexMax = topology->getIndexMax(currentInteraction);

        NewtonImpactLawNSL * newton = static_cast<NewtonImpactLawNSL*>(nslaw);
        double e = newton->getE();
        LLR->computeFreeOutput(time);
        if (topology->isTimeInvariant())
        {
          yFree = new SimpleVector(currentInteraction -> getY(1)); // copy, no pointer equality
          *yFree += e * currentInteraction -> getYOld(1);
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
          *yFree += e * *tmp;
          delete tmp;
          for (unsigned int i = 0; i < effectiveSize; i++)
            (*q)(i + pos) = (*yFree)(i);
        }
      }
      else
        RuntimeException::selfThrow("LCP::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else
      RuntimeException::selfThrow("LCP::computeQ not yet implemented for relation of type " + R->getType());

    delete yFree;
  }
  OUT("LCP::computeQ(void)\n");
}

void LCP::compute(const double& time)
{
  // --- Prepare data for LCP computing ---
  preLCP(time);

  // --- Call Numerics solver ---
  if (nLcp != 0)
  {
    int info;
    int Nlcp = (int)nLcp;
    method solvingMethod = *(solver->getSolvingMethodPtr());
    info = lcp_solver(M->getArray(), q->getArray(), &Nlcp, &solvingMethod, z->getArray(), w->getArray());

    // \warning : info value and signification depends on solver type ...
    check_solver(info);
    // --- Recover the desired variables from LCP output ---
    postLCP(*w, *z);
  }
}

void LCP::postLCP(const SimpleVector& w, const SimpleVector &z)
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
    // Get the 'position' of vector corresponding to the current interaction, in the LCP output:
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


void LCP::display() const
{
  cout << "======= LCP display ======" << endl;
  cout << "| nLcp : " << nLcp << endl;
  cout << "| LCP Matrix M  : " << endl;
  if (M != NULL) M->display();
  else cout << "-> NULL" << endl;
  cout << "| LCP vector q : " << endl;
  if (q != NULL) q->display();
  else cout << "-> NULL" << endl;
  cout << "==========================" << endl;

}

void LCP::saveNSProblemToXML()
{
  IN("LCP::saveNSProblemToXML\n");
  OneStepNSProblem::saveNSProblemToXML();
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
    (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("LCP::saveNSProblemToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveNSProblemToXML\n");
}

void LCP::saveMToXML()
{
  IN("LCP::saveMToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setM(*M);
  }
  else RuntimeException::selfThrow("LCP::saveMToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveMToXML\n");
}

void LCP::saveQToXML()
{
  IN("LCP::saveQToXML\n");
  if (onestepnspbxml != NULL)
  {
    (static_cast<LCPXML*>(onestepnspbxml))->setQ(*q);
  }
  else RuntimeException::selfThrow("LCP::saveQToXML - OneStepNSProblemXML object not exists");
  OUT("LCP::saveQToXML\n");
}

LCP* LCP::convert(OneStepNSProblem* osnsp)
{
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}



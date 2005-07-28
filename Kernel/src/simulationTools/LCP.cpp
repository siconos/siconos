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
  nspbType = LCP_OSNSP;
}

// xml constructor
LCP::LCP(OneStepNSProblemXML* osNsPbXml, Strategy* newStrat):
  OneStepNSProblem(osNsPbXml, newStrat), nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = LCP_OSNSP;
  if (osNsPbXml != NULL)
  {
    LCPXML * xmllcp = (static_cast<LCPXML*>(osNsPbXml));

    // If both q and M are given in xml file, check if sizes are consistent
    if (xmllcp->hasQ() && xmllcp->hasM() && ((xmllcp->getM()).size(0) != (xmllcp->getQ()).size()))
      RuntimeException::selfThrow("LCP: xml constructor, inconsistent sizes between given q and M");

    // first nlcp is given by M matrix size in xml file
    if (xmllcp->hasM())
    {
      nLcp = (xmllcp->getM()).size(0);
      M = new SiconosMatrix(xmllcp->getM());
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
LCP::LCP(Strategy* newStrat, const string& newSolver, const string& newSolvingMethod,
         const int& MaxIter, const double & Tolerance, const string & NormType,
         const double & SearchDirection):
  OneStepNSProblem(newStrat, newSolver, newSolvingMethod, MaxIter, Tolerance, NormType, SearchDirection),
  nLcp(0), w(NULL), z(NULL), M(NULL), q(NULL),
  isWAllocatedIn(false), isZAllocatedIn(false), isMAllocatedIn(false), isQAllocatedIn(false)
{
  nspbType = LCP_OSNSP;
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
    M = new SiconosMatrix(nLcp, nLcp);
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
  // update topology if necessary
  OneStepNSProblem::initialize();

  // get topology to set nLcp
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  nLcp = topology->getEffectiveSizeOutput();

  computeAllBlocks();
  // if relative degree is equal to 0 or 1
  if (topology->isTimeInvariant())
    assembleM();
}

void LCP::computeAllBlocks()
{
  // --- get topology ---
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  // --- get linked-Interactions map
  map< Interaction*, vector<InteractionLink*> > linkedInteractionMap = topology->getLinkedInteractionMap();
  map< Interaction*, vector<InteractionLink*> >::const_iterator itMapInter;
  // --- get interactions list ---
  vector<Interaction*> listInteractions = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
  vector<Interaction*>::iterator iter;
  // --- get time step ---
  double h = strategy->getTimeDiscretisationPtr()->getH(); // time step

  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // the current interaction and its size
    Interaction *currentInteraction = *iter;
    unsigned int sizeInteraction = currentInteraction->getNInteraction();

    // --- DIAGONAL BLOCKS MANAGEMENT ---

    // matrix that corresponds to diagonal block for the current interaction
    SiconosMatrix * currentMatrixBlock = new SiconosMatrix(sizeInteraction, sizeInteraction);
    diagonalBlocksMap[ currentInteraction ] = currentMatrixBlock;

    // get DS list of the current interaction
    vector<DynamicalSystem*> vDS = currentInteraction ->getDynamicalSystems();
    unsigned int nbDS = vDS.size(); // number of DS

    // sum of all DS sizes
    unsigned int globalDSSize = 0;
    vector<DynamicalSystem*>::iterator itDS;
    for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
      globalDSSize += (*itDS)->getN();

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

    //  !!!!!!!! At the time, we only consider LinearDS-LinearTIR  !!!!!!!
    // other cases have not been tested

    // get the relation of the current interaction
    Relation * R = currentInteraction -> getRelationPtr();

    if (R->getType() == LINEARTIRELATION)
    {
      // convert to linearTIR
      LinearTIR *LTIR = static_cast<LinearTIR*>(R);
      // For 1 relation, we need:
      // - the corresponding row of D
      // - a block-row of C, that corresponds to a specific DS
      // - a block of B, that corresponds to a specific DS
      SiconosMatrix* D = LTIR->getDPtr();
      SiconosMatrix * C, *B;

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
          *Drow = D->getRow(i);
          *currentLine = *Drow;
          // loop over the DS of the current interaction
          for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
          {
            unsigned int sizeDS = (*itDS)->getN();
            // get blocks corresponding to the current DS
            C = new SiconosMatrix(sizeInteraction, sizeDS);
            B = new SiconosMatrix(sizeDS, sizeInteraction);
            LTIR->getCBlockDSPtr(*itDS, *C);
            LTIR->getBBlockDSPtr(*itDS, *B);
            // get row i of C
            Crow = new SimpleVector(sizeDS);
            *Crow = C->getRow(i);

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
          unsigned int sizeDS = (*itDS)->getN();
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
    else RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for relation of type " + R->getType());

    // --- EXTRA-DIAGONAL BLOCKS MANAGEMENT ---

    // check if there are linked interactions with current one, and if so get them.
    itMapInter = linkedInteractionMap.find(currentInteraction);
    if (itMapInter != linkedInteractionMap.end())
    {
      vector<InteractionLink*> ILV = linkedInteractionMap[currentInteraction];
      vector<InteractionLink*>::iterator itLink;

      map<Interaction*, SiconosMatrix*> tmpMap;
      SiconosMatrix * coupledInteractionsBlock;

      // loop over LinkedInteractions
      for (itLink = ILV.begin(); itLink != ILV.end(); itLink++)
      {
        Interaction * linkedInteraction = (*itLink)->getLinkedInteractionPtr();
        unsigned int linkedInteractionSize = linkedInteraction->getNInteraction();

        coupledInteractionsBlock = new SiconosMatrix(sizeInteraction, linkedInteractionSize);
        tmpMap[linkedInteraction] = coupledInteractionsBlock;

        // get the list of common DS and their number
        vector<DynamicalSystem*> commonDS = (*itLink)->getCommonDS();
        nbDS = commonDS.size();

        // get the relation of the current interaction
        Relation * RCurrent = currentInteraction -> getRelationPtr();
        Relation * RLinked = linkedInteraction -> getRelationPtr();

        if (RCurrent->getType() == LINEARTIRELATION && RLinked->getType() == LINEARTIRELATION)
        {
          // convert to linearTIR
          LinearTIR *LTIR1 = static_cast<LinearTIR*>(RCurrent);
          LinearTIR *LTIR2 = static_cast<LinearTIR*>(RLinked);

          SiconosMatrix *C, *B;
          coupledInteractionsBlock->zero();
          bool isHomogeneous = false; // \todo to be defined with relative degrees
          if (!isHomogeneous)
          {
            // the resulting row-block
            SimpleVector * currentLine = new SimpleVector(sizeInteraction);
            // a row of C corresponding to a specific DS (its size depends on the DS size)
            SimpleVector * Crow ;

            // compute currentLine =  h sum(j) Crow,j Wj  Bj, with j the list of common DS
            // C belongs to current Interaction and B to linked interaction
            // loop over the relations of the current interaction
            for (unsigned int i = 0; i < sizeInteraction; i++)
            {
              // loop over common DS
              for (itDS = commonDS.begin(); itDS != commonDS.end(); itDS++)
              {
                unsigned int sizeDS = (*itDS)->getN();
                // get blocks corresponding to the current DS
                C = new SiconosMatrix(sizeInteraction, sizeDS);
                B = new SiconosMatrix(sizeDS, linkedInteractionSize);
                LTIR1->getCBlockDSPtr(*itDS, *C);
                LTIR2->getBBlockDSPtr(*itDS, *B);
                // get row i of C
                Crow = new SimpleVector(sizeDS);
                *Crow = C->getRow(i);
                // compute currentLine
                *currentLine +=  h* *Crow * (*W[*itDS]* *B);
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
            for (itDS = vDS.begin(); itDS != vDS.end(); itDS++)
            {
              unsigned int sizeDS = (*itDS)->getN();
              // get blocks corresponding to the current DS
              C = new SiconosMatrix(sizeInteraction, sizeDS);
              B = new SiconosMatrix(sizeDS, linkedInteractionSize);
              LTIR1->getCBlockDSPtr(*itDS, *C);
              LTIR2->getBBlockDSPtr(*itDS, *B);
              *coupledInteractionsBlock += h * *C * (*W[*itDS] * *B);
              delete C;
              delete B;
            }
          }
        }
        else RuntimeException::selfThrow("LCP::computeAllBlocks not yet implemented for relation of type " + R->getType());
      } // end of loop over linked interactions

      // put tmpMap, ie blocks corresponding to interactions linked with current interaction into the map
      extraDiagonalBlocksMap[currentInteraction] = tmpMap;
    }
  } // end of loop over interaction -> increment current interaction
}

void LCP::preLCP(const double& time)
{
  IN("LCP::preLCP()\n");
  // compute M and q operators for LCP problem

  // if relative degree is not equal to 0 or 1
  // get topology
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();
  if (! topology->isTimeInvariant())
  {
    nLcp = topology->getEffectiveSizeOutput();
    // \todo a step to update "effective" topology
    assembleM();
  }

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

  OUT("LCP::preLCP()\n");
}

void LCP::assembleM()
{

  if (M == NULL)
  {
    M = new SiconosMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }
  else if (M->size(0) != nLcp || M->size(1) != nLcp)
  {
    // reset M matrix if it has a wrong size
    if (isMAllocatedIn) delete M;
    M = new SiconosMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }
  M->zero();

  // get topology
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  map< Interaction* , SiconosMatrix*>::iterator itDiago;
  map< Interaction*, unsigned int>  interactionPositionMap =  topology->getInteractionPositionMap();

  unsigned int pos;
  unsigned int col;
  // --- loop over all the interactions ---
  for (itDiago = diagonalBlocksMap.begin(); itDiago != diagonalBlocksMap.end(); itDiago++)
  {
    SiconosMatrix * blockDiago = itDiago->second;
    Interaction * currentInteraction = itDiago->first;

    pos = interactionPositionMap[currentInteraction];

    // diagonal blocks
    M->blockMatrixCopy(*blockDiago, pos, pos); // \todo avoid copy

    // extra-diagonal blocks
    SiconosMatrix * coupledBlock;
    map< Interaction* , map<Interaction *, SiconosMatrix*> >::iterator itCoupled  ;
    itCoupled = extraDiagonalBlocksMap.find(currentInteraction);
    if (itCoupled != extraDiagonalBlocksMap.end())
    {
      map<Interaction *, SiconosMatrix*>::iterator it;
      for (it = extraDiagonalBlocksMap[currentInteraction].begin(); it != extraDiagonalBlocksMap[currentInteraction].end(); it++)
      {
        coupledBlock = (*it).second;
        col = interactionPositionMap[(*it).first ];
        M->blockMatrixCopy(*coupledBlock, pos, col); // \todo avoid copy
      }
    }
  }
}

void LCP::compute(const double& time)
{
  IN("LCP::compute(void)\n");
  int res;

  // --- Prepare data for LCP computing ---
  preLCP(time);

  // --- Call Numerics solver ---
  // !!! Warning: solve_lcp solve Mz - w =q !!!
  if (nLcp != 0)
    res = solve_lcp(M->getArray(), q->getArray(), (int*)&nLcp, &solvingMethod, z->getArray(), w->getArray());

  // --- Recover the desired variables from LCP output ---
  postLCP(*w, *z);

  OUT("LCP::compute(void)\n");
}

void LCP::postLCP(const SimpleVector& w, const SimpleVector &z)
{
  SimpleVector *yDot, *lambda, *y;

  // --- get topology ---
  Topology * topology = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getTopologyPtr();

  // --- get interactions list ---
  vector<Interaction*> listInteractions = strategy->getModelPtr()->getNonSmoothDynamicalSystemPtr()->getInteractions();
  vector<Interaction*>::iterator iter;
  // get Interaction position map
  map< Interaction* , SiconosMatrix*>::iterator itDiago;
  map< Interaction*, unsigned int>  interactionPositionMap =  topology->getInteractionPositionMap();
  unsigned int pos;

  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // the current interaction and its size
    Interaction *currentInteraction = *iter;
    unsigned int sizeInteraction = currentInteraction->getNInteraction();
    // position of current yFree in global (including all interactions) y vector
    pos = interactionPositionMap[currentInteraction];

    yDot = currentInteraction->getYPtr(1);
    y = currentInteraction -> getYPtr(0);
    lambda = currentInteraction ->getLambdaPtr();
    Relation *R = currentInteraction->getRelationPtr();
    string relationType = R->getType();
    for (unsigned int i = 0; i < sizeInteraction; i++)
    {
      if (relationType == LAGRANGIANLINEARRELATION)
        (*yDot)(i) = w(i + pos); // \todo compute y
      else if (relationType == LINEARTIRELATION)
        (*y)(i) = w(i + pos);
      else
        RuntimeException::selfThrow("LCP::compute not yet implemented for relation of type " + R->getType());

      (*lambda)(i) = z(i + pos);
    }
  } // end of loop over interactions
}

void LCP::computeM()
{
  IN("LCP::computeM(void)\n");
  int number, orgDSRank, connectedDSRank;
  int currentActiveInteraction = 0;
  int interConnectedNumber = 0;
  vector<Connection*> vCo;
  SiconosMatrix *H, *WW, *Mtmp ;
  bool isWWAllocatedIn = false;
  SiconosMatrix orgH, connectedH, wTmp;
  Relation *R, *RConnected;
  LagrangianLinearR *LLR ;
  unsigned int i;

  unsigned int currentPosition = 0;

  if (M == NULL)
  {
    M = new SiconosMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }
  else if (M->size(0) != nLcp || M->size(1) != nLcp)
  {
    // reset M matrix if it has a wrong size
    if (isMAllocatedIn) delete M;
    M = new SiconosMatrix(nLcp, nLcp);
    isMAllocatedIn = true;
  }



  // \todo improve WW organization (WW block-structured)

  // --- For each interaction in the Map (ie active interaction) ... ---
  map<Interaction* , vector <Connection*> >::iterator iter;
  for (iter = connectedInteractionMap.begin(); iter != connectedInteractionMap.end(); ++iter)
  {
    Interaction *currentInteraction = iter->first;

    // At the time, W is inversed in initialize of Moreau.cpp
    // -> \todo improve this step using forward-backward to compute inverse of Mlcp blocks rather than W
    vector<DynamicalSystem*> vDS = currentInteraction ->getDynamicalSystems();
    unsigned int sizeDS = vDS.size();

    // Get W matrix of each DS concerned by the interaction and save in into v
    vector<SiconosMatrix*> v;
    vector<OneStepIntegrator*> OsiV;
    OsiV.reserve(sizeDS);
    for (i = 0; i < sizeDS; i++)
    {
      number = vDS[i]->getNumber();
      OsiV[i] = strategy->getIntegratorOfDSPtr(number);
      if (OsiV[i]->getType() == MOREAU_INTEGRATOR)
        v.push_back((static_cast<Moreau*>(OsiV[i]))->getWPtr());
      else
        RuntimeException::selfThrow("LCP::computeM not yet implemented for Integrator of type " + OsiV[i]->getType());
    }

    // Built block matrix WW with all the W previously saved in v..
    // \todo do not assemble WW but work directly with blocks
    unsigned int size = v[0]->size(0);
    for (i = 1; i < sizeDS; i++)
      size += v[i]->size(0);
    WW = new SiconosMatrix(size, size);
    isWWAllocatedIn = true;
    *WW = BlockMatrixAssemble(v);

    // --- Get the relation parameters and compute M ---
    R = currentInteraction->getRelationPtr();
    if (R->getType() == LAGRANGIANLINEARRELATION)
    {
      //  compute H W Ht
      LLR = static_cast<LagrangianLinearR*>(R);
      H = LLR->getHPtr();
      unsigned int sizeH = H->size(0);
      Mtmp = new SiconosMatrix(sizeH, sizeH);
      *Mtmp = *H * WW->multTranspose(*H);
    }
    else if (R->getType() == LINEARTIRELATION)
    {
      LinearTIR *LTIR = static_cast<LinearTIR*>(R);
      SiconosMatrix* Bloc = LTIR->getBPtr();
      SiconosMatrix* Cloc = LTIR->getCPtr();
      SiconosMatrix* Dloc = LTIR->getDPtr();
      double h = strategy->getTimeDiscretisationPtr()->getH(); // time step
      unsigned int sizeMtmp = Dloc->size(0);
      Mtmp = new SiconosMatrix(sizeMtmp, sizeMtmp);
      *Mtmp = (h * (*Cloc * *WW * *Bloc)) + *Dloc;
    }
    else RuntimeException::selfThrow("LCP::computeM [level1] not yet implemented for relation of type " + R->getType());

    // M assembly
    M->blockMatrixCopy(*Mtmp, currentPosition, currentPosition);
    delete Mtmp;

    // --- Compute M for connected interactions ---
    interConnectedNumber = 0;
    if (iter ->second[0] != NULL)
    {
      // get from the map the connexion vector of the current interaction
      vCo = iter -> second ;
      for (unsigned int k = 0; k < vCo.size(); k++)
      {
        orgDSRank = vCo[k]->originInteractionDSRank;
        connectedDSRank = vCo[k]->connectedInteractionDSRank;

        // get W(Moreau) of common DS
        wTmp = *v[orgDSRank];

        // get H matrix of the common DS
        // /!\ we supose that all the DS have the same size !!!!
        if (R->getType() == LAGRANGIANLINEARRELATION)
        {
          LLR = static_cast<LagrangianLinearR*>(R);
          // /!\ copy of matrices !!!
          orgH = LLR->getHRelatingToDS(orgDSRank);
        }
        else
          RuntimeException::selfThrow("LCP::computeM [level2] not yet implemented for relation of type " + R->getType());

        // get H matrix of the connected DS
        RConnected = vCo[k]->connected->getRelationPtr();
        if (RConnected->getType() == LAGRANGIANLINEARRELATION)
        {
          LLR = static_cast<LagrangianLinearR*>(RConnected);
          // /!\ copy of matrices !!!
          connectedH = LLR->getHRelatingToDS(connectedDSRank);
        }
        else
          RuntimeException::selfThrow("LCP::computeM [level3] not yet implemented for relation of type " + RConnected->getType());
        unsigned int sizeOrgH = orgH.size(0);
        unsigned int sizeConnectedH = connectedH.size(0);
        Mtmp = new SiconosMatrix(sizeOrgH, sizeConnectedH);
        *Mtmp = orgH * wTmp.multTranspose(connectedH);

        if (interConnectedNumber == currentActiveInteraction)
          interConnectedNumber++;
        //        M->blockMatrixCopy(Mtmp, currentPosition, connectedPosition);
        delete Mtmp;
        interConnectedNumber++;
      }
    }
    // incrementation of the number of active interaction
    currentActiveInteraction++;
    currentPosition += currentInteraction->getNInteraction();
  }
  if (isWWAllocatedIn) delete WW;
  OUT("LCP::computeM(void)\n");
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
  map< Interaction*, unsigned int>  interactionPositionMap =  topology->getInteractionPositionMap();
  unsigned int pos;
  Relation *R;
  NonSmoothLaw *nslaw;

  // --- loop over all the interactions ---
  for (iter = listInteractions.begin(); iter != listInteractions.end(); ++iter)
  {
    // get current interaction, its relation and its nslaw
    Interaction *currentInteraction = *iter;
    unsigned int sizeInteraction = currentInteraction->getNInteraction();
    R = currentInteraction->getRelationPtr();
    nslaw = currentInteraction->getNonSmoothLawPtr();
    SimpleVector * yFree;
    // position of current yFree in global (including all interactions) y vector
    pos = interactionPositionMap[currentInteraction];

    if (R->getType() == LINEARTIRELATION)
    {
      LinearTIR *LTIR = static_cast<LinearTIR*>(R);
      string nslawType = nslaw->getType() ;
      if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
      {
        LTIR->computeFreeOutput(); // free output is saved in y
        yFree = new SimpleVector(currentInteraction->getY(0));   // copy, no pointer equality
        for (unsigned int i = 0; i < sizeInteraction; i++)
          (*q)(i + pos) = -(*yFree)(i); // - because of LCP solver, Mz - w =q
      }
      else
        RuntimeException::selfThrow("LCP::computeQ not yet implemented for NSlaw of type " + nslaw->getType() + "and for relation of type " + R->getType());
    }
    else if (R->getType() == LAGRANGIANLINEARRELATION)
    {
      LagrangianLinearR *LLR = static_cast<LagrangianLinearR*>(R);
      if (nslaw->getType() == NEWTONIMPACTLAWNSLAW)
      {
        NewtonImpactLawNSL * newton = static_cast<NewtonImpactLawNSL*>(nslaw);
        double e = newton->getE();
        LLR->computeFreeOutput(time);
        yFree = new SimpleVector(currentInteraction -> getY(1)); // copy, no pointer equality
        *yFree += e * currentInteraction -> getYOld(1);
        for (unsigned int i = 0; i < sizeInteraction; i++)
          (*q)(i + pos) = -(*yFree)(0);
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
  cout << "LCP::convert (OneStepNSProblem* osnsp)" << endl;
  LCP* lcp = dynamic_cast<LCP*>(osnsp);
  return lcp;
}



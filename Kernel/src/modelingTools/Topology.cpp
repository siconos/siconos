#include "Topology.h"

using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

// default
Topology::Topology():
  sizeOutput(0), effectiveSizeOutput(sizeOutput), isTopologyUpToDate(false), isTopologyTimeInvariant(true), nsds(NULL)
{}

// with NonSmoothDynamicalSystem
Topology::Topology(NonSmoothDynamicalSystem* newNsds):
  sizeOutput(0), effectiveSizeOutput(sizeOutput), isTopologyUpToDate(false), isTopologyTimeInvariant(true), nsds(newNsds)
{
  updateTopology();
}

// destructor
Topology::~Topology()
{
  map< Interaction*, std::vector<InteractionLink*> >::iterator mapIt;
  // for all elements in the map ...
  for (mapIt = linkedInteractionMap.begin(); mapIt != linkedInteractionMap.end(); mapIt++)
  {
    // get vector of interactionLink
    vector<InteractionLink*> linkVector = (*mapIt).second;
    vector<InteractionLink*>::iterator linkVectorIt;
    // delete each interactionLink
    for (linkVectorIt = linkVector.begin(); linkVectorIt != linkVector.end(); linkVectorIt++)
    {
      delete(*linkVectorIt);
      *linkVectorIt = NULL ;
    }
    linkVector.clear();
  }
  linkedInteractionMap.clear();
  nsds = NULL;
}

void Topology::updateTopology()
{
  //-- Get all the interactions --
  vector<Interaction*> listInteractions = nsds->getInteractions();
  // -- Compute sizeOutput --
  vector<Interaction*>::iterator it;
  // loop through interactions list
  unsigned int size;
  sizeOutput = 0;
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
  {
    size = (*it)->getNInteraction();
    sizeOutput += size;
    indexMaxMap[*it].resize(size, 0); // by default, all indexMax are set to 0.
  }

  // fill linkedInteractionMap and interactionPositionMap in:
  computeLinkedInteractionMap();
  computeInteractionPositionMap();

  // fill RelativeDegreesMaps in
  computeRelativeDegreesMap();

  // fill IndexMinMap
  computeIndexMinMap();

  if (isTopologyTimeInvariant)
    effectiveSizeOutput = sizeOutput;
  // else this value will be computed in OneStepNSProblem
  // and so for indexMap.

  isTopologyUpToDate = true;

}

void Topology::computeEffectiveSizeOutput()
{
  effectiveSizeOutput = 0;

  // Get all the interactions
  vector<Interaction*> listInteractions = nsds->getInteractions();

  vector<Interaction*>::iterator it;
  // loop over interactions list
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
  {
    // loop over ouput vector
    vector<unsigned int>::iterator it2;
    for (it2 = indexMaxMap[*it].begin(); it2 != indexMaxMap[*it].end(); it2++)
      effectiveSizeOutput += indexMaxMap[*it][*it2] - indexMinMap[*it][*it2] + 1 ;
  }
}

// Linked interactions map computing:
void Topology::computeLinkedInteractionMap()
{
  vector<DynamicalSystem*> dsOrig, dsLinked;

  // Get all the interactions
  vector<Interaction*> listInteractions = nsds->getInteractions();

  linkedInteractionMap.clear();

  vector<Interaction*>::iterator itOrig;

  // -- loop over all the interactions of the NonSmoothDynamicalSystem --
  for (itOrig = listInteractions.begin(); itOrig != listInteractions.end(); itOrig++)
  {
    dsOrig = (*itOrig)->getDynamicalSystems();

    vector<Interaction*>::iterator itLinked;
    // -- check all other interactions --
    for (itLinked = listInteractions.begin(); itLinked != listInteractions.end() && itLinked != itOrig; itLinked++)
    {
      dsLinked = (*itLinked)->getDynamicalSystems();

      vector<DynamicalSystem*> commonDS;  // list of ds common to both interactions
      vector<DynamicalSystem*>::iterator itDS;

      // compare list of DS of the 2 interactions
      for (unsigned int k = 0; k < dsLinked.size(); k++)
      {
        itDS = find(dsOrig.begin(), dsOrig.end(), dsLinked[k]);
        if (itDS != dsOrig.end()) commonDS.push_back(*itDS);
      }

      // built linkedInteractionMap
      if (commonDS.size() != 0)
        linkedInteractionMap[(*itOrig) ].push_back(new InteractionLink((*itOrig), (*itLinked), commonDS));
    }
  }
}

// Linked interactions map computing:
void Topology::computeInteractionPositionMap()
{
  interactionPositionMap.clear();

  // Get all the interactions
  vector<Interaction*> listInteractions = nsds->getInteractions();

  vector<Interaction*>::iterator it;
  unsigned int currentPosition = 0;
  // loop through interactions list
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
  {
    interactionPositionMap[*it] = currentPosition;
    currentPosition += (*it)->getNInteraction();
  }
}

// Compute relative degrees map
void Topology::computeRelativeDegreesMap()
{
  // for each interaction relative degree vector depends on NonSmooth Law and relations

  relativeDegreesMap.clear();

  // Get all the interactions
  vector<Interaction*> listInteractions = nsds->getInteractions();

  // loop through interactions list
  vector<Interaction*>::iterator it;
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
    relativeDegreesMap[(*it) ] = computeRelativeDegrees(*it) ;
}

// Compute relative degree of a single interaction
vector<unsigned int> Topology::computeRelativeDegrees(Interaction * inter)
{
  // get NonSmooth law, relation and their types
  NonSmoothLaw * nslaw = inter ->getNonSmoothLawPtr();
  string nslawType = nslaw -> getType();
  Relation * relation  = inter ->getRelationPtr();
  string relationType = relation->getType();

  vector<unsigned int> relativeDegree;
  unsigned int sizeInter = inter->getNInteraction();

  // loop over various non smooth law types
  // \todo check isTimeInvariant correctly
  if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
  {
    relativeDegree.resize(sizeInter, 0);
    // \todo compute correctly
  }
  else if (nslawType == NEWTONIMPACTLAWNSLAW)
  {
    relativeDegree.resize(sizeInter, 2);
    isTopologyTimeInvariant = false;
  }
  else
    RuntimeException::selfThrow("Topology::computeRelativeDegrees, not yet implemented for non smooth law of type" + nslawType);

  return relativeDegree;
}

// Compute indexMin map
void Topology::computeIndexMinMap()
{
  // for each interaction indexMin depends on NonSmooth Law

  indexMinMap.clear();

  // Get all the interactions
  vector<Interaction*> listInteractions = nsds->getInteractions();

  // loop through interactions list
  vector<Interaction*>::iterator it;
  for (it = listInteractions.begin(); it != listInteractions.end(); it++)
    indexMinMap[(*it) ] = computeIndexMin(*it);

}

// Compute relative degree of a single interaction
vector<unsigned int> Topology::computeIndexMin(Interaction * inter)
{
  // get NonSmooth law, relation and their types
  NonSmoothLaw * nslaw = inter ->getNonSmoothLawPtr();
  string nslawType = nslaw -> getType();
  Relation * relation  = inter ->getRelationPtr();
  string relationType = relation->getType();

  vector<unsigned int> indexMin;
  unsigned int sizeInter = inter->getNInteraction();

  // loop over various non smooth law types
  if (nslawType == COMPLEMENTARITYCONDITIONNSLAW)
  {
    indexMin.resize(sizeInter, 0);
  }
  else if (nslawType == NEWTONIMPACTLAWNSLAW)
    indexMin.resize(sizeInter, 1);

  else
    RuntimeException::selfThrow("Topology::computeIndexMin, not yet implemented for non smooth law of type" + nslawType);

  return indexMin;
}


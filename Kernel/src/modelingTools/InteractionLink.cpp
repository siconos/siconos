#include "InteractionLink.h"

using namespace std;

// --- CONSTRUCTORS/DESTRUCTOR ---

InteractionLink::InteractionLink(Interaction* itOrig, Interaction* itLinked, const std::vector<DynamicalSystem*>& dsList):
  originInteraction(itOrig), linkedInteraction(itLinked), commonDS(dsList)
{}

InteractionLink::~InteractionLink()
{
  originInteraction = NULL;
  linkedInteraction = NULL;
  commonDS.resize(1, NULL);
}

void InteractionLink::display() const
{
  cout << " ===== interactionLink display ===== " << endl;
  cout << " Origin interaction number: " << endl;
  if (originInteraction != NULL)
    cout << originInteraction->getNumber();
  else
    cout << "-> NULL" << endl;
  cout << " linked interaction number: " << endl;
  if (linkedInteraction != NULL)
    cout << linkedInteraction->getNumber();
  else
    cout << "-> NULL" << endl;

  cout << " common dynamical systems: " << endl;
  for (unsigned int i = 0; i < commonDS.size(); i++)
    commonDS[i]->display();
  cout << " =================================== " << endl;
}

// default (private) constructor
InteractionLink::InteractionLink(): originInteraction(NULL), linkedInteraction(NULL)
{}


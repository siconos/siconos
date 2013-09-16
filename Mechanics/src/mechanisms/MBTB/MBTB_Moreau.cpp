#include "MBTB_Moreau.hpp"
#define DEBUG_STDOUT
//#define DEBUG_MESSAGES
//#define DEBUG_WHERE_MESSAGES
#include <debug.h>
#include "CxxStd.hpp"

//#define STANDARD_ACTIVATION
#define FIRSTWAY_ACTIVATION


MBTB_Moreau::MBTB_Moreau(SP::DynamicalSystem ds , double d):Moreau(ds ,d)
{
  _deactivateYPosThreshold= 1e-4;
  _deactivateYVelThreshold=0;
  _activateYPosThreshold=0;
  _activateYVelThreshold=100; 
}

#ifdef STANDARD_ACTIVATION
bool MBTB_Moreau::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  

  assert(i==1);
  double h = simulationLink->timeStep();
  double y = (inter->y(i-1))->getValxxue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity 
  
  double gamma = 1.0/2.0;
  if (_useGamma)
  {
    gamma = _gamma;
  } 
  DEBUG_PRINTF("MBTB_Moreau::addInteractionInIndexSet yref=%e, yDot=%e, y_estimated=%e.\n", y, yDot, y+gamma*h*yDot);
  y += gamma*h*yDot;
  assert(!isnan(y));
  if(y<=0)
    DEBUG_PRINT("MBTB_Moreau::addInteractionInIndexSet ACTIVATE.\n");
  return (y<=0);
}
bool MBTB_Moreau::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  assert(i==1);
  double h = simulationLink->timeStep();
  double y = (inter->y(i-1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity 
  double gamma = 1.0/2.0;
  if (_useGamma)
  {
    gamma = _gamma;
  } 
  DEBUG_PRINTF("MBTB_Moreau::addInteractionInIndexSet yref=%e, yDot=%e .\n", y, yDot, y+gamma*h*yDot);
  y += gamma*h*yDot;
  assert(!isnan(y));
  if(y>0)
    DEBUG_PRINT("MBTB_Moreau::removeInteractionInIndexSet DEACTIVATE.\n");
  return (y>0);
}
#endif

#ifdef FIRSTWAY_ACTIVATION
bool MBTB_Moreau::addInteractionInIndexSet(SP::Interaction inter, unsigned int i)
{
  
  assert(i==1);
  double y = (inter->y(i-1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity 
 
  DEBUG_PRINTF("MBTB_Moreau::addInteractionInIndexSet yref=%e, yDot=%e.\n", y, yDot);

  DEBUG_PRINTF("MBTB_Moreau::addInteractionInIndexSet  _activateYPosThreshold =%e, _activateYVelThreshold=%e\n",
               _activateYPosThreshold ,
               _activateYVelThreshold );
  
  assert(!isnan(y));
  
  if(y<=_activateYPosThreshold)
    DEBUG_PRINT("MBTB_Moreau::addInteractionInIndexSet ACTIVATE.\n");
  return (y<=_activateYPosThreshold);
}

bool MBTB_Moreau::removeInteractionInIndexSet(SP::Interaction inter, unsigned int i)

{
  assert(i==1);
  double h = simulationLink->timeStep();
  double y = (inter->y(i-1))->getValue(0); // for i=1 y(i-1) is the position
  double yDot = (inter->y(i))->getValue(0); // for i=1 y(i) is the velocity 

  DEBUG_PRINTF("MBTB_Moreau::removeInteractionInIndexSet yref=%e, yDot=%e .\n", y, yDot);

  DEBUG_PRINTF("MBTB_Moreau::removeInteractionInIndexSet  _deactivateYPosThreshold =%e, _deactivateYVelThreshold=%e\n",
               _deactivateYPosThreshold ,
               _deactivateYVelThreshold );

  assert(!isnan(y));
  if(y>_deactivateYPosThreshold && yDot>=_deactivateYVelThreshold)
    DEBUG_PRINT("MBTB_Moreau::removeInteractionInIndexSet DEACTIVATE.\n");
  return (y>_deactivateYPosThreshold && yDot>=_deactivateYVelThreshold);
}

#endif

#ifndef OccSpaceFilter_hpp
#define OccSpaceFilter_hpp

#include "SpaceFilter.hpp"

class OccSpaceFilter : public SpaceFilter
{
public:
  OccSpaceFilter() : SpaceFilter() {};

  virtual void updateInteractions(SP::Simulation) {};

};

#endif

#ifndef OccSpaceFilter_hpp
#define OccSpaceFilter_hpp

#include "SpaceFilter.hpp"

class OccSpaceFilter : public SpaceFilter
{
public:
  OccSpaceFilter(SP::NonSmoothDynamicalSystem nsds) : SpaceFilter(nsds) {};

  virtual void buildInteractions(double time) {};

};

#endif

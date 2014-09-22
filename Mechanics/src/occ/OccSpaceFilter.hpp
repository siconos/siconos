#ifndef OccSpaceFilter_hpp
#define OccSpaceFilter_hpp

#include "SpaceFilter.hpp"

class OccSpaceFilter : public SpaceFilter
{
public:
  OccSpaceFilter(SP::Model model) : SpaceFilter(model) {};

  virtual void buildInteractions(double time) {};

};

#endif

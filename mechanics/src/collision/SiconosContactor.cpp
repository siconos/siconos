
#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"

#include <boost/make_shared.hpp>

SiconosContactor::SiconosContactor(SP::SiconosShape _shape,
                                   SP::SiconosVector _offset,
                                   int _collision_group)
  : shape(_shape), offset(_offset), collision_group(_collision_group)
{
  if (!offset)
  {
    offset = std11::make_shared<SiconosVector>(7);
    offset->zero();
    (*offset)(3) = 1.0;
  }
}


#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"
SiconosContactor::SiconosContactor(SP::SiconosShape _shape,
                                   SP::SiconosVector _offset,
                                   int _collision_group)
  : shape(_shape), offset(_offset), collision_group(_collision_group)
{
  // First strategy: fill offset with the identity if not provided
  // if(!offset)
  // {
  //   offset = std::make_shared<SiconosVector>(7);
  //   offset->zero();
  //   (*offset)(3) = 1.0;
  // }

  // Second strategy: leave offset as a null pointer if identity is provided
  if(offset)
  {
    if ((*offset)(3) == 1.0
        && (*offset)(0) == 0.0
        && (*offset)(1) == 0.0
        && (*offset)(2) == 0.0
        && (*offset)(4) == 0.0
        && (*offset)(5) == 0.0
        && (*offset)(6) == 0.0
      )
    {
      offset.reset();
    }
  }
}


#include "SiconosContactor.hpp"
#include "SiconosShape.hpp"

void SiconosContactor::addShape(SP::SiconosShape shape,
                                SP::SiconosVector offset)
{
  if (!offset) { offset.reset(new SiconosVector(7)); offset->zero(); }
  _shapes.push_back(std::make_pair(shape, offset));
}

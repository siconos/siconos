
#include "Contactor.hpp"

void Contactor::addShape(SP::SiconosShape shape)
{
  _shapes.push_back(shape);
}

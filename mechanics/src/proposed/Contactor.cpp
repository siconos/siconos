
#include "Contactor.hpp"
#include "SiconosShape.hpp"

void Contactor::addShape(SP::SiconosShape shape)
{
  _shapes.push_back(shape);
}

void Contactor::setPosition(const SP::SiconosVector position)
{
  std::vector<SP::SiconosShape>::iterator s;
  for (s = _shapes.begin();
       s != _shapes.end(); s++)
  {
    (*s)->setPosition(position);
  }
}

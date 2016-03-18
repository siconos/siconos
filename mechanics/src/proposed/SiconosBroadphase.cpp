
#include "SiconosBroadphase.hpp"
#include "BodyDS.hpp"
#include "Model.hpp"

void  SiconosBroadphase::visit(SP::BodyDS body)
{
  printf("BodyDS: %p\n", *body);
  //contactor->acceptSP(shared_from_this());
}

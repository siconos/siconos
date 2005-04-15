
















#include "SiconosObject.h"
#include <iostream>


void BaseObjectSiconos::display()
{
  cout << this->i << endl;
  cout << this->type << endl;
}


void ObjectSiconos::display()
{
  //Do we need to call it inside the braces ?
  BaseObjectSiconos::display();
  cout << this->att << endl;
  cout << this->type << endl;
}


int main(void)

{

  ObjectSiconos o1;
  o1.display();



  ObjectSiconos o2(10, "derived_TYPE");
  o2.display();
}

#include "MBTB_Contact.hpp"
#include "MBTB_ContactRelation.hpp"
#include "MBTB_FC3DContactRelation.hpp"

MBTB_Contact::MBTB_Contact() {}

MBTB_Contact::MBTB_Contact(unsigned int id,const std::string& ContactName, unsigned int indexBody1, int indexBody2,unsigned int indexCAD1,unsigned int indexCAD2, int withFriction)
{
  strcpy(_ContactName,ContactName.c_str());
  _id=id;
  _indexBody1=indexBody1;
  _indexBody2=indexBody2;
  _indexCAD1=indexCAD1;
  _indexCAD2=indexCAD2;
  _withFriction=0;
  _en=0.5;
  _et=0.0;
  _curTimeh=-1.0;
  _Offset=0.01;
  _OffsetP1=1;
  _withFriction=withFriction;
  _normalFromFace1=1;
  if(_withFriction)
    _Relation.reset(new MBTB_FC3DContactRelation(this));
  else
    _Relation.reset(new MBTB_ContactRelation(this));
}
 void MBTB_Contact::setInteraction(SP::Interaction newInteraction)
  {

    std::cout << "MBTB_Contact::setInteraction(SP::Interaction newInteraction)"<< std::endl;
    std::cout << "_interaction before"<< _interaction << std::endl;
    std::cout << "newinteraction "<< newInteraction << std::endl;
    newInteraction->display();
    
    _interaction = newInteraction ;
    std::cout << "_interaction after"<< _interaction << std::endl;
    _interaction->display();

  }

#include "MBTB_FC3DContactRelation.hpp"
#include "CADMBTB_API.hpp"
#include "MBTB_DATA.hpp"
#include "ace.h"


//#define DEBUG_MESSAGES
//#define DEBUG_STDOUT
#include "debug.h"

#define OUTPUT_H_IN_FILE

MBTB_FC3DContactRelation::~MBTB_FC3DContactRelation()
{
  ;
}

MBTB_FC3DContactRelation::MBTB_FC3DContactRelation(MBTB_Contact * pC)
{
  _pContact=pC;
}

MBTB_FC3DContactRelation::MBTB_FC3DContactRelation()
{
  _pContact=nullptr;
}
/*This function has to compute the distance between the objects*/
void MBTB_FC3DContactRelation::computeh(double time, const BlockVector& q0, SiconosVector& y)
{
//  DSIterator itDS=_pContact->interaction()->dynamicalSystemsBegin();
//  SP::DynamicalSystem aux = *itDS;
// if(sPrintDist)
//  {
//    printf("MBTB_FC3DContactRelation::computeh Start display for contact name %s\n",_pContact->_ContactName);
//  }
//  if(sDS[_pContact->_indexBody1] != aux)
//  {
//    printf("MBTB_FC3DContactRelation::computeh wrong short of DS\n");
//    exit(1);

//  }

  DEBUG_PRINT("MBTB_FC3DContactRelation::computeh(double time, BlockVector& q0, SiconosVector& y )\n");
  DEBUG_EXPR(_pContact->interaction()->y(0)->display(););
  DEBUG_EXPR(y.display(););
  //SP::SiconosVector y = _pContact->interaction()->y(0);




  //if (_pContact->_curTimeh + 1e-9 < time){
  ACE_times[ACE_TIMER_DIST].start();
  double X1,X2,Y1,Y2,Z1,Z2,n1x,n1y,n1z;
  CADMBTB_getMinDistance(_pContact->_id,_pContact->_indexCAD1,_pContact->_indexCAD2,
                         X1,Y1,Z1,
                         X2,Y2,Z2,
                         n1x,n1y,n1z,_pContact->_normalFromFace1,
                         _pContact->_dist);

  if(sPrintDist)
  {
    printf("    Minimal distance computed from CAD and n2qn1 : %lf \n",_pContact->_dist);
    printf("    Proximal point 1 computed from CAD :  X1=%lf, Y1=%lf, Z1=%lf \n",X1,Y1,Z1);
    printf("    Proximal point 2 computed from CAD :  X2=%lf, Y2=%lf, Z2=%lf \n",X2,Y2,Z2);
    if(_pContact->_normalFromFace1)
      printf("    Normal vector computed from CAD taken from  Object 1 :  nx=%lf, ny=%lf, nz=%lf \n",n1x,n1y,n1z);
    else
      printf("    Normal vector computed from CAD taken from  Object 2 :  nx=%lf, ny=%lf, nz=%lf \n",n1x,n1y,n1z);
  }

  //_Pc1->setValue(0,X1); _Pc1->setValue(1,Y1); _Pc1->setValue(2,Z1);
  if(_pContact->_OffsetP1)
  {
    _Pc1->setValue(0,X1+_pContact->_Offset*n1x);
    _Pc1->setValue(1,Y1+_pContact->_Offset*n1y);
    _Pc1->setValue(2,Z1+_pContact->_Offset*n1z);
    _Pc2->setValue(0,X2);
    _Pc2->setValue(1,Y2);
    _Pc2->setValue(2,Z2);
    if(sPrintDist)
    {
      printf("    OffSet is added from contact point PC1 : newPC1 =  PC1 + Offset.n ");
      printf("    OffSet %lf\n ",_pContact->_Offset);
    }
  }
  else
  {
    _Pc1->setValue(0,X1);
    _Pc1->setValue(1,Y1);
    _Pc1->setValue(2,Z1);
    _Pc2->setValue(0,X2-_pContact->_Offset*n1x);
    _Pc2->setValue(1,Y2-_pContact->_Offset*n1y);
    _Pc2->setValue(2,Z2-_pContact->_Offset*n1z);
    if(sPrintDist)
    {
      printf("    OffSet is substracted from contact point PC2 : newPC1 =  PC2 - Offset.n ");
      printf("    point PC2 : X2=%lf,Y2=%lf,Z2=%lf",X2-_pContact->_Offset*n1x,Y2-_pContact->_Offset*n1y,Z2-_pContact->_Offset*n1z);
      printf("    OffSet %lf\n",_pContact->_Offset);
    }
  }
  /*Because in CAD model, the normal is going outside of the body.*/
  _Nc->setValue(0,-n1x);
  _Nc->setValue(1,-n1y);
  _Nc->setValue(2,-n1z);

  DEBUG_EXPR(_Nc->display(););
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););



  ACE_times[ACE_TIMER_DIST].stop();
  _pContact->_dist-=_pContact->_Offset;
  _pContact->_curTimeh=time;


  if(sPrintDist)
  {
    printf("MBTB_FC3DContactRelation compute h of %s: %14e \n",_pContact->_ContactName,_pContact->_dist);
    printf("MBTB_FC3DContactRelation compute h the normal rentrante : nx=%lf, ny=%lf, nz=%lf \n",-n1x,-n1y,-n1z);
  }

  //}
  // if (_pContact->_curTimeh + 1e-9 < time){
  //   ACE_times[ACE_TIMER_DIST].start();
  //   double X1,X2,Y1,Y2,Z1,Z2,n1x,n1y,n1z;
  //   CADMBTB_getMinDistance(_pContact->_id,_pContact->_indexCAD1,_pContact->_indexCAD2 ,
  // 			   _pContact->_X1,_pContact->_Y1,_pContact->_Z1,
  // 			   _pContact->_X2,_pContact->_Y2,_pContact->_Z2,
  // 			   _pContact->_n1X,_pContact->_n1Y,_pContact->_n1Z,_pContact->_dist);
  //   ACE_times[ACE_TIMER_DIST].stop();
  //   _pContact->_dist-=_pContact->sOffset;
  //   _pContact->_curTimeh=time;
  // }
  y.setValue(0,_pContact->_dist);
  DEBUG_EXPR(y.display(););

  //SP::NewtonEulerR ner =(boost::static_pointer_cast<NewtonEulerR>(interaction()->relation()));
  //ner->yProj()->setValue(0,_pContact->_dist);
  // _Pc1->setValue(0,0.5*(_pContact->_X1+_pContact->_X2));
  // _Pc1->setValue(1,0.5*(_pContact->_Y1+_pContact->_Y2));
  // _Pc1->setValue(2,0.5*(_pContact->_Z1+_pContact->_Z2));
  // _Nc->setValue(0,-_pContact->_n1X);
  // _Nc->setValue(1,-_pContact->_n1Y);
  // _Nc->setValue(2,-_pContact->_n1Z);
  if(sPrintDist)
  {
    printf("MBTB_ContactRelation::computeh End display for contact name %s\n",_pContact->_ContactName);
  }

}

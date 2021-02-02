#include "MBTB_ContactRelation.hpp"
#include "MBTB_DATA.hpp"
#include "CADMBTB_API.hpp"
#include <boost/math/quaternion.hpp>
#include "ace.h"

#define DEBUG_MESSAGES
#include "debug.h"

MBTB_ContactRelation::~MBTB_ContactRelation()
{
  ;
}



MBTB_ContactRelation::MBTB_ContactRelation(): NewtonEuler1DR()
{
  _pContact=nullptr;
}
MBTB_ContactRelation::MBTB_ContactRelation(MBTB_Contact *pC): NewtonEuler1DR()
{
  _pContact=pC;
}

void MBTB_ContactRelation::computeh(double time, BlockVector& q0, SiconosVector& y)
{


  DEBUG_PRINT("MBTB_ContactRelation::computeh(double time, BlockVector& q0, SiconosVector& y)\n");


  printf("sPrintDist=%d\t",sPrintDist);
  if(sPrintDist)
  {
    printf("MBTB_ContactRelation::computeh Start display for contact name %s\n",_pContact->_ContactName);
  }



  //printf("contactName :%s\n",_ContactName);
  //if (_pContact->_curTimeh + 1e-9 < time){
  ACE_times[ACE_TIMER_DIST].start();
  double X1,X2,Y1,Y2,Z1,Z2,nx,ny,nz;
  CADMBTB_getMinDistance(_pContact->_id,_pContact->_indexCAD1,_pContact->_indexCAD2,
                         X1,Y1,Z1,
                         X2,Y2,Z2,
                         nx,ny,nz,_pContact->_normalFromFace1,
                         _pContact->_dist);
  if(sPrintDist)
  {
    printf("    Minimal distance computed from CAD and n2qn1 : %lf \n",_pContact->_dist);
    printf("    Proximal point 1 computed from CAD :  X1=%lf, Y1=%lf, Z1=%lf \n",X1,Y1,Z1);
    printf("    Proximal point 2 computed from CAD :  X2=%lf, Y2=%lf, Z2=%lf \n",X2,Y2,Z2);
    if(_pContact->_normalFromFace1)
      printf("    Normal vector computed from CAD taken from  Object 1 :  nx=%lf, ny=%lf, nz=%lf \n",nx,ny,nz);
    else
      printf("    Normal vector computed from CAD taken from  Object 2 :  nx=%lf, ny=%lf, nz=%lf \n",nx,ny,nz);
  }

  // // Hack when we deal with FaceEdge and reverted=1 !!!
  // nx=-nx;
  // ny=-ny;
  // nz=-nz;



  if(_pContact->_OffsetP1)
  {
    _Pc1->setValue(0,X1+_pContact->_Offset*nx);
    _Pc1->setValue(1,Y1+_pContact->_Offset*ny);
    _Pc1->setValue(2,Z1+_pContact->_Offset*nz);
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
    _Pc2->setValue(0,X2-_pContact->_Offset*nx);
    _Pc2->setValue(1,Y2-_pContact->_Offset*ny);
    _Pc2->setValue(2,Z2-_pContact->_Offset*nz);
    if(sPrintDist)
    {
      printf("    OffSet is substracted from contact point PC2 : newPC1 =  PC2 - Offset.n ");
      printf("    OffSet %lf\n",_pContact->_Offset);
    }
  }
  /** Because in CAD model, the normal is going outside of the body .*/
  /** V.A. 11/01/2012 Which Body ?  */
  // _Nc->setValue(0,-nx);
  // _Nc->setValue(1,-ny);
  // _Nc->setValue(2,-nz);


  DEBUG_EXPR(display(););
  DEBUG_EXPR(_Nc->display(););
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););


  ACE_times[ACE_TIMER_DIST].stop();
  /** V.A. Is the following lin always true ? */

  if(_pContact->_OffsetP1)
    _pContact->_dist -= _pContact->_Offset;
  else
    _pContact->_dist += _pContact->_Offset;

  _pContact->_curTimeh=time;



  //y->setValue(0,_pContact->_dist);
  if(sPrintDist)
  {
    printf("    Gap : Signed contact distance after Offset correction : %lf \n",_pContact->_dist);
    printf("    Coordinates of contact points  after Offset corrections :\n");
    _Pc1->display();
    _Pc2->display();
    SP::SiconosVector deltaPC(new SiconosVector(3));
    *deltaPC = *_Pc1-*_Pc2;
    double realdist = (deltaPC->norm2());
    printf("    Distance between contact points =%lf \n",realdist);
    printf("    Normal vector: nx=%lf, ny=%lf, nz=%lf \n",_Nc->getValue(0),_Nc->getValue(1),_Nc->getValue(2));
  }
  //}
  y.setValue(0,_pContact->_dist);
  DEBUG_EXPR(y.display(););
  if(sPrintDist)
  {
    printf("MBTB_ContactRelation::computeh End display for contact name %s\n",_pContact->_ContactName);
  }

}





/*This function has to compute a Normal, at the extremal point. It is also the direction (P1,P2) between the two extremal points*/
// void MBTB_ContactRelation::computeJachq(double time ){

//   _jachq->setValue(0,0,-_pContact->_nX);
//   _jachq->setValue(0,1,-_pContact->_nY);
//   _jachq->setValue(0,2,-_pContact->_nZ);
//   if (_pContact->_indexBody2!=-1){
//     _jachq->setValue(0,7,_pContact->_nX);
//     _jachq->setValue(0,8,_pContact->_nY);
//     _jachq->setValue(0,9,_pContact->_nZ);
//   }
//   SP::BlockVector BlockX =boost::static_pointer_cast<BlockVector>((data[q0]));
//   for (int iDS=0;iDS<2;iDS++){
//     if (_pContact->_indexBody2==-1 && iDS ==1)
//       continue;
//     double sign = -1.0;
//     SP::SiconosVector q=(BlockX->getAllVect())[iDS];
//     //      printf("ds%d->q :",iDS);q->display();
//     ::boost::math::quaternion<float>    quatGP;
//     if (iDS==0){
//       ::boost::math::quaternion<float>    quatAux(0,_pContact->_X1-q->getValue(0),_pContact->_Y1-q->getValue(1),_pContact->_Z1-q->getValue(2));
//       quatGP=quatAux;
//     }else{
//       sign=1.0;
//       ::boost::math::quaternion<float>    quatAux(0,_pContact->_X2-q->getValue(0),_pContact->_Y2-q->getValue(1),_pContact->_Z2-q->getValue(2));
//       quatGP=quatAux;
//     }
//     //      printf("GP :%lf, %lf, %lf\n",quatGP.R_component_2(),quatGP.R_component_3(),quatGP.R_component_4());
//     ::boost::math::quaternion<float>    quatQ(q->getValue(3),q->getValue(4),q->getValue(5),q->getValue(6));
//     ::boost::math::quaternion<float>    quatcQ(q->getValue(3),-q->getValue(4),-q->getValue(5),-q->getValue(6));
//     ::boost::math::quaternion<float>    quat0(1,0,0,0);
//     ::boost::math::quaternion<float>    quatBuff;
//     quatBuff = quat0*(quatcQ*quatGP*quatQ)*quatcQ+quatQ*(quatcQ*quatGP*quatQ)*quat0;
//     _jachq->setValue(0,7*iDS+3, sign*(quatBuff.R_component_2()*_pContact->_nX + quatBuff.R_component_3()*_pContact->_nY + quatBuff.R_component_4()*_pContact->_nZ));
//     for (int i=1;i<4;i++){
//       ::boost::math::quaternion<float>    quatei(0,(i==1)?1:0,(i==2)?1:0,(i==3)?1:0);
//       quatBuff = quatei*(quatcQ*quatGP*quatQ)*quatcQ-quatQ*(quatcQ*quatGP*quatQ)*quatei;
//       //	printf("P1P1(ei)=%lf, %lf, %lf (PS=%lf)\n",quatBuff.R_component_2(),quatBuff.R_component_3(),quatBuff.R_component_4(),
//       //	       (quatBuff.R_component_2()*sOCCContacts[_indexContact]._nX + quatBuff.R_component_3()*sOCCContacts[_indexContact]._nY + quatBuff.R_component_4()*sOCCContacts[_indexContact]._nZ));
// 	_jachq->setValue(0,7*iDS+3+i, sign*(quatBuff.R_component_2()*_pContact->_nX + quatBuff.R_component_3()*_pContact->_nY + quatBuff.R_component_4()*_pContact->_nZ));
//     }
//   }
//   //    printf("computeJachq :");_jachq->display();
//   //    printf("q1dot : ");sOCCContacts[_indexContact]._DS1->dotq()->display();
//   //    printf("q2dot : ");sOCCContacts[_indexContact]._DS2->dotq()->display();
// }

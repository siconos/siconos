#include "MBTB_JointR.hpp"
#include "NewtonEulerJointR.hpp"
#include <boost/math/quaternion.hpp>
#include "op3x3.h"
//#define MBTB_JOINTR_DEBUG

MBTB_JointR::MBTB_JointR()
{
  _M.reset(new SimpleMatrix(6,6));
  _F.reset(new SiconosVector(6));
}

/*
 *
 *
 * ML_A_abs = M_FA_A + M_FB_A = F2/\BA
 * the system is :
 * FL=F1+F2
 * F1.e0 = F2.e0
 * ML_A_abs . ei = F2 /\ BA .ei for i=1,2
 *
 */

void MBTB_JointR::computeEquivalentForces()
{
  if(!_G0C1)
    return;
  double q1=_ds1->q()->getValue(3);
  double q2=_ds1->q()->getValue(4);
  double q3=_ds1->q()->getValue(5);
  double q4=_ds1->q()->getValue(6);
  ::boost::math::quaternion<double>    quattrf(q1,q2,q3,q4);
  ::boost::math::quaternion<double>    cquattrf(q1,-q2,-q3,-q4);
  ::boost::math::quaternion<double>    quatbuff;
  SP::SiconosVector Blambda=_jointR->contactForce();
  SiconosVector FL(3);
  FL.setValue(0,Blambda->getValue(0));
  FL.setValue(1,Blambda->getValue(1));
  FL.setValue(2,Blambda->getValue(2));

  SiconosVector ML_G(3);
  SiconosVector ML_G_abs(3);
  ML_G.setValue(0,Blambda->getValue(3));
  ML_G.setValue(1,Blambda->getValue(4));
  ML_G.setValue(2,Blambda->getValue(5));
  SP::SiconosVector spML_G_abs(new SiconosVector(3));
  *spML_G_abs = ML_G;
  changeFrameBodyToAbs(_ds1->q(),spML_G_abs );
  ML_G_abs = * spML_G_abs;
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces Blambda\n");
  Blambda->display();
  printf("MBTB_JointR::computeEquivalentForces ML_G_abs (in abs frame):\n");
  ML_G_abs.display();
#endif

  SiconosVector AB0(3);
  AB0=(*_G0C2)-(*_G0C1);
  ::boost::math::quaternion<double>    quatAB0(0,AB0.getValue(0),AB0.getValue(1),AB0.getValue(2));
  quatbuff = quattrf*quatAB0*cquattrf;
  SiconosVector AB(3);
  AB.setValue(0,quatbuff.R_component_2());
  AB.setValue(1,quatbuff.R_component_3());
  AB.setValue(2,quatbuff.R_component_4());
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces AB0:\n");
  AB0.display();
  printf("MBTB_JointR::computeEquivalentForces AB:\n");
  AB.display();
#endif

  ::boost::math::quaternion<double>    quatG0C1(0,_G0C1->getValue(0),_G0C1->getValue(1),_G0C1->getValue(2));
  quatbuff = quattrf*quatG0C1*cquattrf;
  SiconosVector GA(3);
  GA.setValue(0,quatbuff.R_component_2());
  GA.setValue(1,quatbuff.R_component_3());
  GA.setValue(2,quatbuff.R_component_4());
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces GA:\n");
  GA.display();
#endif
  SiconosVector ML_A_abs(3);
  SiconosVector GA_FL(3);
  cross_product(GA,FL,GA_FL);
  ML_A_abs=ML_G_abs - GA_FL;
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces GA_FL:\n");
  GA_FL.display();
  printf("MBTB_JointR::computeEquivalentForces ML_A_abs :\n");
  ML_A_abs.display();
#endif
  double normAB=1.0/AB.norm2();
  double e0_x =normAB*AB.getValue(0);
  double e0_y =normAB*AB.getValue(1);
  double e0_z =normAB*AB.getValue(2);
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces e0_x . ML_A_abs (must be zero):%e\n",
         e0_x*ML_A_abs.getValue(0)+e0_y*ML_A_abs.getValue(1)+e0_z*ML_A_abs.getValue(2));
#endif
  double  e1_x,e1_y,e1_z,e2_x,e2_y,e2_z;
  orthoBaseFromVector(&e0_x,&e0_y,&e0_z,&e1_x,&e1_y,&e1_z,&e2_x,&e2_y,&e2_z);

  /*
   * ML_A_abs = M_FA_A + M_FB_A = F2/\BA
   * the system is :
   * FL=F1+F2
   * F1.e0 = F2.e0
   * ML_A_abs . ei = F2 /\ BA .ei for i=1,2
   *
   */
  /*Fill _M*/
  _M->zero();
  //FL= F1+F2.
  _M->setValue(0,0,1);
  _M->setValue(0,3,1);
  _M->setValue(1,1,1);
  _M->setValue(1,4,1);
  _M->setValue(2,2,1);
  _M->setValue(2,5,1);
  /*F1.e0=F2.e0*/
  _M->setValue(3,0,e0_x);
  _M->setValue(3,1,e0_y);
  _M->setValue(3,2,e0_z);
  _M->setValue(3,3,-e0_x);
  _M->setValue(3,4,-e0_y);
  _M->setValue(3,5,-e0_z);

  /*
   *             F2x  BAx    F2y*BAz-F2z*BAy
   * FB/\BA . ei=F2y/\BAy.ei=F2z*BAx-F2x*BAz .ei = F2x*(BAy *ei_z-BAz*ei_y)+F2y*(BAz*ei_x-BAx*ei_z)+F2z*(BAx*ei_y-BAy*ei_x)
   *             F2z  BAz    F2x*BAy-F2y*BAx
   *
   */
  _M->setValue(4,3,-AB.getValue(1) *e1_z+AB.getValue(2)*e1_y);
  _M->setValue(4,4,-AB.getValue(2) *e1_x+AB.getValue(0)*e1_z);
  _M->setValue(4,5,-AB.getValue(0) *e1_y+AB.getValue(1)*e1_x);

  _M->setValue(5,3,-AB.getValue(1) *e2_z+AB.getValue(2)*e2_y);
  _M->setValue(5,4,-AB.getValue(2) *e2_x+AB.getValue(0)*e2_z);
  _M->setValue(5,5,-AB.getValue(0) *e2_y+AB.getValue(1)*e2_x);
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces the sytem M is:\n");
  _M->display();
#endif
  _F->setValue(0,FL.getValue(0));
  _F->setValue(1,FL.getValue(1));
  _F->setValue(2,FL.getValue(2));
  _F->setValue(3,0);
  _F->setValue(4,e1_x * ML_A_abs.getValue(0) + e1_y * ML_A_abs.getValue(1) +e1_z * ML_A_abs.getValue(2));
  _F->setValue(5,e2_x * ML_A_abs.getValue(0) + e2_y * ML_A_abs.getValue(1) +e2_z * ML_A_abs.getValue(2));
#ifdef MBTB_JOINTR_DEBUG
  printf("MBTB_JointR::computeEquivalentForces Mx=b, b:");
  _F->display();
#endif

  /*Solve the system.*/
  try
  {
    _M->PLUFactorizationInPlace();
    _M->PLUForwardBackwardInPlace(*_F);
#ifdef MBTB_JOINTR_DEBUG
    printf("MBTB_JointR::computeEquivalentForces Forces equivalent:");
    _F->display();
    printf("MBTB_JointR::computeEquivalentForces checking ML_G_abs = MF1_G+MF2_G:");
    SiconosVector Maux1(3),Maux2(3),F1(3),F2(3),GC1(3),GC2(3),dif(3);
    F1.setValue(0,_F->getValue(0));
    F1.setValue(1,_F->getValue(1));
    F1.setValue(2,_F->getValue(2));
    F2.setValue(0,_F->getValue(3));
    F2.setValue(1,_F->getValue(4));
    F2.setValue(2,_F->getValue(5));
    //FL=F1+F2
    dif = FL-F1-F2;
    printf("MBTB_JointR::computeEquivalentForces  FL-F1-F2(must be zero):\n");
    dif.display();
    //MF1_G=F1/\C1G
    cross_product(GA,F1,Maux1);
    ::boost::math::quaternion<double>    quatG0C2(0,_G0C2->getValue(0),_G0C2->getValue(1),_G0C2->getValue(2));
    quatbuff = quattrf*quatG0C2*cquattrf;
    SiconosVector GB(3);
    GB.setValue(0,quatbuff.R_component_2());
    GB.setValue(1,quatbuff.R_component_3());
    GB.setValue(2,quatbuff.R_component_4());
    cross_product(GB,F2,Maux2);
    dif=Maux1+Maux2;
    printf("MBTB_JointR::computeEquivalentForces  momentum (must be ML_G_abs):\n");
    dif.display();
    dif=Maux1+Maux2 - ML_G_abs;
    printf("MBTB_JointR::computeEquivalentForces  dif momentum(must be zero):\n");
    dif.display();
#endif
  }
  catch(const std::exception& e)
  {
    printf("MBTB_JointR: exception caught.\n");
  }


}

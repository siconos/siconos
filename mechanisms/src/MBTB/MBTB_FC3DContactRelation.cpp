/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "MBTB_FC3DContactRelation.hpp"

#include "CADMBTB_API.hpp"
#include "MBTB_Contact.hpp"
#include "MBTB_DATA.hpp"
#include "SiconosVector.hpp"
#include "ace.h"
// #define DEBUG_MESSAGES
// #define DEBUG_STDOUT
#include "siconos_debug.h"

#define OUTPUT_H_IN_FILE

siconos::mechanisms::MBTB_FC3DContactRelation::MBTB_FC3DContactRelation(
    std::shared_ptr<MBTB_Contact> pc)
    : _pContact{pc} {}

/*This function has to compute the distance between the objects*/
void siconos::mechanisms::MBTB_FC3DContactRelation::computeh(
    double time, const siconos::algebra::BlockVector& q0,
    siconos::algebra::SiconosVector& y) {
  //  DSIterator itDS=_pContact->interaction()->dynamicalSystemsBegin();
  //  auto aux = *itDS;
  // if(mbtb::data::sPrintDist)
  //  {
  //    printf("siconos::mechanisms::MBTB_FC3DContactRelation::computeh Start
  //    display for contact name %s\n",_pContact->_ContactName);
  //  }
  //  if(mbtb::data::sDS[_pContact->_indexBody1] != aux)
  //  {
  //    printf("siconos::mechanisms::MBTB_FC3DContactRelation::computeh wrong
  //    short of DS\n"); exit(1);

  //  }

  DEBUG_PRINT(
      "siconos::mechanisms::MBTB_FC3DContactRelation::computeh(double time, "
      "BlockVector& q0, SiconosVector& y )\n");
  DEBUG_EXPR(_pContact->interaction()->y(0)->display(););
  DEBUG_EXPR(y.display(););
  // auto y = _pContact->interaction()->y(0);

  // if (_pContact->_curTimeh + 1e-9 < time){
  ACE_times[ACE_TIMER_DIST].start();
  double X1, X2, Y1, Y2, Z1, Z2, n1x, n1y, n1z;
  CADMBTB_getMinDistance(
      _pContact->_id, _pContact->_indexCAD1, _pContact->_indexCAD2, X1, Y1, Z1,
      X2, Y2, Z2, n1x, n1y, n1z, _pContact->_normalFromFace1, _pContact->_dist);

  if (mbtb::data::sPrintDist) {
    printf("    Minimal distance computed from CAD and n2qn1 : %lf \n",
           _pContact->_dist);
    printf(
        "    Proximal point 1 computed from CAD :  X1=%lf, Y1=%lf, Z1=%lf \n",
        X1, Y1, Z1);
    printf(
        "    Proximal point 2 computed from CAD :  X2=%lf, Y2=%lf, Z2=%lf \n",
        X2, Y2, Z2);
    if (_pContact->_normalFromFace1)
      printf(
          "    Normal vector computed from CAD taken from  Object 1 :  nx=%lf, "
          "ny=%lf, nz=%lf \n",
          n1x, n1y, n1z);
    else
      printf(
          "    Normal vector computed from CAD taken from  Object 2 :  nx=%lf, "
          "ny=%lf, nz=%lf \n",
          n1x, n1y, n1z);
  }

  //_Pc1->setValue(0,X1); _Pc1->setValue(1,Y1); _Pc1->setValue(2,Z1);
  if (_pContact->_OffsetP1) {
    _Pc1->setValue(0, X1 + _pContact->_Offset * n1x);
    _Pc1->setValue(1, Y1 + _pContact->_Offset * n1y);
    _Pc1->setValue(2, Z1 + _pContact->_Offset * n1z);
    _Pc2->setValue(0, X2);
    _Pc2->setValue(1, Y2);
    _Pc2->setValue(2, Z2);
    if (mbtb::data::sPrintDist) {
      printf(
          "    OffSet is added from contact point PC1 : newPC1 =  PC1 + "
          "Offset.n ");
      printf("    OffSet %lf\n ", _pContact->_Offset);
    }
  } else {
    _Pc1->setValue(0, X1);
    _Pc1->setValue(1, Y1);
    _Pc1->setValue(2, Z1);
    _Pc2->setValue(0, X2 - _pContact->_Offset * n1x);
    _Pc2->setValue(1, Y2 - _pContact->_Offset * n1y);
    _Pc2->setValue(2, Z2 - _pContact->_Offset * n1z);
    if (mbtb::data::sPrintDist) {
      printf(
          "    OffSet is substracted from contact point PC2 : newPC1 =  PC2 - "
          "Offset.n ");
      printf("    point PC2 : X2=%lf,Y2=%lf,Z2=%lf",
             X2 - _pContact->_Offset * n1x, Y2 - _pContact->_Offset * n1y,
             Z2 - _pContact->_Offset * n1z);
      printf("    OffSet %lf\n", _pContact->_Offset);
    }
  }
  /*Because in CAD model, the normal is going outside of the body.*/
  _Nc->setValue(0, -n1x);
  _Nc->setValue(1, -n1y);
  _Nc->setValue(2, -n1z);

  DEBUG_EXPR(_Nc->display(););
  DEBUG_EXPR(_Pc1->display(););
  DEBUG_EXPR(_Pc2->display(););

  ACE_times[ACE_TIMER_DIST].stop();
  _pContact->_dist -= _pContact->_Offset;
  _pContact->_curTimeh = time;

  if (mbtb::data::sPrintDist) {
    std::cout << "MBTB_FC3DContactRelation compute h of "
              << _pContact->_ContactName << ": " << _pContact->_dist << "\n";

    printf(
        "MBTB_FC3DContactRelation compute h the normal rentrante : nx=%lf, "
        "ny=%lf, nz=%lf \n",
        -n1x, -n1y, -n1z);
  }

  //}
  // if (_pContact->_curTimeh + 1e-9 < time){
  //   ACE_times[ACE_TIMER_DIST].start();
  //   double X1,X2,Y1,Y2,Z1,Z2,n1x,n1y,n1z;
  //   CADMBTB_getMinDistance(_pContact->_id,_pContact->_indexCAD1,_pContact->_indexCAD2
  //   ,
  // 			   _pContact->_X1,_pContact->_Y1,_pContact->_Z1,
  // 			   _pContact->_X2,_pContact->_Y2,_pContact->_Z2,
  // 			   _pContact->_n1X,_pContact->_n1Y,_pContact->_n1Z,_pContact->_dist);
  //   ACE_times[ACE_TIMER_DIST].stop();
  //   _pContact->_dist-=_pContact->sOffset;
  //   _pContact->_curTimeh=time;
  // }
  y.setValue(0, _pContact->_dist);
  DEBUG_EXPR(y.display(););

  // auto ner
  // =(boost::static_pointer_cast<NewtonEulerR>(interaction()->relation()));
  // ner->yProj()->setValue(0,_pContact->_dist);
  //  _Pc1->setValue(0,0.5*(_pContact->_X1+_pContact->_X2));
  //  _Pc1->setValue(1,0.5*(_pContact->_Y1+_pContact->_Y2));
  //  _Pc1->setValue(2,0.5*(_pContact->_Z1+_pContact->_Z2));
  //  _Nc->setValue(0,-_pContact->_n1X);
  //  _Nc->setValue(1,-_pContact->_n1Y);
  //  _Nc->setValue(2,-_pContact->_n1Z);
  if (mbtb::data::sPrintDist) {
    std::cout << "MBTB_ContactRelation::computeh End display for contact name "
              << _pContact->_ContactName << "\n";
  }
}

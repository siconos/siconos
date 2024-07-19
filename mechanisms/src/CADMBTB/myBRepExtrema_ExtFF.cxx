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

#include "myBRepExtrema_ExtFF.hxx"

#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepExtrema_ExtCF.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <Extrema_POnSurf.hxx>
#include <Geom_Curve.hxx>
#include <Geom_Surface.hxx>
#include <Precision.hxx>
#include <Standard_Failure.hxx>
#include <StdFail_NotDone.hxx>
#include <TopAbs_ShapeEnum.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <gp_Pnt2d.hxx>

#include "CADMBTB_API.hpp"
#include "ace.h"
//=======================================================================
// function : myBRepExtrema_ExtFF
// purpose  :
//=======================================================================

myBRepExtrema_ExtFF::myBRepExtrema_ExtFF() {}

//=======================================================================
// function : myBRepExtrema_ExtFF
// purpose  :
//=======================================================================
myBRepExtrema_ExtFF::myBRepExtrema_ExtFF(const TopoDS_Face& F1, const TopoDS_Face& F2, int id1,
                                         int id2) {
  _id1 = id1;
  _id2 = id2;
  Initialize(F2);

  Perform(F1, F2);
}
myBRepExtrema_ExtFF::myBRepExtrema_ExtFF(const TopoDS_Face& F1, const TopoDS_Face& F2) {
  _id1 = -1;
  _id2 = -1;
  Initialize(F2);

  Perform(F1, F2);
}
//=======================================================================
// function : Initialize
// purpose  :
//=======================================================================

void myBRepExtrema_ExtFF::Initialize(const TopoDS_Face& F2) {
  BRepAdaptor_Surface Surf(F2);
  myHS = new BRepAdaptor_Surface(Surf);

  Standard_Real Tol = BRep_Tool::Tolerance(F2);
  Standard_Real U1, U2, V1, V2;
  BRepTools::UVBounds(F2, U1, U2, V1, V2);
  myExtrem.Initialize(myHS->Surface(), U1, U2, V1, V2, Tol);
}

//=======================================================================
// function : Perform
// purpose  :
//=======================================================================

void myBRepExtrema_ExtFF::Perform(const TopoDS_Face& F1, const TopoDS_Face& F2) {
  ACE_times[ACE_TIMER_CAD_1].start();
  ACE_times[ACE_TIMER_CAD_OK].start();
  Standard_Real U1, U2, V1, V2;
  Standard_Integer i;
  mydist.Clear();
  myPointsOnS1.Clear();
  myPointsOnS2.Clear();
  ACE_times[ACE_TIMER_CAD_OK].stop();

  ACE_times[ACE_TIMER_CAD_14].start();
  BRepAdaptor_Surface Surf1(F1); /*call UVBounds*/
  ACE_times[ACE_TIMER_CAD_14].stop();

  ACE_times[ACE_TIMER_CAD_OK].start();
  Handle(BRepAdaptor_Surface) HS1 = new BRepAdaptor_Surface(Surf1);
  // ACE_times[ACE_TIMER_CAD_OK].start();
  Standard_Real Tol1 = BRep_Tool::Tolerance(F1);
  ACE_times[ACE_TIMER_CAD_OK].stop();
  ACE_times[ACE_TIMER_CAD_13].start();
  printf("myBRepExtrema_ExtFF::Perform id1=%d, id2=%d\n", _id1, _id2);
  if (_id1 >= 0)
    CADMBTB_getUVBounds(_id1, U1, U2, V1, V2);
  else
    BRepTools::UVBounds(F1, U1, U2, V1, V2);
  ACE_times[ACE_TIMER_CAD_13].stop();
  ACE_times[ACE_TIMER_CAD_OK].start();
  // ACE_times[ACE_TIMER_CAD_2].start();
  myExtrem.Perform(HS1->Surface(), U1, U2, V1, V2, Tol1);
  // ACE_times[ACE_TIMER_CAD_2].stop();
  // ACE_times[ACE_TIMER_CAD_3].stop();
  //  exploration des points et classification:
  BRepClass_FaceClassifier classifier;
  gp_Pnt2d Puv;
  TopAbs_State state1, state2;
  Standard_Real Tol2 = BRep_Tool::Tolerance(F2);
  Extrema_POnSurf P1, P2;
  mynbext = 0;
  if (myExtrem.IsParallel()) {
    ACE_times[ACE_TIMER_CAD_OK].stop();
    mydist.Append(myExtrem.Value(1));
    mynbext = 1;
  } else {
    ACE_times[ACE_TIMER_CAD_OK].stop();
    for (i = 1; i <= myExtrem.NbExt(); i++) {
      ACE_times[ACE_TIMER_CAD_OK].start();
      myExtrem.Points(i, P1, P2);
      P1.Parameter(U1, U2);
      Puv.SetCoord(U1, U2);
      ACE_times[ACE_TIMER_CAD_OK].stop();
      ACE_times[ACE_TIMER_CAD_15].start();
      classifier.Perform(F1, Puv, Tol1); /*Call UVBounds*/
      ACE_times[ACE_TIMER_CAD_15].stop();
      ACE_times[ACE_TIMER_CAD_OK].start();
      state1 = classifier.State();
      P2.Parameter(U1, U2);
      Puv.SetCoord(U1, U2);
      classifier.Perform(F2, Puv, Tol2);
      state2 = classifier.State();
      //      if(true || (state1 == TopAbs_ON || state1 == TopAbs_IN) &&
      if (state2 == TopAbs_ON || state2 == TopAbs_IN)
      //    )
      {
        mynbext++;
        mydist.Append(myExtrem.Value(i));
        myPointsOnS1.Append(P1);
        myPointsOnS2.Append(P2);
      } else {
        cout << "myBRepExtrema_ExtFF::Perform Point out of Face\n";
      }
      ACE_times[ACE_TIMER_CAD_OK].stop();
      // ACE_times[ACE_TIMER_CAD_16].stop();
    }
  }
  // ACE_times[ACE_TIMER_CAD_3].stop();
  ACE_times[ACE_TIMER_CAD_1].stop();
}

//=======================================================================
// function : IsDone
// purpose  :
//=======================================================================

Standard_Boolean myBRepExtrema_ExtFF::IsDone() const { return myExtrem.IsDone(); }

//=======================================================================
// function : IsParallel
// purpose  :
//=======================================================================

Standard_Boolean myBRepExtrema_ExtFF::IsParallel() const { return myExtrem.IsParallel(); }

//=======================================================================
// function : NbExt
// purpose  :
//=======================================================================

Standard_Integer myBRepExtrema_ExtFF::NbExt() const {
  if (!myExtrem.IsDone()) StdFail_NotDone::Raise();
  return mynbext;
}

//=======================================================================
// function : Value
// purpose  :
//=======================================================================

Standard_Real myBRepExtrema_ExtFF::Value(const Standard_Integer N) const {
  if (!myExtrem.IsDone()) StdFail_NotDone::Raise();
  if ((N < 1) || (N > mynbext)) Standard_OutOfRange::Raise();
  return mydist.Value(N);
}

//=======================================================================
// function : ParameterOnFace1
// purpose  :
//=======================================================================

void myBRepExtrema_ExtFF::ParameterOnFace1(const Standard_Integer N, Standard_Real& U,
                                           Standard_Real& V) const {
  if (!myExtrem.IsDone()) StdFail_NotDone::Raise();
  if ((N < 1) || (N > mynbext)) Standard_OutOfRange::Raise();
  myPointsOnS1.Value(N).Parameter(U, V);
}

//=======================================================================
// function : PointOnFace1
// purpose  :
//=======================================================================

gp_Pnt myBRepExtrema_ExtFF::PointOnFace1(const Standard_Integer N) const {
  if (!myExtrem.IsDone()) StdFail_NotDone::Raise();
  if ((N < 1) || (N > mynbext)) Standard_OutOfRange::Raise();
  gp_Pnt P = myPointsOnS1.Value(N).Value();
  return P;
}

//=======================================================================
// function : ParameterOnFace2
// purpose  :
//=======================================================================

void myBRepExtrema_ExtFF::ParameterOnFace2(const Standard_Integer N, Standard_Real& U,
                                           Standard_Real& V) const {
  if (!myExtrem.IsDone()) StdFail_NotDone::Raise();
  if ((N < 1) || (N > mynbext)) Standard_OutOfRange::Raise();
  myPointsOnS2.Value(N).Parameter(U, V);
}

//=======================================================================
// function : PointOnFace1
// purpose  :
//=======================================================================

gp_Pnt myBRepExtrema_ExtFF::PointOnFace2(const Standard_Integer N) const {
  if (!myExtrem.IsDone()) StdFail_NotDone::Raise();
  if ((N < 1) || (N > mynbext)) Standard_OutOfRange::Raise();
  gp_Pnt P = myPointsOnS2.Value(N).Value();
  return P;
}

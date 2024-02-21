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
#include "cadmbtb.hpp"

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <TopoDS.hxx>
#include <gp_Dir.hxx>
#include <gp_Quaternion.hxx>
#include <gp_Vec.hxx>
#include <iostream>

#include "OccContactEdge.hpp"
#include "OccContactFace.hpp"
#include "SiconosConfig.h"
#include "SiconosException.hpp"
// #define DEBUG_MESSAGES 1
#include "SiconosFortran.h" // For n2qn1
#include "siconos_debug.h"

namespace {  // anonymous namespace for local functions

void myf_FaceFace(double* x, double* fx, double* gx, const TopoDS_Face& face1,
                  const TopoDS_Face& face2) {
  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Vec aVP2P1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  gp_Vec aV2u;
  gp_Vec aV2v;
  BRepAdaptor_Surface SF1(face1);
  BRepAdaptor_Surface SF2(face2);

  SF1.D1(x[0], x[1], aP1, aV1u, aV1v);
  SF2.D1(x[2], x[3], aP2, aV2u, aV2v);

  aVP2P1.SetX(aP1.X() - aP2.X());
  aVP2P1.SetY(aP1.Y() - aP2.Y());
  aVP2P1.SetZ(aP1.Z() - aP2.Z());
  *fx = aVP2P1.X() * aVP2P1.X() + aVP2P1.Y() * aVP2P1.Y() + aVP2P1.Z() * aVP2P1.Z();
  DEBUG_PRINTF("myf %e %e %e %e --> %e\n", x[0], x[1], x[2], x[3], *fx);
  gx[0] = 2 * aV1u.Dot(aVP2P1);
  gx[1] = 2 * aV1v.Dot(aVP2P1);
  gx[2] = -2 * aV2u.Dot(aVP2P1);
  gx[3] = -2 * aV2v.Dot(aVP2P1);
}

void myf_FaceEdge(double* x, double* fx, double* gx, const TopoDS_Face& face1,
                  const TopoDS_Edge& edge2) {
  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Vec aVP2P1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  gp_Vec aV2u;
  //  gp_Vec aV2v;/*here, zero*/
  BRepAdaptor_Surface SF1(face1);
  BRepAdaptor_Curve SC(edge2);

  SF1.D1(x[0], x[1], aP1, aV1u, aV1v);
  SC.D1(x[2], aP2, aV2u);

  aVP2P1.SetX(aP1.X() - aP2.X());
  aVP2P1.SetY(aP1.Y() - aP2.Y());
  aVP2P1.SetZ(aP1.Z() - aP2.Z());
  *fx = aVP2P1.X() * aVP2P1.X() + aVP2P1.Y() * aVP2P1.Y() + aVP2P1.Z() * aVP2P1.Z();
  DEBUG_PRINTF("myf %e %e %e %e --> %e\n", x[0], x[1], x[2], x[3], *fx);
  gx[0] = 2 * aV1u.Dot(aVP2P1);
  gx[1] = 2 * aV1v.Dot(aVP2P1);
  gx[2] = -2 * aV2u.Dot(aVP2P1);
  //  gx[3]=-2*aV2v.Dot(aVP2P1);
}
}  // namespace

gp_Pnt siconos::mechanics::occ::cadmbtb::FacePoint(const TopoDS_Face& face, Standard_Real u,
                                                   Standard_Real v) {
  // get bounds of face
  /* Standard_Real umin, umax, vmin, vmax;
   BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
   if (u < umin) u = umin;
  if (v < vmin) v = vmin;
  if (u > umax) u = umax;
  if (v > vmax) v = vmax;*/

  BRepAdaptor_Surface SF(face);
  // gp_Vec VecU,VecV;
  gp_Pnt aPaux;
  SF.D0(u, v, aPaux);  //,VecU,VecV);// compute point on surface
  return aPaux;
}

gp_Pnt siconos::mechanics::occ::cadmbtb::EdgePoint(const TopoDS_Edge& edge, Standard_Real u) {
  // get bounds of face
  // Standard_Real umin, umax, vmin, vmax;
  // BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
  // if (u < umin) u = umin;
  // if (v < vmin) v = vmin;
  // if (u > umax) u = umax;
  // if (v > vmax) v = vmax;

  BRepAdaptor_Curve SC(edge);
  //  gp_Vec VecU,VecV;
  gp_Pnt aPaux;
  SC.D0(u, aPaux);

  return aPaux;
}

gp_Dir siconos::mechanics::occ::cadmbtb::FaceNormal(const TopoDS_Face& face, Standard_Real u,
                                                    Standard_Real v) {
  // get bounds of face
  //  Standard_Real umin, umax, vmin, vmax;
  //  BRepTools::UVBounds(face, umin, umax, vmin, vmax);

  // create surface
  Handle(Geom_Surface) surf = BRep_Tool::Surface(face);

  // get surface properties
  GeomLProp_SLProps props(surf, u, v, 1, 0.01);

  // get surface normal
  gp_Dir norm = props.Normal();

  // check orientation
  //  if(face.Orientation()==TopAbs_REVERSED) norm.Reverse();
  return norm;
}

// adapted from _CADMBTB_getMinDistanceFace*_using_n2qn1  (Olivier Bonnefon)
void siconos::mechanics::occ::cadmbtb::distanceFaceFace(
    const OccContactFace& csh1, const OccContactFace& csh2, Standard_Real& X1,
    Standard_Real& Y1, Standard_Real& Z1, Standard_Real& X2, Standard_Real& Y2,
    Standard_Real& Z2, Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
    Standard_Real& MinDist) {
  // need the 2 sp pointers to keep memory

  auto pface1 = csh1.contact();
  auto pface2 = csh2.contact();

  const TopoDS_Face& face1 = *pface1;
  const TopoDS_Face& face2 = *pface2;

  double x[4];
  double f = 0;
  double g[4];
  double dxim[4];
  double df1 = 0;
  double epsabs = 0;
  double binf[4];
  double bsup[4];
  double rz[4 * (4 + 9) / 2 + 1];
  int iz[2 * 4 + 1];

  binf[0] = csh1.binf1[0];
  binf[1] = csh1.binf1[1];
  bsup[0] = csh1.bsup1[0];
  bsup[1] = csh1.bsup1[1];

  binf[2] = csh2.binf1[0];
  binf[3] = csh2.binf1[1];
  bsup[2] = csh2.bsup1[0];
  bsup[3] = csh2.bsup1[1];

  dxim[0] = 1e-6 * (bsup[0] - binf[0]);
  dxim[1] = 1e-6 * (bsup[1] - binf[1]);
  dxim[2] = 1e-6 * (bsup[2] - binf[2]);
  dxim[3] = 1e-6 * (bsup[3] - binf[3]);

  x[0] = (binf[0] + bsup[0]) * 0.5;
  x[1] = (binf[1] + bsup[1]) * 0.5;
  x[2] = (binf[2] + bsup[2]) * 0.5;
  x[3] = (binf[3] + bsup[3]) * 0.5;

  myf_FaceFace(x, &f, g, face1, face2);

  df1 = f;

  int n = 4;
  int mode = 1;

  //      DEBUG_PRINTF("call n2qn1: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,x[3]=%e,fx=%e \n
  //      g[0]=%e,g[1]=%e,g[2]=%e,g[3]=%e \n
  //      dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,dxim[3]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d
  //      \n binf[0]=%e,binf[1]=%e,binf[2]=%e,binf[3]=%e \n
  //      bsup[0]=%e,bsup[1]=%e,bsup[2]=%e,bsup[3]=%e \n
  //      sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],x[3],f,g[0],g[1],g[2],g[3],dxim[0],dxim[1],dxim[2],dxim[3],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],binf[3],bsup[0],bsup[1],bsup[2],bsup[3],sizeD,sizeI);
  siconos::fortran::optim::n2qn1(&n, x, &f, g, dxim, &df1, &epsabs, &mode, binf, bsup, iz, rz);
  while (mode > 7) {
    myf_FaceFace(x, &f, g, face1, face2);
    siconos::fortran::optim::n2qn1(&n, x, &f, g, dxim, &df1, &epsabs, &mode, binf, bsup, iz,
                                   rz);
  }

  DEBUG_PRINTF("mode=%d and min value at u=%e,v=%e f=%e\n", mode, x[0], x[1], sqrt(f));
  DEBUG_PRINTF("_CADMBTB_getMinDistanceFaceFace_using_n2qn1 dist = %e\n", sqrt(f));

  MinDist = sqrt(f);

  gp_Dir normal = FaceNormal(face2, x[2], x[3]);
  normal.Coord(nX, nY, nZ);

  /**check orientation of normal from face 2**/
  gp_Pnt aPaux1 = FacePoint(face1, x[0], x[1]);
  aPaux1.Coord(X1, Y1, Z1);
  gp_Pnt aPaux2 = FacePoint(face2, x[2], x[3]);
  aPaux2.Coord(X2, Y2, Z2);
  if (((X1 - X2) * nX + (Y1 - Y2) * nY + (Z1 - Z2) * nZ) < 0) {
    normal.Reverse();
  }
  normal.Coord(nX, nY, nZ);
}

/*idContact useful for the memory management of n2qn1.*/
void siconos::mechanics::occ::cadmbtb::distanceFaceEdge(
    const OccContactFace& csh1, const OccContactEdge& csh2, Standard_Real& X1,
    Standard_Real& Y1, Standard_Real& Z1, Standard_Real& X2, Standard_Real& Y2,
    Standard_Real& Z2, Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
    Standard_Real& MinDist) {
  // need the 2 sp pointers to keep memory
  auto pface1 = csh1.contact();
  auto pedge2 = csh2.contact();

  const TopoDS_Face& face1 = *pface1;
  const TopoDS_Edge& edge2 = *pedge2;

  int n = 3;
  double x[3];
  double f = 0;
  double g[3];
  double dxim[3];
  double df1 = 0;
  double epsabs = 0;
  double binf[3];
  double bsup[3];
  double rz[4 * (4 + 9) / 2 + 1];
  int iz[2 * 4 + 1];

  binf[0] = csh1.binf1[0];
  binf[1] = csh1.binf1[1];
  bsup[0] = csh1.bsup1[0];
  bsup[1] = csh1.bsup1[1];

  binf[2] = csh2.binf1[0];
  bsup[2] = csh2.bsup1[0];

  dxim[0] = 1e-6 * (bsup[0] - binf[0]);
  dxim[1] = 1e-6 * (bsup[1] - binf[1]);
  dxim[2] = 1e-6 * (bsup[2] - binf[2]);

  x[0] = (binf[0] + bsup[0]) * 0.5;
  x[1] = (binf[1] + bsup[1]) * 0.5;
  x[2] = (binf[2] + bsup[2]) * 0.5;
  myf_FaceEdge(x, &f, g, face1, edge2);

  df1 = f;

  int mode = 1;

  //    DEBUG_PRINTF("call n2qn1: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,fx=%e \n
  //    g[0]=%e,g[1]=%e,g[2]=%e \n
  //    dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n
  //    binf[0]=%e,binf[1]=%e,binf[2]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e \n
  //    sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],f,g[0],g[1],g[2],dxim[0],dxim[1],dxim[2],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],bsup[0],bsup[1],bsup[2],sizeD,sizeI);

  siconos::fortran::optim::n2qn1(&n, x, &f, g, dxim, &df1, &epsabs, &mode, binf, bsup, iz, rz);
  while (mode > 7) {
    myf_FaceEdge(x, &f, g, face1, edge2);

    siconos::fortran::optim::n2qn1(&n, x, &f, g, dxim, &df1, &epsabs, &mode, binf, bsup, iz,
                                   rz);
  }

  MinDist = sqrt(f);

  DEBUG_PRINTF("mode=%d and min value at u=%e,v=%e f=%e\n", mode, x[0], x[1], MinDist);
  DEBUG_PRINTF("cadmbtb_getMinDistanceFaceEdge_using_n2qn1 dist = %e\n", MinDist);

  /** V.A. Normal is always computed form the surface which is safer  */
  gp_Dir normal = FaceNormal(face1, x[0], x[1]);
  gp_Pnt aPaux = FacePoint(face1, x[0], x[1]);

  /** Coordinate of the contact point on the surface */
  aPaux.Coord(X1, Y1, Z1);
  normal.Coord(nX, nY, nZ);

  /** Coordinate of the contact point on the edge  */
  aPaux = EdgePoint(edge2, x[2]);
  aPaux.Coord(X2, Y2, Z2);

  /** check orientation of normal*/
  if (((X1 - X2) * nX + (Y1 - Y2) * nY + (Z1 - Z2) * nZ) > 0) {
    normal.Reverse();
    normal.Coord(nX, nY, nZ);
  }
}

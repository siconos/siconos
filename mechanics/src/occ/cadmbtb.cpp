#include "RuntimeException.hpp"
#include "SiconosConfig.h"
#include "cadmbtb.hpp"

#include "OccContactFace.hpp"
#include "OccContactEdge.hpp"

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <gp_Dir.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TopLoc_Location.hxx>

#include <gp_Vec.hxx>
#include <gp_Quaternion.hxx>

#include <iostream>
//#define DEBUG_MESSAGES 1
#include <debug.h>

gp_Pnt cadmbtb_FacePoint(const TopoDS_Face &face,Standard_Real u, Standard_Real v)
{
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
  SF.D0(u,v,aPaux);//,VecU,VecV);// compute point on surface
  return aPaux;
}

gp_Pnt cadmbtb_EdgePoint(const TopoDS_Edge &edge,Standard_Real u)
  {
  // get bounds of face
  //Standard_Real umin, umax, vmin, vmax;
  //BRepTools::UVBounds(face, umin, umax, vmin, vmax);          // create surface
  //if (u < umin) u = umin;
  //if (v < vmin) v = vmin;
  //if (u > umax) u = umax;
  //if (v > vmax) v = vmax;

  BRepAdaptor_Curve SC(edge);
  //  gp_Vec VecU,VecV;
  gp_Pnt aPaux;
  SC.D0(u,aPaux);

  return aPaux;
}


gp_Dir cadmbtb_FaceNormal(const TopoDS_Face &face,Standard_Real u, Standard_Real v)
{
  // get bounds of face
  //  Standard_Real umin, umax, vmin, vmax;
  //  BRepTools::UVBounds(face, umin, umax, vmin, vmax);

  // create surface
  Handle(Geom_Surface) surf=BRep_Tool::Surface(face);

  // get surface properties
  GeomLProp_SLProps props(surf, u, v, 1, 0.01);

  // get surface normal
   gp_Dir norm=props.Normal();

  // check orientation
  //  if(face.Orientation()==TopAbs_REVERSED) norm.Reverse();
   return norm;

}
#ifdef HAS_FORTRAN
extern "C"
{
  void n2qn1_(int* n, double* x, double* f, double* g, double* dxmin, double* df1, double* epsabs, int* imp,
              int *io,int* mode, int* iter, int * nsim, double* binf, double* bsup, int* iz, double* rz, int * reverse);
}
#endif
void cadmbtb_myf_FaceFace(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Face& face2);
void cadmbtb_myf_FaceFace(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Face& face2)
{
  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Vec aVP2P1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  gp_Vec aV2u;
  gp_Vec aV2v;
  BRepAdaptor_Surface SF1(face1);
  BRepAdaptor_Surface SF2(face2);

  SF1.D1(x[0],x[1],aP1, aV1u, aV1v);
  SF2.D1(x[2],x[3],aP2, aV2u, aV2v);

  aVP2P1.SetX(aP1.X()-aP2.X());
  aVP2P1.SetY(aP1.Y()-aP2.Y());
  aVP2P1.SetZ(aP1.Z()-aP2.Z());
  *fx = aVP2P1.X()*aVP2P1.X()+aVP2P1.Y()*aVP2P1.Y()+aVP2P1.Z()*aVP2P1.Z();
  DEBUG_PRINTF("myf %e %e %e %e --> %e\n",x[0],x[1],x[2],x[3],*fx);
  gx[0]=2*aV1u.Dot(aVP2P1);
  gx[1]=2*aV1v.Dot(aVP2P1);
  gx[2]=-2*aV2u.Dot(aVP2P1);
  gx[3]=-2*aV2v.Dot(aVP2P1);

}

void cadmbtb_myf_FaceEdge(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Edge& edge2);
void cadmbtb_myf_FaceEdge(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Edge& edge2)
{

  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Vec aVP2P1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  gp_Vec aV2u;
  //  gp_Vec aV2v;/*here, zero*/
  BRepAdaptor_Surface SF1(face1);
  BRepAdaptor_Curve SC(edge2);

  SF1.D1(x[0],x[1],aP1, aV1u, aV1v);
  SC.D1(x[2],aP2, aV2u);

  aVP2P1.SetX(aP1.X()-aP2.X());
  aVP2P1.SetY(aP1.Y()-aP2.Y());
  aVP2P1.SetZ(aP1.Z()-aP2.Z());
  *fx = aVP2P1.X()*aVP2P1.X()+aVP2P1.Y()*aVP2P1.Y()+aVP2P1.Z()*aVP2P1.Z();
  DEBUG_PRINTF("myf %e %e %e %e --> %e\n",x[0],x[1],x[2],x[3],*fx);
  gx[0]=2*aV1u.Dot(aVP2P1);
  gx[1]=2*aV1v.Dot(aVP2P1);
  gx[2]=-2*aV2u.Dot(aVP2P1);
  //  gx[3]=-2*aV2v.Dot(aVP2P1);

}



// adapted from _CADMBTB_getMinDistanceFace*_using_n2qn1  (Olivier Bonnefon)
void cadmbtb_distanceFaceFace(const OccContactFace& csh1,
                              const OccContactFace& csh2,
                              Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                              Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                              Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                              Standard_Real& MinDist)
{
  // need the 2 sp pointers to keep memory
  SPC::TopoDS_Face pface1 = csh1.contact();
  SPC::TopoDS_Face pface2 = csh2.contact();

  const TopoDS_Face& face1 = *pface1;
  const TopoDS_Face& face2 = *pface2;

  double x[4];
  double f = 0;
  double g[4];
  double dxim[4];
  double df1 =0;
  double epsabs=0;
  double binf[4];
  double bsup[4];
  double rz[4*(4+9)/2 + 1];
  int iz[2*4+1];

  binf[0] = csh1.binf1[0];
  binf[1] = csh1.binf1[1];
  bsup[0] = csh1.bsup1[0];
  bsup[1] = csh1.bsup1[1];

  binf[2] = csh2.binf1[0];
  binf[3] = csh2.binf1[1];
  bsup[2] = csh2.bsup1[0];
  bsup[3] = csh2.bsup1[1];

  dxim[0]=1e-6*(bsup[0]-binf[0]);
  dxim[1]=1e-6*(bsup[1]-binf[1]);
  dxim[2]=1e-6*(bsup[2]-binf[2]);
  dxim[3]=1e-6*(bsup[3]-binf[3]);

  x[0]=(binf[0]+bsup[0])*0.5;
  x[1]=(binf[1]+bsup[1])*0.5;
  x[2]=(binf[2]+bsup[2])*0.5;
  x[3]=(binf[3]+bsup[3])*0.5;

  cadmbtb_myf_FaceFace(x,&f,g,face1,face2);

  df1=f;


  int mode =1;
  int imp=0;
  int io=16;
  int iter=500;
  int nsim=3*iter;
  int reverse=1;

  int n = 4;

//      DEBUG_PRINTF("call n2qn1_: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,x[3]=%e,fx=%e \n g[0]=%e,g[1]=%e,g[2]=%e,g[3]=%e \n dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,dxim[3]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n binf[0]=%e,binf[1]=%e,binf[2]=%e,binf[3]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e,bsup[3]=%e \n sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],x[3],f,g[0],g[1],g[2],g[3],dxim[0],dxim[1],dxim[2],dxim[3],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],binf[3],bsup[0],bsup[1],bsup[2],bsup[3],sizeD,sizeI);
#ifdef HAS_FORTRAN
  n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
#else
  RuntimeException::selfThrow("_CADMBTB_getMinDistanceFaceFace_using_n2qn1, Fortran Language is not enabled in siconos mechanisms. Compile with fortran if you need n2qn1");
#endif
  while(mode > 7)
  {
    cadmbtb_myf_FaceFace(x,&f,g,face1,face2);
#ifdef HAS_FORTRAN
    n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
#else
  RuntimeException::selfThrow("_CADMBTB_getMinDistanceFaceFace_using_n2qn1, Fortran Language is not enabled in siconos mechanisms. Compile with fortran if you need n2qn1");
#endif
  }

  DEBUG_PRINTF("mode=%d and min value at u=%e,v=%e f=%e\n",mode,x[0],x[1],sqrt(f));
  DEBUG_PRINTF("_CADMBTB_getMinDistanceFaceFace_using_n2qn1 dist = %e\n",sqrt(f));

  MinDist=sqrt(f);

  gp_Dir normal = cadmbtb_FaceNormal(face2,x[2],x[3]);
  normal.Coord(nX,nY,nZ);

  /**check orientation of normal from face 2**/
  gp_Pnt aPaux1 = cadmbtb_FacePoint(face1,x[0],x[1]);
  aPaux1.Coord(X1, Y1, Z1);
  gp_Pnt aPaux2 = cadmbtb_FacePoint(face2,x[2],x[3]);
  aPaux2.Coord(X2, Y2, Z2);
  if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)<0)
  {
    normal.Reverse();
  }
  normal.Coord(nX,nY,nZ);

}



/*idContact useful for the memory management of n2qn1.*/
void cadmbtb_distanceFaceEdge(
  const OccContactFace& csh1, const OccContactEdge& csh2,
  Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
  Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
  Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
  Standard_Real& MinDist)
{

  // need the 2 sp pointers to keep memory
  SPC::TopoDS_Face pface1 = csh1.contact();
  SPC::TopoDS_Edge pedge2 = csh2.contact();

  const TopoDS_Face& face1 = *pface1;
  const TopoDS_Edge& edge2 = *pedge2;

  int n = 3;
  double x[3];
  double f = 0;
  double g[3];
  double dxim[3];
  double df1 =0;
  double epsabs=0;
  double binf[3];
  double bsup[3];
  double rz[4*(4+9)/2 + 1];
  int iz[2*4+1];

  binf[0] = csh1.binf1[0];
  binf[1] = csh1.binf1[1];
  bsup[0] = csh1.bsup1[0];
  bsup[1] = csh1.bsup1[1];

  binf[2] = csh2.binf1[0];
  bsup[2] = csh2.bsup1[0];

  dxim[0]=1e-6*(bsup[0]-binf[0]);
  dxim[1]=1e-6*(bsup[1]-binf[1]);
  dxim[2]=1e-6*(bsup[2]-binf[2]);

  x[0]=(binf[0]+bsup[0])*0.5;
  x[1]=(binf[1]+bsup[1])*0.5;
  x[2]=(binf[2]+bsup[2])*0.5;
  cadmbtb_myf_FaceEdge(x,&f,g,face1,edge2);

  df1=f;

  int mode =1;
  int imp=0;
  int io=16;
  int iter=500;
  int nsim=3*iter;
  int reverse=1;

//    DEBUG_PRINTF("call n2qn1_: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,fx=%e \n g[0]=%e,g[1]=%e,g[2]=%e \n dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n binf[0]=%e,binf[1]=%e,binf[2]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e \n sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],f,g[0],g[1],g[2],dxim[0],dxim[1],dxim[2],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],bsup[0],bsup[1],bsup[2],sizeD,sizeI);

#ifdef HAS_FORTRAN
  n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
#else
  RuntimeException::selfThrow("_CADMBTB_getMinDistanceFaceFace_using_n2qn1, Fortran Language is not enabled in siconos mechanisms. Compile with fortran if you need n2qn1");
#endif
  while(mode > 7)
  {
    cadmbtb_myf_FaceEdge(x,&f,g,face1,edge2);
#ifdef HAS_FORTRAN
    n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
#else
  RuntimeException::selfThrow("_CADMBTB_getMinDistanceFaceFace_using_n2qn1, Fortran Language is not enabled in siconos mechanisms. Compile with fortran if you need n2qn1");
#endif
  }

  MinDist=sqrt(f);

  DEBUG_PRINTF("mode=%d and min value at u=%e,v=%e f=%e\n",mode,x[0],x[1],MinDist);
  DEBUG_PRINTF("cadmbtb_getMinDistanceFaceEdge_using_n2qn1 dist = %e\n",MinDist);

  /** V.A. Normal is always computed form the surface which is safer  */
  gp_Dir normal = cadmbtb_FaceNormal(face1,x[0],x[1]);
  gp_Pnt aPaux = cadmbtb_FacePoint(face1,x[0],x[1]);

  /** Coordinate of the contact point on the surface */
  aPaux.Coord(X1, Y1, Z1);
  normal.Coord(nX,nY,nZ);

  /** Coordinate of the contact point on the edge  */
  aPaux = cadmbtb_EdgePoint(edge2,x[2]);
  aPaux.Coord(X2, Y2, Z2);

  /** check orientation of normal*/
  if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)>0)
  {
    normal.Reverse();
    normal.Coord(nX,nY,nZ);
  }

}

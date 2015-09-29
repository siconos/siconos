#include "CADMBTB_internalTools.hpp"

#include "TopExp_Explorer.hxx"
#include "TopoDS_Iterator.hxx"
#include "Geom_Surface.hxx"
#include "Geom_SphericalSurface.hxx"
#include "Geom_ToroidalSurface.hxx"
#include "Geom_RectangularTrimmedSurface.hxx"
#include "BRep_Builder.hxx"
#include "TopoDS_Face.hxx"
#include "ShapeFix_Shape.hxx"
#include "gp_Ax3.hxx"
#include "Geom_Circle.hxx"
#include "BRep_Tool.hxx"
#include "BRepExtrema_DistShapeShape.hxx"
#include "BRepAdaptor_Surface.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "GeomLProp_SLProps.hxx"
#include "TopoDS.hxx"
#include "BRepTools.hxx"
#include "TopoDS_Vertex.hxx"
#include "TopoDS_Edge.hxx"
#include "ShapeAnalysis_Surface.hxx"

#include "BRepExtrema_ExtFF.hxx"
#include "myBRepExtrema_ExtFF.hxx"
#include "ace.h"
#include "CADMBTB_DATA.hpp"
#include "gp_Pnt.hxx"
#include "gp_Vec.hxx"
#include "CADMBTB_API.hpp"
#include <assert.h>
#include "BRepClass_FaceClassifier.hxx"
#include "TopAbs_State.hxx"

#include <Standard_Version.hxx>
#if (OCC_VERSION_MAJOR >= 6 && OCC_VERSION_MINOR >= 8) || (OCC_VERSION_MAJOR == 6 && OCC_VERSION_MINOR == 7 && OCC_VERSION_MAINTENANCE >= 1)
#include <Visual3d_View.hxx>
#else
#include <Visual3d_ViewOrientation.hxx>
#endif
#include "gp_Lin.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "AIS_InteractiveContext.hxx"
#include "TopoDS_Compound.hxx"
#include "BRep_Builder.hxx"
unsigned int sCADPrintDist=0;

//#define DEBUG_USING_N2QN1
//#define CADMBTB_PRINT_DIST
gp_Pnt _CADMBTB_FacePoint(const TopoDS_Face &face,Standard_Real u, Standard_Real v)
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

gp_Pnt _CADMBTB_EdgePoint(const TopoDS_Edge &edge,Standard_Real u)
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

gp_Dir _CADMBTB_FaceNormal(const TopoDS_Face &face,Standard_Real u, Standard_Real v)
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
extern "C"
{
  void n2qn1_(int* n, double* x, double* f, double* g, double* dxmin, double* df1, double* epsabs, int* imp,
              int *io,int* mode, int* iter, int * nsim, double* binf, double* bsup, int* iz, double* rz, int * reverse);
}


void _myf_FaceFace(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Face& face2)
{
  ACE_times[ACE_TIMER_CAD_13].start();
  ACE_times[ACE_TIMER_CAD_15].start();
  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Vec aVP2P1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  gp_Vec aV2u;
  gp_Vec aV2v;
  BRepAdaptor_Surface SF1(face1);
  BRepAdaptor_Surface SF2(face2);
  ACE_times[ACE_TIMER_CAD_15].stop();
  ACE_times[ACE_TIMER_CAD_16].start();
  SF1.D1(x[0],x[1],aP1, aV1u, aV1v);
  SF2.D1(x[2],x[3],aP2, aV2u, aV2v);
  ACE_times[ACE_TIMER_CAD_16].stop();
  aVP2P1.SetX(aP1.X()-aP2.X());
  aVP2P1.SetY(aP1.Y()-aP2.Y());
  aVP2P1.SetZ(aP1.Z()-aP2.Z());
  *fx = aVP2P1.X()*aVP2P1.X()+aVP2P1.Y()*aVP2P1.Y()+aVP2P1.Z()*aVP2P1.Z();
  //printf("myf %e %e %e %e --> %e\n",x[0],x[1],x[2],x[3],*fx);
  gx[0]=2*aV1u.Dot(aVP2P1);
  gx[1]=2*aV1v.Dot(aVP2P1);
  gx[2]=-2*aV2u.Dot(aVP2P1);
  gx[3]=-2*aV2v.Dot(aVP2P1);
  ACE_times[ACE_TIMER_CAD_13].stop();
}
void _myf_FaceEdge(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Edge& edge2)
{
  ACE_times[ACE_TIMER_CAD_13].start();
  ACE_times[ACE_TIMER_CAD_15].start();
  gp_Pnt aP1;
  gp_Pnt aP2;
  gp_Vec aVP2P1;
  gp_Vec aV1u;
  gp_Vec aV1v;
  gp_Vec aV2u;
  //  gp_Vec aV2v;/*here, zero*/
  BRepAdaptor_Surface SF1(face1);
  BRepAdaptor_Curve SC(edge2);
  ACE_times[ACE_TIMER_CAD_15].stop();
  ACE_times[ACE_TIMER_CAD_16].start();
  SF1.D1(x[0],x[1],aP1, aV1u, aV1v);
  SC.D1(x[2],aP2, aV2u);
  ACE_times[ACE_TIMER_CAD_16].stop();
  aVP2P1.SetX(aP1.X()-aP2.X());
  aVP2P1.SetY(aP1.Y()-aP2.Y());
  aVP2P1.SetZ(aP1.Z()-aP2.Z());
  *fx = aVP2P1.X()*aVP2P1.X()+aVP2P1.Y()*aVP2P1.Y()+aVP2P1.Z()*aVP2P1.Z();
  //printf("myf %e %e %e %e --> %e\n",x[0],x[1],x[2],x[3],*fx);
  gx[0]=2*aV1u.Dot(aVP2P1);
  gx[1]=2*aV1v.Dot(aVP2P1);
  gx[2]=-2*aV2u.Dot(aVP2P1);
  //  gx[3]=-2*aV2v.Dot(aVP2P1);
  ACE_times[ACE_TIMER_CAD_13].stop();
}

void _CADMBTB_getMinDistanceFaceFace_using_n2qn1(unsigned int idContact, unsigned int idFace1, unsigned int idFace2,
                                                 Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                                 Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                                 Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ, 
                                                 unsigned int normalFromFace1, Standard_Real& MinDist)
{
  unsigned int idFace11=sNumberOfObj+(2*idContact-2*sNumberOfContacts);
  unsigned int idFace21=sNumberOfObj+(2*idContact+1-2*sNumberOfContacts);
  assert(idFace11==idFace1);
  assert(idFace21==idFace2);
  //AIS_Shape * sAISArtefacts[3];
  int n = 4;
  int sizeD=
    /*x*/n+
    /*dxmin*/n+
    /*df1*/n+
    /*binf*/n+
    /*bsup*/n+
    /*rz*/n*(n+9)/2 +
    /*because fortran ?*/1;

  double * sauceD = sWorkD + idContact*sizeD;
  double * pInSauceD = sauceD+1;

  int sizeI= 2*n+1 +1;
  int * sauceI = sWorkInt + idContact*sizeI;

  double * x = pInSauceD;
  pInSauceD +=n;
  double  f =0;
  double * g =pInSauceD;
  pInSauceD+=n;
  double * dxim = pInSauceD;
  pInSauceD+=n;
  double  df1 =0;
  double  epsabs=0;
  double *binf= pInSauceD;
  pInSauceD+=n;
  double *bsup= pInSauceD;
  pInSauceD+=n;
  double * rz = pInSauceD;
  int * iz =sauceI+1;

  // double start_x[4];
  // double cand_x[4];
  // /*Only the optimaling value of x has been saved. if possible, try it.*/
  // start_x[0]=x[0];
  // start_x[1]=x[1];
  // start_x[2]=x[2];
  // start_x[3]=x[3];
  // /*The current optimal will be copy in x at the end.*/
  // cand_x[0]=x[0];
  // cand_x[1]=x[1];
  // cand_x[2]=x[2];
  // cand_x[3]=x[3];

  TopoDS_Shape& sh1= sTopoDS[idFace11];
  TopoDS_Shape& sh2= sTopoDS[idFace21];
  TopExp_Explorer Ex1;
  Ex1.Init(sh1,TopAbs_FACE);
  int firstFace1=1;
  while(Ex1.More())
  {
    const TopoDS_Face& face1 = TopoDS::Face(Ex1.Current());
    TopExp_Explorer Ex2;
    Ex2.Init(sh2,TopAbs_FACE);
    int firstFace2=1;
    while(Ex2.More())
    {
      const TopoDS_Face& face2 = TopoDS::Face(Ex2.Current());
      //      ACE_times[ACE_TIMER_CAD_1].start();

      if(firstFace1)
        CADMBTB_getUVBounds(idFace1,binf[0],bsup[0],binf[1],bsup[1]);
      else
        CADMBTB_getUVBounds2(idFace1,binf[0],bsup[0],binf[1],bsup[1]);
      if(firstFace2)
        CADMBTB_getUVBounds(idFace2,binf[2],bsup[2],binf[3],bsup[3]);
      else
        CADMBTB_getUVBounds2(idFace2,binf[2],bsup[2],binf[3],bsup[3]);

      // start_x[0]=x[0];
      // start_x[1]=x[1];
      // start_x[2]=x[2];
      // start_x[3]=x[3];


      /*dxim*/

      // /*dxim*/
      dxim[0]=1e-6*(bsup[0]-binf[0]);
      dxim[1]=1e-6*(bsup[1]-binf[1]);
      dxim[2]=1e-6*(bsup[2]-binf[2]);
      dxim[3]=1e-6*(bsup[3]-binf[3]);

      // if (!(start_x[0] >= binf[0] && start_x[0] <= bsup[0] &&
      // 	  start_x[1] >= binf[1] && start_x[1] <= bsup[1] &&
      // 	  start_x[2] >= binf[2] && start_x[2] <= bsup[2] &&
      // 	    start_x[3] >= binf[3] && start_x[3] <= bsup[3] )){
      x[0]=(binf[0]+bsup[0])*0.5;
      x[1]=(binf[1]+bsup[1])*0.5;
      x[2]=(binf[2]+bsup[2])*0.5;
      x[3]=(binf[3]+bsup[3])*0.5;
      // }else{
      // 	x[0]=start_x[0];
      // 	x[1]=start_x[1];
      // 	x[2]=start_x[2];
      // 	x[3]=start_x[3];
      // }



      // epsabs=1e-4;

      // x[0]=(binf[0]+bsup[0])*0.5;
      // x[1]=(binf[1]+bsup[1])*0.5;
      // x[2]=(binf[2]+bsup[2])*0.5;
      // x[3]=(binf[3]+bsup[3])*0.5;
      _myf_FaceFace(x,&f,g,face1,face2);

      df1=f;


      int mode =1;
      int imp=0;
      int io=16;
      int iter=500;
      int nsim=3*iter;
      int reverse=1;
#ifdef DEBUG_USING_N2QN1
      printf("call n2qn1_: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,x[3]=%e,fx=%e \n g[0]=%e,g[1]=%e,g[2]=%e,g[3]=%e \n dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,dxim[3]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n binf[0]=%e,binf[1]=%e,binf[2]=%e,binf[3]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e,bsup[3]=%e \n sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],x[3],f,g[0],g[1],g[2],g[3],dxim[0],dxim[1],dxim[2],dxim[3],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],binf[3],bsup[0],bsup[1],bsup[2],bsup[3],sizeD,sizeI);
#endif
      //      ACE_times[ACE_TIMER_CAD_12].start();
      n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
      //      ACE_times[ACE_TIMER_CAD_12].stop();
      while(mode > 7)
      {
        _myf_FaceFace(x,&f,g,face1,face2);
        //	ACE_times[ACE_TIMER_CAD_12].start();
        n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
        //	ACE_times[ACE_TIMER_CAD_12].stop();
      }
      //      ACE_times[ACE_TIMER_CAD_12].stop();
      //      ACE_times[ACE_TIMER_CAD_14].start();
#ifdef DEBUG_USING_N2QN1
      printf("mode=%d and min value at u=%e,v=%e f=%e\n",mode,x[0],x[1],sqrt(f));
      printf("_CADMBTB_getMinDistanceFaceFace_using_n2qn1 dist = %e\n",sqrt(f));
#endif
      double sqrt_f=sqrt(f);
      if(MinDist>sqrt_f)
      {
        MinDist=sqrt_f;

        if(normalFromFace1)
        {
          gp_Dir normal = _CADMBTB_FaceNormal(face1,x[0],x[1]);
          normal.Coord(nX,nY,nZ);
          //check orientation of normal from face 1
	     gp_Pnt aPaux1 = _CADMBTB_FacePoint(face1,x[0],x[1]);
          aPaux1.Coord(X1, Y1, Z1);
	    gp_Pnt aPaux2 = _CADMBTB_FacePoint(face2,x[2],x[3]);
          aPaux2.Coord(X2, Y2, Z2);
	    if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)>0)
	   normal.Reverse();
	    normal.Coord(nX,nY,nZ);
	 
       	}
	else
        {
	gp_Dir normal = _CADMBTB_FaceNormal(face2,x[2],x[3]);
	normal.Coord(nX,nY,nZ); 
	/**check orientation of normal from face 2**/
	 gp_Pnt aPaux1 = _CADMBTB_FacePoint(face1,x[0],x[1]);
          aPaux1.Coord(X1, Y1, Z1);
          gp_Pnt aPaux2 = _CADMBTB_FacePoint(face2,x[2],x[3]);
          aPaux2.Coord(X2, Y2, Z2);
	   if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)<0)
	      { normal.Reverse();
	      normal.Coord(nX,nY,nZ);}
     
	    
	  /*   gp_Vec norm(normal);
           gp_Lin norm1(aPaux2,norm);
    
          BRepBuilderAPI_MakeEdge MakeEdge(norm1);
          TopoDS_Edge aVert1= BRepBuilderAPI_MakeEdge (norm1); 
	
           sAISArtefacts[idFace21] = new AIS_Shape(aVert1);
          pAIS_InteractiveContext->Display(sAISArtefacts[idFace21],false);*/ 
	  nX=-nX; 
	  nY=-nY;
	  nZ=-nZ;

         

	}
        if (sCADPrintDist)
        {
	  
          
	  /*  gp_Pnt aPaux1=_CADMBTB_FacePoint(face1,x[0],x[1]);
          aPaux1.Coord(X1, Y1, Z1);
          gp_Pnt aPaux2 = _CADMBTB_FacePoint(face2,x[2],x[3]);
          aPaux2.Coord(X2, Y2, Z2);
          BRepClass_FaceClassifier classify1(face1, aPaux1, Precision::Confusion());
          TopAbs_State st1 = classify1.State();
          BRepClass_FaceClassifier classify2(face2, aPaux2, Precision::Confusion());
          TopAbs_State st2 = classify2.State();*/
          printf("   _CADMBTB_getMinDistanceFaceFace_using_n2qn1:  Minimal distance computed from CAD and n2qn1 : %lf \n",MinDist);
          if (normalFromFace1)
	    { printf("    Normal vector computed from CAD taken from  Object 1 :  nx=%lf, ny=%lf, nz=%lf \n",nX,nY,nZ);
            printf("    First contact point taken from  Object 1 :  x1=%lf, y1=%lf, z1=%lf \n",X1,Y1,Z1);
	    printf("    Second contact point taken from  Object 2 :  x2=%lf, y2=%lf, z2=%lf \n",X2,Y2,Z2);}
         else
           { printf("    Normal vector computed from CAD taken from  Object 2 :  nx=%lf, ny=%lf, nz=%lf \n",nX,nY,nZ);
             printf("    First contact point taken from  Object 1 :  x1=%lf, y1=%lf, z1=%lf \n",X1,Y1,Z1);
	     printf("    Second contact point taken from  Object 2 :  x2=%lf, y2=%lf, z2=%lf \n",X2,Y2,Z2);}
        
	  // check status of point of face1
	  /* if (st1 == TopAbs_ON) {
	    printf("point1 sur la surface\n");
	   }else if(st1 == TopAbs_IN){
	     printf("ponit1 est interieur\n");
	   }else if(st1 == TopAbs_OUT){
	     printf("point1 est exterieur\n");
	    }else {
	     printf("unkown 4\n");
	   } 
	 //check status of point of face2
	   if (st2 == TopAbs_ON) {
	     printf("point2 sur la surface\n");
	    }else if(st2 == TopAbs_IN){
	     printf("ponit2 est interieur\n");
	   }else if(st2 == TopAbs_OUT){
	    printf("point2 est exterieur\n");
	   }else {
	    printf("unkown 4\n");
	   }*/
        }
        gp_Pnt aPaux1 = _CADMBTB_FacePoint(face1,x[0],x[1]);
        aPaux1.Coord(X1, Y1, Z1);
	gp_Pnt aPaux2 = _CADMBTB_FacePoint(face2,x[2],x[3]);
        aPaux2.Coord(X2, Y2, Z2);       
	/*	printf("    second display \n");
	printf("    First contact point taken from  Object 1 :  x1=%lf, y1=%lf, z1=%lf \n",X1,Y1,Z1);
	printf("    Second contact point taken from  Object 2 :  x2=%lf, y2=%lf, z2=%lf \n",X2,Y2,Z2);*/
        // cand_x[0]=x[0];
        // cand_x[1]=x[1];
        // cand_x[2]=x[2];
        // cand_x[3]=x[3];
      }
      //      ACE_times[ACE_TIMER_CAD_14].stop();
      //      ACE_times[ACE_TIMER_CAD_1].stop();
      firstFace2=0;
      Ex2.Next();
    }
    firstFace1=0;
    Ex1.Next();
  }
  // x[0]=cand_x[0];
  // x[1]=cand_x[1];
  // x[2]=cand_x[2];
  // x[3]=cand_x[3];

}

/*idContact useful for the memory management of n2qn1.*/
void _CADMBTB_getMinDistanceFaceEdge_using_n2qn1(
  unsigned int idContact,
  unsigned int idFace1, unsigned int idFace2,
  Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
  Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
  Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
  unsigned int normalFromFace1,
  Standard_Real& MinDist)
{
  unsigned int idFace11=sNumberOfObj+(2*idContact-2*sNumberOfContacts);
  unsigned int idFace21=sNumberOfObj+(2*idContact+1-2*sNumberOfContacts);
  unsigned int reverted=0;
  assert(idFace11==idFace1);
  assert(idFace21==idFace2);
  /** V.A.  If Face1 is the Edge, we revert Face1 and Face2.*/
  if(sTopoDSType[idFace11] == CADMBTB_TYPE_EDGE)
  {
    reverted=1;
    unsigned int aux = idFace11;
    idFace11=idFace21;
    idFace21=aux;
    if(sCADPrintDist)
    {
      printf("_CADMBTB_getMinDistanceFaceEdge_using_n2qn1 : Face1 is of type edge. Face1 <--> Face2. reverted=1 \n");
    }
  }

  TopoDS_Shape& sh1= sTopoDS[idFace11];
  TopoDS_Shape& sh2= sTopoDS[idFace21];
  TopExp_Explorer Ex1;
  Ex1.Init(sh1,TopAbs_FACE);
  TopExp_Explorer Ex2;
  Ex2.Init(sh2,TopAbs_EDGE);
  //const TopoDS_Face& face1 = TopoDS::Face(Ex1.Current());
  const TopoDS_Edge& edge2 = TopoDS::Edge(Ex2.Current());

  //    ACE_times[ACE_TIMER_CAD_1].start();
  int n = 4;
  int sizeD=
    /*x*/n+
    /*dxmin*/n+
    /*df1*/n+
    /*binf*/n+
    /*bsup*/n+
    /*rz*/n*(n+9)/2 +
    /*because fortran ?*/1;

  double * sauceD = sWorkD + idContact*sizeD;
  double * pInSauceD = sauceD+1;

  int sizeI= 2*n+1 +1;
  int * sauceI = sWorkInt + idContact*sizeI;

  double * x = pInSauceD;
  pInSauceD +=n;
  double  f =0;
  double * g =pInSauceD;
  pInSauceD+=n;
  double * dxim = pInSauceD;
  pInSauceD+=n;
  double  df1 =0;
  double  epsabs=0;
  double *binf= pInSauceD;
  pInSauceD+=n;
  double *bsup= pInSauceD;
  pInSauceD+=n;
  double * rz = pInSauceD;
  int * iz =sauceI+1;
  int firstFace1=1;

  //const TopoDS_Face& face1 = TopoDS::Face(Ex1.Current());
  while(Ex1.More())
  {
    const TopoDS_Face& face1 = TopoDS::Face(Ex1.Current());


    if(firstFace1)
      CADMBTB_getUVBounds(idFace11,binf[0],bsup[0],binf[1],bsup[1]);
    else
      CADMBTB_getUVBounds2(idFace11,binf[0],bsup[0],binf[1],bsup[1]);

    CADMBTB_getUVBounds(idFace21,binf[2],bsup[2],binf[3],bsup[3]);

    /*dxim*/

    // /*parameter domaine*/
    // CADMBTB_getUVBounds(idFace1,binf[0],bsup[0],binf[1],bsup[1]);
    // CADMBTB_getUVBounds(idFace2,binf[2],bsup[2],binf[3],bsup[3]);
    // /*dxim*/
    // dxim[0]=1e-6*(bsup[0]-binf[0]);
    // dxim[1]=1e-6*(bsup[1]-binf[1]);
    // dxim[2]=1e-6*(bsup[2]-binf[2]);
    // dxim[3]=1e-6*(bsup[3]-binf[3]);


    // epsabs=1e-4;
    /*Because of the case of cylinder, we chose to start from the midle point.*/
    x[0]=(binf[0]+bsup[0])*0.5;
    x[1]=(binf[1]+bsup[1])*0.5;
    x[2]=(binf[2]+bsup[2])*0.5;
    // x[3]=(binf[3]+bsup[3])*0.5;
    _myf_FaceEdge(x,&f,g,face1,edge2);

    df1=f;
    /*n=3 because of Face, edge.*/
    n=3;

    int mode =1;
    int imp=0;
    int io=16;
    int iter=500;
    int nsim=3*iter;
    int reverse=1;
#ifdef DEBUG_USING_N2QN1
    printf("call n2qn1_: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,fx=%e \n g[0]=%e,g[1]=%e,g[2]=%e \n dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n binf[0]=%e,binf[1]=%e,binf[2]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e \n sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],f,g[0],g[1],g[2],dxim[0],dxim[1],dxim[2],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],bsup[0],bsup[1],bsup[2],sizeD,sizeI);
#endif
    //    ACE_times[ACE_TIMER_CAD_12].start();
    n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
    //    ACE_times[ACE_TIMER_CAD_12].stop();
    while(mode > 7)
    {
      _myf_FaceEdge(x,&f,g,face1,edge2);
      //      ACE_times[ACE_TIMER_CAD_12].start();
      n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);
      //      ACE_times[ACE_TIMER_CAD_12].stop();
    }
    //    ACE_times[ACE_TIMER_CAD_12].stop();
    //    ACE_times[ACE_TIMER_CAD_14].start();
    double sqrt_f=sqrt(f);
#ifdef DEBUG_USING_N2QN1
    printf("mode=%d and min value at u=%e,v=%e f=%e\n",mode,x[0],x[1],sqrt_f);
    printf("_CADMBTB_getMinDistanceFaceEdge_using_n2qn1 dist = %e\n",sqrt_f);
#endif
    if(MinDist>sqrt_f)
    {
      /** V.A. Normal is always computed form the surface which is safer  */
      MinDist=sqrt_f;
      gp_Dir normal = _CADMBTB_FaceNormal(face1,x[0],x[1]);
      gp_Pnt aPaux = _CADMBTB_FacePoint(face1,x[0],x[1]);
      /** Coodinate of the contact point on the surface */
      aPaux.Coord(X1, Y1, Z1);
      normal.Coord(nX,nY,nZ);
      /** Coordinate of the contact point on the edge  */
      aPaux = _CADMBTB_EdgePoint(edge2,x[2]);
      aPaux.Coord(X2, Y2, Z2);
      /** check orientation of normal*/
       if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)>0) 
	{normal.Reverse();
	normal.Coord(nX,nY,nZ);}
      
    }

    Ex1.Next();
    firstFace1=0;
  }
  if(reverted)
  {    
    Standard_Real Raux;
    Raux=X1;
    X1=X2;
    X2=Raux;
    Raux=Y1;
    Y1=Y2;
    Y2=Raux;
    Raux=Z1;
    Z1=Z2;
    Z2=Raux;
    //if (normalFromFace1)
    //{
    nX=-nX;
    nY=-nY;
    nZ=-nZ;
    //}
  }
  if (sCADPrintDist)
  {
     printf(" min value at u=%e,v=%e \n",x[0],x[1]);
    printf("   _CADMBTB_getMinDistanceFaceEdge_using_n2qn1:  Minimal distance computed from CAD and n2qn1 : %lf \n",MinDist);
    if (reverted)
      printf("    Normal vector computed from CAD (reverted) :  nx=%lf, ny=%lf, nz=%lf \n",nX,nY,nZ);
    else
      printf("    Normal vector computed from CAD (not reverted) :  nx=%lf, ny=%lf, nz=%lf \n",nX,nY,nZ);
    printf("    First contact point  :  X1=%lf, Y1=%lf, Z1=%lf \n",X1,Y1,Z1);
    printf("    Second contact point  :  X2=%lf, Y2=%lf, Z2=%lf \n",X2,Y2,Z2);
    
  }
  //  ACE_times[ACE_TIMER_CAD_14].stop();
  //  ACE_times[ACE_TIMER_CAD_1].stop();
}




void _CADMBTB_getMinDistanceFaceFace(unsigned int idFace1, unsigned int idFace2,
                                     Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                     Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                     Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ, Standard_Real& MinDist)
{
  TopoDS_Shape& sh1= sTopoDS[idFace1];
  TopoDS_Shape& sh2= sTopoDS[idFace2];
  TopExp_Explorer Ex1,Ex2;
  Ex1.Init(sh1,TopAbs_FACE);
  Ex2.Init(sh2,TopAbs_FACE);
  const TopoDS_Face& aFace1 = TopoDS::Face(Ex1.Current());
  const TopoDS_Face& aFace2 = TopoDS::Face(Ex2.Current());
//Compute the parameters
  try
  {
    //ACE_times[ACE_TIMER_CAD_1].start();
    myBRepExtrema_ExtFF dst(aFace1, aFace2,idFace1, idFace2);
    printf("_CADMBTB_getMinDistanceFaceFace %d %d\n",dst._id1,dst._id2);
    printf("_CADMBTB_getMinDistanceFaceFace %d %d\n",dst._id1,dst._id2);
    //BRepExtrema_ExtFF dst (aFace1, aFace2);
    //ACE_times[ACE_TIMER_CAD_1].stop();

    Standard_Integer aN=-1;
    if(dst.IsDone())
    {

#ifdef CADMBTB_PRINT_DIST
      printf("nbSolution _CADMBTB_ getMinDistanceFaceFace %d :",dst.NbExt());
      if(dst.IsParallel())
        printf("_CADMBTB_getMinDistanceFaceFace parallel\n");
#endif
      for(int i = 1; i <= dst.NbExt(); i++)
      {



        if(MinDist > dst.Value(i))
        {
          aN = i;
          MinDist = dst.Value(i);
        }
      }
    }
    if(aN<0)
      return ;
#ifdef CADMBTB_PRINT_DIST
    printf("MinDist=%e\n",MinDist);
#endif
    Standard_Real u,v;
    dst.ParameterOnFace1(aN,u,v);
#ifdef CADMBTB_PRINT_DIST
    printf("compute Normal using _CADMBTB_FaceNormal\n");
#endif
    gp_Dir normal = _CADMBTB_FaceNormal(aFace1,u,v);
    gp_Pnt aPaux = _CADMBTB_FacePoint(aFace1,u,v);
    aPaux.Coord(X1, Y1, Z1);
    normal.Coord(nX,nY,nZ);

    dst.ParameterOnFace2(aN,u,v);
    aPaux = _CADMBTB_FacePoint(aFace2,u,v);
    aPaux.Coord(X2, Y2, Z2);
  
      
#ifdef CADMBTB_PRINT_DIST
    printf("_CADMBTB_getMinDistanceFaceFace nx %e, ny %e, nz %e .\n",nX,nY,nZ);
#endif

  }
  catch(Standard_Failure)
  {
    printf("Error in _CADMBTB_getMinDistanceFaceFace\n");
    return ;
  }



}

void _CADMBTB_getMinDistanceFaceFace(const TopoDS_Face& aFace1, const TopoDS_Face& aFace2,
                                     Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                                     Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                                     Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ, Standard_Real& MinDist)
{




  //Compute the parameters
  try
  {
    //ACE_times[ACE_TIMER_CAD_1].start();
    myBRepExtrema_ExtFF dst(aFace1, aFace2);
    //BRepExtrema_ExtFF dst (aFace1, aFace2);
    //ACE_times[ACE_TIMER_CAD_1].stop();

    Standard_Integer aN=-1;
    if(dst.IsDone())
    {

#ifdef CADMBTB_PRINT_DIST
      printf("nbSolution _CADMBTB_ getMinDistanceFaceFace %d :",dst.NbExt());
      if(dst.IsParallel())
        printf("_CADMBTB_getMinDistanceFaceFace parallel\n");
#endif
      for(int i = 1; i <= dst.NbExt(); i++)
      {



        if(MinDist > dst.Value(i))
        {
          aN = i;
          MinDist = dst.Value(i);
        }
      }
    }
    if(aN<0)
      return ;
#ifdef CADMBTB_PRINT_DIST
    printf("MinDist=%e\n",MinDist);
#endif
    Standard_Real u,v;
    dst.ParameterOnFace1(aN,u,v);
#ifdef CADMBTB_PRINT_DIST
    printf("compute Normal using _CADMBTB_FaceNormal\n");
#endif
    gp_Dir normal = _CADMBTB_FaceNormal(aFace1,u,v);
    gp_Pnt aPaux = _CADMBTB_FacePoint(aFace1,u,v);
    aPaux.Coord(X1, Y1, Z1);
    normal.Coord(nX,nY,nZ);

    dst.ParameterOnFace2(aN,u,v);
    aPaux = _CADMBTB_FacePoint(aFace2,u,v);
    aPaux.Coord(X2, Y2, Z2);
#ifdef CADMBTB_PRINT_DIST
    printf("_CADMBTB_getMinDistanceFaceFace nx %e, ny %e, nz %e .\n",nX,nY,nZ);
#endif

  }
  catch(Standard_Failure)
  {
    printf("Error in _CADMBTB_getMinDistanceFaceFace\n");
    return ;
  }

}

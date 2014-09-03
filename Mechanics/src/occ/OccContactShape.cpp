#include "OccContactShape.hpp"
#include "ContactShapeDistance.hpp"

#include <RuntimeException.hpp>

#include <TopoDS.hxx>
#include <TopExp_Explorer.hxx>
#include <BRepTools.hxx>
#include <BRep_Builder.hxx>
#include <gp_Dir.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRep_Tool.hxx>
#include <GeomLProp_SLProps.hxx>
#include <TopoDS.hxx>
#include <TopLoc_Location.hxx>

#include <gp_Vec.hxx>
#include <gp_Quaternion.hxx>

#include <limits>

OccContactShape::ContactTypeValue OccContactShape::contactType() const
{
  switch (this->_shape->ShapeType())
  {
  case TopAbs_EDGE: { return OccContactShape::Edge; }
  case TopAbs_FACE: { return OccContactShape::Face; }
  default:
    return OccContactShape::Unknown;
  };


};

void OccContactShape::computeUVBounds()
{
  TopExp_Explorer Ex1;
  TopoDS_Shape& aShape1= this->data();
  if(this->contactType() == OccContactShape::Face)
  {
    Ex1.Init(aShape1, TopAbs_FACE);
    BRepTools::UVBounds(TopoDS::Face(Ex1.Current()),
                        this->binf1[0],
                        this->bsup1[0],
                        this->binf1[1],
                        this->bsup1[1]);
    Ex1.Next();
    if(Ex1.More())
      BRepTools::UVBounds(TopoDS::Face(Ex1.Current()),
                          this->binf2[0],
                          this->bsup2[0],
                          this->binf2[1],
                          this->bsup2[1]);
  }
  else if(this->contactType() == OccContactShape::Edge)
  {
    Ex1.Init(aShape1,TopAbs_EDGE);
    const TopoDS_Edge& edge = TopoDS::Edge(Ex1.Current());
    BRepAdaptor_Curve SC(edge);
    this->binf1[0]=SC.FirstParameter();
    this->bsup1[0]=SC.LastParameter();
    this->binf1[1]=0.;
    this->bsup1[1]=0.;
  }
  else
  {
    RuntimeException::selfThrow(
      "OccContactShape::computeUVBounds() : cannot compute UV bounds for this contact shape"
      );
  }
}


//#define DEBUG_USING_N2QN1
//#define CADMBTB_PRINT_DIST
gp_Pnt cadmbtb_FacePoint(const TopoDS_Face &face,Standard_Real u, Standard_Real v);
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

gp_Pnt cadmbtb_EdgePoint(const TopoDS_Edge &edge,Standard_Real u);
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

gp_Dir cadmbtb_FaceNormal(const TopoDS_Face &face,Standard_Real u, Standard_Real v);
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
extern "C"
{
  void n2qn1_(int* n, double* x, double* f, double* g, double* dxmin, double* df1, double* epsabs, int* imp,
              int *io,int* mode, int* iter, int * nsim, double* binf, double* bsup, int* iz, double* rz, int * reverse);
}

void _myf_FaceFace(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Face& face2);
void _myf_FaceFace(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Face& face2)
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
  printf("myf %e %e %e %e --> %e\n",x[0],x[1],x[2],x[3],*fx);
  gx[0]=2*aV1u.Dot(aVP2P1);
  gx[1]=2*aV1v.Dot(aVP2P1);
  gx[2]=-2*aV2u.Dot(aVP2P1);
  gx[3]=-2*aV2v.Dot(aVP2P1);

}

void _myf_FaceEdge(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Edge& edge2);
void _myf_FaceEdge(double *x, double * fx, double * gx,const TopoDS_Face& face1,const TopoDS_Edge& edge2)
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
  //printf("myf %e %e %e %e --> %e\n",x[0],x[1],x[2],x[3],*fx);
  gx[0]=2*aV1u.Dot(aVP2P1);
  gx[1]=2*aV1v.Dot(aVP2P1);
  gx[2]=-2*aV2u.Dot(aVP2P1);
  //  gx[3]=-2*aV2v.Dot(aVP2P1);

}

void cadmbtb_distanceFaceFace(const OccContactShape& csh1,
                              const OccContactShape& csh2,
                              Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                              Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                              Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                              bool normalFromFace1,
                              Standard_Real& MinDist);

// adapted from _CADMBTB_getMinDistanceFace*_using_n2qn1  (Olivier Bonnefon)
void cadmbtb_distanceFaceFace(const OccContactShape& csh1,
                              const OccContactShape& csh2,
                              Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
                              Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
                              Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
                              bool normalFromFace1,
                              Standard_Real& MinDist)
{

  TopoDS_Shape& sh1 = csh1.data();
  TopoDS_Shape& sh2 = csh2.data();

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

  TopExp_Explorer Ex1;
  Ex1.Init(sh1, TopAbs_FACE);
  int firstFace1=1;
  while(Ex1.More())
  {
    const TopoDS_Face& face1 = TopoDS::Face(Ex1.Current());
    TopExp_Explorer Ex2;
    Ex2.Init(sh2, TopAbs_FACE);
    int firstFace2=1;
    while(Ex2.More())
    {
      const TopoDS_Face& face2 = TopoDS::Face(Ex2.Current());

      if(firstFace1)
      {
        binf[0] = csh1.binf1[0];
        binf[1] = csh1.binf1[1];
        bsup[0] = csh1.bsup1[0];
        bsup[1] = csh1.bsup1[1];
      }
      else
      {
        binf[0] = csh1.binf2[0];
        binf[1] = csh1.binf2[1];
        bsup[0] = csh1.bsup2[0];
        bsup[1] = csh1.bsup2[1];
      }
      if(firstFace2)
      {
        binf[2] = csh2.binf1[0];
        binf[3] = csh2.binf1[1];
        bsup[2] = csh2.bsup1[0];
        bsup[3] = csh2.bsup1[1];
      }
      else
      {
        binf[2] = csh2.binf2[0];
        binf[3] = csh2.binf2[1];
        bsup[2] = csh2.bsup2[0];
        bsup[3] = csh2.bsup2[1];
      }

      dxim[0]=1e-6*(bsup[0]-binf[0]);
      dxim[1]=1e-6*(bsup[1]-binf[1]);
      dxim[2]=1e-6*(bsup[2]-binf[2]);
      dxim[3]=1e-6*(bsup[3]-binf[3]);

      x[0]=(binf[0]+bsup[0])*0.5;
      x[1]=(binf[1]+bsup[1])*0.5;
      x[2]=(binf[2]+bsup[2])*0.5;
      x[3]=(binf[3]+bsup[3])*0.5;

      _myf_FaceFace(x,&f,g,face1,face2);

      df1=f;


      int mode =1;
      int imp=0;
      int io=16;
      int iter=500;
      int nsim=3*iter;
      int reverse=1;

      int n = 4;

//      DEBUG_PRINTF("call n2qn1_: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,x[3]=%e,fx=%e \n g[0]=%e,g[1]=%e,g[2]=%e,g[3]=%e \n dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,dxim[3]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n binf[0]=%e,binf[1]=%e,binf[2]=%e,binf[3]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e,bsup[3]=%e \n sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],x[3],f,g[0],g[1],g[2],g[3],dxim[0],dxim[1],dxim[2],dxim[3],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],binf[3],bsup[0],bsup[1],bsup[2],bsup[3],sizeD,sizeI);

      n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);

      while(mode > 7)
      {
        _myf_FaceFace(x,&f,g,face1,face2);

        n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);

      }

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
          gp_Dir normal = cadmbtb_FaceNormal(face1,x[0],x[1]);
          normal.Coord(nX,nY,nZ);
          //check orientation of normal from face 1
          gp_Pnt aPaux1 = cadmbtb_FacePoint(face1,x[0],x[1]);
          aPaux1.Coord(X1, Y1, Z1);
          gp_Pnt aPaux2 = cadmbtb_FacePoint(face2,x[2],x[3]);
          aPaux2.Coord(X2, Y2, Z2);
          if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)>0)
            normal.Reverse();
          normal.Coord(nX,nY,nZ);
       	}
        else
        {
          gp_Dir normal = cadmbtb_FaceNormal(face2,x[2],x[3]);
          normal.Coord(nX,nY,nZ);
          /**check orientation of normal from face 2**/
          gp_Pnt aPaux1 = cadmbtb_FacePoint(face1,x[0],x[1]);
          aPaux1.Coord(X1, Y1, Z1);
          gp_Pnt aPaux2 = cadmbtb_FacePoint(face2,x[2],x[3]);
          aPaux2.Coord(X2, Y2, Z2);
          if(((X1-X2)*nX+(Y1-Y2)*nY+(Z1-Z2)*nZ)<0)
          { normal.Reverse();
            normal.Coord(nX,nY,nZ);}


          nX=-nX;
          nY=-nY;
          nZ=-nZ;
        }

      }
      gp_Pnt aPaux1 = cadmbtb_FacePoint(face1,x[0],x[1]);
      aPaux1.Coord(X1, Y1, Z1);
      gp_Pnt aPaux2 = cadmbtb_FacePoint(face2,x[2],x[3]);
      aPaux2.Coord(X2, Y2, Z2);
      firstFace2=0;
      Ex2.Next();
    }
    firstFace1=0;
    Ex1.Next();
  }
}

void cadmbtb_distanceFaceEdge(
  const OccContactShape& sh1, const OccContactShape& sh2,
  Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
  Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
  Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
  bool normalFromFace1,
  Standard_Real& MinDist);

/*idContact useful for the memory management of n2qn1.*/
void cadmbtb_distanceFaceEdge(
  const OccContactShape& csh1, const OccContactShape& csh2,
  Standard_Real& X1, Standard_Real& Y1, Standard_Real& Z1,
  Standard_Real& X2, Standard_Real& Y2, Standard_Real& Z2,
  Standard_Real& nX, Standard_Real& nY, Standard_Real& nZ,
  bool normalFromFace1,
  Standard_Real& MinDist)
{
  TopoDS_Shape& sh1 = csh1.data();
  TopoDS_Shape& sh2 = csh2.data();

  TopExp_Explorer Ex1;
  Ex1.Init(sh1,TopAbs_FACE);
  TopExp_Explorer Ex2;
  Ex2.Init(sh2,TopAbs_EDGE);

  const TopoDS_Edge& edge2 = TopoDS::Edge(Ex2.Current());


  int n = 4;
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

  int firstFace1=1;

  while(Ex1.More())
  {
    const TopoDS_Face& face1 = TopoDS::Face(Ex1.Current());

    if(firstFace1)
    {
      binf[0] = csh1.binf1[0];
      binf[1] = csh1.binf1[1];
      bsup[0] = csh1.bsup1[0];
      bsup[1] = csh1.bsup1[1];
    }
    else
    {
      binf[0] = csh1.binf2[0];
      binf[1] = csh1.binf2[1];
      bsup[0] = csh1.bsup2[0];
      bsup[1] = csh1.bsup2[1];
    }

    binf[2] = csh2.binf1[0];
    binf[3] = csh2.binf1[1];
    bsup[2] = csh2.bsup1[0];
    bsup[3] = csh2.bsup1[1];

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

//    DEBUG_PRINTF("call n2qn1_: n=%d,x[0]=%e,x[1]=%e,x[2]=%e,fx=%e \n g[0]=%e,g[1]=%e,g[2]=%e \n dxim[0]=%e,dxim[1]=%e,dxim[2]=%e,epsabs=%e,imp=%d,io=%d,mode=%d,iter=%d,nsim=%d \n binf[0]=%e,binf[1]=%e,binf[2]=%e \n bsup[0]=%e,bsup[1]=%e,bsup[2]=%e \n sizeD=%d,sizeI=%d\n",n,x[0],x[1],x[2],f,g[0],g[1],g[2],dxim[0],dxim[1],dxim[2],epsabs,imp,io,mode,iter,nsim,binf[0],binf[1],binf[2],bsup[0],bsup[1],bsup[2],sizeD,sizeI);


    n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);

    while(mode > 7)
    {
      _myf_FaceEdge(x,&f,g,face1,edge2);

      n2qn1_(&n, x, &f, g, dxim, &df1, &epsabs, &imp, &io,&mode, &iter, &nsim, binf, bsup, iz, rz, &reverse);

    }

    double sqrt_f=sqrt(f);

#ifdef DEBUG_USING_N2QN1
    printf("mode=%d and min value at u=%e,v=%e f=%e\n",mode,x[0],x[1],sqrt_f);
    printf("cadmbtb_getMinDistanceFaceEdge_using_n2qn1 dist = %e\n",sqrt_f);
#endif

    if(MinDist>sqrt_f)
    {
      /** V.A. Normal is always computed form the surface which is safer  */
      MinDist=sqrt_f;
      gp_Dir normal = cadmbtb_FaceNormal(face1,x[0],x[1]);
      gp_Pnt aPaux = cadmbtb_FacePoint(face1,x[0],x[1]);
      /** Coodinate of the contact point on the surface */
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

    Ex1.Next();
    firstFace1=0;
  }

}



ContactShapeDistance OccContactShape::distance(
  const OccContactShape& sh2, bool normalFromFace1) const
{

  ContactShapeDistance dist;

  dist.value = std::numeric_limits<double>::infinity();

  if(this->contactType() == OccContactShape::Face &&
     sh2.contactType() == OccContactShape::Edge)
  {
    cadmbtb_distanceFaceEdge(*this, sh2,
                             dist.x1, dist.y1, dist.z1,
                             dist.x2, dist.y2, dist.z2,
                             dist.nx, dist.ny, dist.nz,
                             normalFromFace1,
                             dist.value);
  }
  else if(this->contactType() == OccContactShape::Edge &&
          sh2.contactType() == OccContactShape::Face)
  {
    cadmbtb_distanceFaceEdge(sh2, *this,
                             dist.x1, dist.y1, dist.z1,
                             dist.x2, dist.y2, dist.z2,
                             dist.nx, dist.ny, dist.nz,
                             normalFromFace1,
                             dist.value);
  }
  else if(this->contactType() == OccContactShape::Face &&
          sh2.contactType() == OccContactShape::Face)
  {
    cadmbtb_distanceFaceFace(*this, sh2,
                             dist.x1, dist.y1, dist.z1,
                             dist.x2, dist.y2, dist.z2,
                             dist.nx, dist.ny, dist.nz,
                             normalFromFace1,
                             dist.value);
  }
  else
  {
    RuntimeException::selfThrow(
      "cadmbtb_distance : cannot compute distance for this geometries");

  }

  return dist;
}

std::string OccContactShape::exportBRepAsString() const
{
  std::stringstream out;

  BRepTools::Write(this->data(), out);

  return out.str();
}

void OccContactShape::importBRepFromString(const std::string& brepstr)
{
  std::stringstream in;
  BRep_Builder brep_builder;

  in << brepstr;

  BRepTools::Read(this->data(), in, brep_builder);

}

#include <SiconosVector.hpp>

TopoDS_Shape& OccContactShape::data()
{
  if(_shape)
  {
    return *_shape;
  }
  else
  {
    _shape.reset(new TopoDS_Shape());
    return *_shape;
  }
}

void OccContactShape::move(const SiconosVector& q)
{

  const gp_Vec& translat = gp_Vec(q(0), q(1), q(2));

  const gp_Quaternion& rota = gp_Quaternion(q(4), q(5), q(6), q(3));

  gp_Trsf transfo;

  transfo.SetRotation(rota);
  transfo.SetTranslationPart(translat);

  this->data().Move(transfo);

  // cf code from Olivier
  // reset Location to avoid accumulation of TopLoc_Datum3D
  const TopLoc_Location& aLoc = this->data().Location();
  const gp_Trsf& T = aLoc.Transformation();
  TopLoc_Location aLocWithoutList(T);
  this->data().Location(aLocWithoutList);

}

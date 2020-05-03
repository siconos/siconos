#include "CADMBTB_API.hpp"
#include "CADMBTB_DATA.hpp"
#include "CADMBTB_internalTools.hpp"
#include <boost/math/quaternion.hpp>
#include "STEPControl_Reader.hxx"
#include "AIS_InteractiveContext.hxx"
#include "TopExp_Explorer.hxx"
#include "BRepTools.hxx"
#include "BRepBuilderAPI_MakeEdge.hxx"
#include "BRepBuilderAPI_MakeVertex.hxx"
#include "BRepPrim_Cylinder.hxx"
#include "BRepAdaptor_Curve.hxx"
#include "TopoDS_Compound.hxx"
#include "BRep_Builder.hxx"
#include "BRepPrim_Cone.hxx"
#include "V3d_Viewer.hxx"
#include "V3d_View.hxx"
//#include "AIS_LengthDimension.hxx"
#include "TopoDS_Edge.hxx"
//#define CADMBTB_LOAD_CONTACT
#include "Prs3d_Presentation.hxx"
#include "Prs3d_Arrow.hxx"
#include "GeomLProp_SLProps.hxx"
#include "BRep_Tool.hxx"
#include "gp_Vec.hxx"
#include "gp_Lin.hxx"
#include "GeomLProp_SLProps.hxx"
#include "TColgp_HArray1OfVec.hxx"
#include "Vrml_Normal.hxx"
unsigned int sNumberOfObj=0;
TopoDS_Shape sTopoDS[NB_OBJ];
unsigned int sTopoDSType[NB_OBJ];
AIS_Shape* spAISToposDS[NB_OBJ];
double spAISTrans[NB_OBJ];
gp_Ax3  sStartTopoDS[NB_OBJ];
gp_Trsf  sTrsfTopoDS[NB_OBJ];
Geom_Transformation  sGeomTrsf[NB_OBJ];
double sTopoDSBinf[2*NB_OBJ];
double sTopoDSBsup[2*NB_OBJ];
double sTopoDSBinf2[2*NB_OBJ];
double sTopoDSBsup2[2*NB_OBJ];
unsigned int sNumberOfArtefacts=0;
AIS_Shape * sAISArtefacts[NB_OBJ];
double sMinLineLength = 0.00001;

double * sWorkD=nullptr;
int * sWorkInt=nullptr;
int sNumberOfContacts=0;
TopoDS_Face  sFaces[2*NB_OBJ];
TopExp_Explorer Ex[2*NB_OBJ];
unsigned int sDumpGraphic=0;
// TopoDS_Shape sTopoDS1;
// gp_Ax3  sStartTopoDS1;
// gp_Trsf  sTrsfTopoDS1;
// Geom_Transformation  sGeomTrsf1;
// AIS_Shape* spAISToposDS1=0;

void CADMBTB_init(unsigned int NumberOfObj,unsigned int NumberOfContacts)
{
  assert(NumberOfObj<NB_OBJ&&"CADMBTB_init NumberOfObj to large, set NB_OBJ in CADMBTB lib");
  sNumberOfObj=NumberOfObj;
  sNumberOfContacts=NumberOfContacts;
  for(int ii=0; ii <NB_OBJ; ii++)
  {
    spAISToposDS[ii]=nullptr;
    spAISTrans[ii]=2.0;
    sTopoDSType[ii]=CADMBTB_TYPE_NONE;
  }
  if(sNumberOfObj>NB_OBJ)
    printf("***********************CADMBTB_init error too much objects\n");
  int n=4;
  int sizeD=
    /*x*/n+
    /*dxmin*/n+
    /*df1*/n+
    /*binf*/n+
    /*bsup*/n+
    /*rz*/n*(n+9)/2 +
    /*because fortran ?*/1;
  int sizeI= 2*n+1 +1;
  sWorkD=(double *) calloc(sizeD*sNumberOfContacts,sizeof(double));
  sWorkInt=(int *) malloc(sizeI*sNumberOfContacts*sizeof(int));

  return;
  // sTopoDS = (TopoDS_Shape*)malloc(sNumberOfObj*sizeof(TopoDS_Shape));
  // spAISToposDS=(AIS_Shape**)malloc(sNumberOfObj*sizeof(AIS_Shape*));
  // sStartTopoDS=(gp_Ax3 *)malloc(sNumberOfObj*sizeof(gp_Ax3));
  // sTrsfTopoDS=(gp_Trsf *)malloc(sNumberOfObj*sizeof(gp_Trsf));
  // sGeomTrsf=(Geom_Transformation *)malloc(sNumberOfObj*sizeof(Geom_Transformation));

  // /*At the begining: no orientation no translation*/
  // gp_Dir DirAuxz(0,0,1);
  // gp_Dir DirAuxx(1,0,0);
  // gp_Pnt PtAux(0,0,0);
  // gp_Ax3 Ax3Aux(PtAux,DirAuxz,DirAuxx);
  // for (int id = 0;id <NumberOfObj; id++)
  //   sStartTopoDS[id]=Ax3Aux;
  // return;

}


void CADMBTB_initContact(unsigned int idContact)
{
  assert((int)idContact<sNumberOfContacts && "CADMBTB_initContact contactId out of range");
  unsigned int idFace1=sNumberOfObj+(2*idContact-2*sNumberOfContacts);
  unsigned int idFace2=sNumberOfObj+(2*idContact+1-2*sNumberOfContacts);
  if(sTopoDSType[idFace1] == CADMBTB_TYPE_EDGE)
  {
    unsigned int aux = idFace1;
    idFace1=idFace2;
    idFace2=aux;
  }

  printf("CADMBTB_initContact id =%d,type1=%d,type2=%d.\n",idContact,sTopoDSType[idFace1],sTopoDSType[idFace2]);



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


  double * x = pInSauceD;
  pInSauceD +=n;
  //double  f =0;
  //double * g =pInSauceD;
  pInSauceD+=n;
  double * dxim = pInSauceD;
  pInSauceD+=n;
  //double  df1 =0;
  double  epsabs=0;
  double *binf= pInSauceD;
  pInSauceD+=n;
  double *bsup= pInSauceD;
  pInSauceD+=n;
  //double * rz = pInSauceD;

  /*dxim*/

  /*parameter domaine*/
  CADMBTB_getUVBounds(idFace1,binf[0],bsup[0],binf[1],bsup[1]);
  CADMBTB_getUVBounds(idFace2,binf[2],bsup[2],binf[3],bsup[3]);
  /*dxim*/
  dxim[0]=1e-6*(bsup[0]-binf[0]);
  dxim[1]=1e-6*(bsup[1]-binf[1]);
  dxim[2]=1e-6*(bsup[2]-binf[2]);
  dxim[3]=1e-6*(bsup[3]-binf[3]);


  epsabs=1e-4;

  x[0]=(binf[0]+bsup[0])*0.5;
  x[1]=(binf[1]+bsup[1])*0.5;
  x[2]=(binf[2]+bsup[2])*0.5;
  x[3]=(binf[3]+bsup[3])*0.5;



}
void CADMBTB_reset()
{
  free(sWorkD);
  free(sWorkInt);
  // free(sTopoDS);sTopoDS=0;
  // free(spAISToposDS);spAISToposDS=0;
  // free(sStartTopoDS);sStartTopoDS=0;
  // free(sTrsfTopoDS);sTrsfTopoDS=0;
  // free(sGeomTrsf);sGeomTrsf=0;
}

void CADMBTB_moveModelFromModel(unsigned int idModel1, int unsigned idModel2)
{
  assert(sNumberOfObj>idModel1&&"CADMBTB_moveModelFromModel idModel1 out of range");
  assert(sNumberOfObj>idModel2&&"CADMBTB_moveModelFromModel idModel2 out of range");
  sTopoDS[idModel1].Move(sTrsfTopoDS[idModel2]);
  /* Move is sufficient, but the following code merge
     the list contains in the TopLoc in an unique item.
     That's necessary because of performance. */
  const TopLoc_Location& aLoc = sTopoDS[idModel1].Location();
  const gp_Trsf& T = aLoc.Transformation();
  TopLoc_Location aLocWithoutList(T);
  sTopoDS[idModel1].Location(aLocWithoutList);


}
void CADMBTB_moveGraphicalModelFromModel(unsigned int idGraphicModel, unsigned int idModel)
{
  if(!pAIS_InteractiveContext || !spAISToposDS[idGraphicModel])
    return;
  assert(sNumberOfObj>idGraphicModel &&"CADMBTB_moveGraphicModelFromModel idGraphicModel out of range");
  assert(sNumberOfObj>idModel &&"CADMBTB_moveGraphicModelFromModel idModel out of range");
  // First tentative :
  // spAISToposDS[idGraphicModel]->SetTransformation(&(sGeomTrsf[idModel]),true,false);

  //Second tentative
  //Handle_AIS_Shape aAS = Handle(AIS_Shape)::DownCast(spAISToposDS[idGraphicModel]);
  //aAS->SetTransformation(&(sGeomTrsf[idModel]),true,false);

  //Third attempt

  const TopLoc_Location& aLoc = sTopoDS[idModel].Location();
  pAIS_InteractiveContext->SetLocation(spAISToposDS[idGraphicModel], aLoc);

  //spAISToposDS1->SetTransformation(&(sGeomTrsf1),false,false);
}

void CADMBTB_moveObjectFromQ(unsigned int id,double& x,double& y, double& z, double& q1,double& q2,double& q3,double& q4)
{
  assert(sNumberOfObj>id &&"CADMBTB_moveGraphicModelFromModel id out of range");
  ::boost::math::quaternion<double>    quattrf(q1,q2,q3,q4);

  //gp_Trsf aTrsf;


  ::boost::math::quaternion<double> quatZ(0,0,0,1);
  ::boost::math::quaternion<double> quatX(0,1,0,0);
  ::boost::math::quaternion<double> quatBuff(0,0,0,0);
  quatBuff=quattrf*quatZ/quattrf;
  //    std::cout<<"Z axis"<<quatBuff<<"\n";
  gp_Dir axeZ(quatBuff.R_component_2(),quatBuff.R_component_3(),quatBuff.R_component_4());

  quatBuff=quattrf*quatX/quattrf;
  //    std::cout<<"X axis"<<quatBuff<<"\n";
  gp_Dir axeX(quatBuff.R_component_2(),quatBuff.R_component_3(),quatBuff.R_component_4());
  gp_Ax3 aDestAx3(gp_Pnt(x,y,z),axeZ,axeX);
  // Set transformation
  sTrsfTopoDS[id].SetDisplacement(sStartTopoDS[id], aDestAx3);

  // Perform transformation
  sTopoDS[id].Move(sTrsfTopoDS[id]);
  /*Move is suffisient, but the following code merge the list contains in the TopLoc in an unique item.
    That's necessary because of performance.*/
  const TopLoc_Location& aLoc = sTopoDS[id].Location();
  const gp_Trsf& T = aLoc.Transformation();
  TopLoc_Location aLocWithoutList(T);
  sTopoDS[id].Location(aLocWithoutList);

  sGeomTrsf[id].SetTrsf(sTopoDS[id].Location());
  sStartTopoDS[id]= aDestAx3;

}

void CADMBTB_loadCADFile(unsigned int id, const char * fileName)
{
  assert(id < sNumberOfObj && "CADMBTB_loadCADFile id out of range");
  bool affected =false;
  printf("CADMBTB_loadCADFile id = %d using file %s.\n",id,fileName);
  STEPControl_Reader aReader;
  IFSelect_ReturnStatus status = aReader.ReadFile(fileName);
  affected=false;
  if(status == IFSelect_RetDone)
  {
    //Interface_TraceFile::SetDefault();
    bool failsonly = false;
    aReader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);
    int nbr = aReader.NbRootsForTransfer();
    aReader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);
    for(Standard_Integer n = 1; n <= nbr; n++)
    {
      bool ok;
      ok = aReader.TransferRoot(n);
      int nbs = aReader.NbShapes();
      printf("importSTEP Solid, nb shapes: %d",nbs);
      if(nbs > 0)
      {
        for(int i = 1; i <= nbs; i++)
        {
          //	TopoDS_Shape shape = aReaderManette.Shape( i );
          sTopoDS[id]=aReader.Shape(i);

          /*type Face or edge ?*/
          TopoDS_Shape& sh1= sTopoDS[id];
          TopExp_Explorer Ex1;
          int nbFaces=0;
          Ex1.Init(sh1,TopAbs_FACE);
          while(Ex1.More())
          {
            Ex1.Current();
            Ex1.Next();
            nbFaces++;
          }
          printf("CADMBTB_loadCADFile number of face=%i.\n",nbFaces);
          Ex1.Init(sh1,TopAbs_FACE);
          if(Ex1.More())
          {
            sTopoDSType[id]=CADMBTB_TYPE_FACE;
          }
          else
          {
            TopExp_Explorer Ex2;
            Ex2.Init(sh1,TopAbs_EDGE);
            if(Ex2.More())
            {
              sTopoDSType[id]=CADMBTB_TYPE_EDGE;
            }
            else
            {
              printf("CADMBTB_loadCADFile TopoDS without faces or edges.\n");
            }
          }

          //sTopoDS1=aReader.Shape( i );
          affected = true;
          //		sSManette=shape;
          printf("CADMBTB_loadCADFile %d affected, type=%d.\n",id, sTopoDSType[id]);
        }
      }
    }
  }
  if(!affected)
  {
    printf("ERRRROR CADMBTB_loadCADFile(%s)  id=%d failed.\n",fileName,id) ;
    exit(1);
  }
}
static int ca=0;
void CADMBTB_buildGraphicalModel(unsigned int id)
{
  if(!pAIS_InteractiveContext)
    return;
  assert(id < sNumberOfObj && "CADMBTB_buildGraphicModel id out of range");
  spAISToposDS[id] = new AIS_Shape(sTopoDS[id]);
  //  spAISToposDS[id]->SetColor(Quantity_NOC_PINK);
  //spAISToposDS[id]->SetColor(Quantity_NOC_BLUE1);
  //spAISToposDS[id]->SetMaterial(Graphic3d_NOM_BRONZE);
  //spAISToposDS[id]->UnsetMaterial ();
  if(ca%6==0)
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_DARKVIOLET);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  else if(ca%6==1)
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_BLUE1);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  else if(ca%6==2)
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_GREEN);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  else if(ca%6==3)
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_RED);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  else if(ca%6==4)
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_ORANGE);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  else if(ca%6==5)
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_SALMON);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  else
  {
    spAISToposDS[id]->SetColor(Quantity_NOC_YELLOW);
    spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  ca++;
  if(spAISTrans[id] > 0 && spAISTrans[id] <1)
    spAISToposDS[id]->SetTransparency(spAISTrans[id]);
  //spAISToposDS1= new AIS_Shape( sTopoDS1 );
  pAIS_InteractiveContext->Display(spAISToposDS[id], false);
}
void CADMBTB_updateGraphic()
{
  if(!pAIS_InteractiveContext)
    return;
  pAIS_InteractiveContext->UpdateCurrentViewer() ;
  if(pV3d_View && sDumpGraphic)
  {
    char file[16];
    sprintf(file,"toto%d.gif",sCmpDump);
    pV3d_View->Dump(file);
    sCmpDump++;
  }
}
void CADMBTB_setLocation(unsigned int id, double& x,double& y, double& z)
{
  assert(id < sNumberOfObj && "CADMBTB_setLocation id out of range");
  gp_Pnt PtAux2(x,y,z);
  //sStartTopoDS[id].
  sStartTopoDS[id].SetLocation(PtAux2);
}
//#define CADMBTB_PRINT_DIST
void CADMBTB_getMinDistance
(unsigned int idContact, unsigned int id1, unsigned int id2,
 double& X1, double& Y1, double& Z1,
 double& X2, double& Y2, double& Z2,
 double& nX, double& nY, double& nZ,
 unsigned int normalFromFace1,double& MinDist)
{
  assert(id1 < sNumberOfObj && "CADMBTB_getMinDistance id1 out of range");
  assert(id2 < sNumberOfObj && "CADMBTB_getMinDistance id2 out of range");
  assert((int)idContact<sNumberOfContacts && "CADMBTB_getMinDistance idContact out of range");
  MinDist = 1.e9;
  if(sTopoDSType[id2] == CADMBTB_TYPE_EDGE || sTopoDSType[id1] == CADMBTB_TYPE_EDGE)
  {
    _CADMBTB_getMinDistanceFaceEdge_using_n2qn1(idContact,id1,id2,
        X1,Y1,Z1,
        X2,Y2,Z2,
        nX,nY,nZ,
        normalFromFace1,
        MinDist);
  }
  else
  {
    _CADMBTB_getMinDistanceFaceFace_using_n2qn1(idContact,id1,id2,
        X1,Y1,Z1,
        X2,Y2,Z2,
        nX,nY,nZ,
        normalFromFace1,
        MinDist);
  }

#ifdef CADMBTB_PRINT_DIST
  printf("  CADMBTB_getMinDistance, P1(%e,%e,%e) P2(%e,%e,%e) n(%e,%e,%e): %e  \n",X1,Y1,Z1,X2,Y2,Z2,nX,nY,nZ,MinDist);
#endif
}
void CADMBTB_computeUVBounds(unsigned int id)
{
  assert(id < sNumberOfObj && "CADMBTB_computeUVBounds id out of range");
  TopExp_Explorer Ex1;
  TopoDS_Shape& aShape1=sTopoDS[id];
  if(sTopoDSType[id] == CADMBTB_TYPE_FACE)
  {
    Ex1.Init(aShape1,TopAbs_FACE);
    BRepTools::UVBounds(TopoDS::Face(Ex1.Current()), sTopoDSBinf[2*id], sTopoDSBsup[2*id], sTopoDSBinf[2*id+1], sTopoDSBsup[2*id+1]);
    Ex1.Next();
    if(Ex1.More())
      BRepTools::UVBounds(TopoDS::Face(Ex1.Current()), sTopoDSBinf2[2*id], sTopoDSBsup2[2*id], sTopoDSBinf2[2*id+1], sTopoDSBsup2[2*id+1]);
  }
  else if(sTopoDSType[id] == CADMBTB_TYPE_EDGE)
  {
    Ex1.Init(aShape1,TopAbs_EDGE);
    const TopoDS_Edge& edge = TopoDS::Edge(Ex1.Current());
    BRepAdaptor_Curve SC(edge);
    sTopoDSBinf[2*id]=SC.FirstParameter();
    sTopoDSBsup[2*id]=SC.LastParameter();
    sTopoDSBinf[2*id+1]=0;
    sTopoDSBsup[2*id+1]=0;
  }
  else
  {
    printf("CADMBTB_computeUVBounds type unknownn\n");
  }
}
/*Could be call even if the case of an edge.*/
void CADMBTB_getUVBounds(unsigned int id, double& U1, double& U2, double& V1, double& V2)
{
  assert(id < sNumberOfObj && "CADMBTB_getUVBounds id out of range");
  U1= sTopoDSBinf[2*id];
  V1= sTopoDSBinf[2*id+1];
  U2= sTopoDSBsup[2*id];
  V2= sTopoDSBsup[2*id+1];
#ifdef CADMBTB_LOAD_CONTACT
  printf("CADMBTB_getUVBounds UVBOUNDS idContact1=%d,U1=%e,U2=%e,V1=%e,V2=%e\n",id,U1,U2,V1,V2);
#endif
}
/*Could be call even if the case of an edge.*/
void CADMBTB_getUVBounds2(unsigned int id, double& U1, double& U2, double& V1, double& V2)
{
  assert(id < sNumberOfObj && "CADMBTB_getUVBounds id out of range");
  U1= sTopoDSBinf2[2*id];
  V1= sTopoDSBinf2[2*id+1];
  U2= sTopoDSBsup2[2*id];
  V2= sTopoDSBsup2[2*id+1];
#ifdef CADMBTB_LOAD_CONTACT
  printf("CADMBTB_getUVBounds2 UVBOUNDS idContact1=%d,U1=%e,U2=%e,V1=%e,V2=%e\n",id,U1,U2,V1,V2);
#endif
}

// void CADMBTB_getUBounds(unsigned int id, double& U1, double& U2){
//   assert( id < sNumberOfObj && "CADMBTB_getUBounds id out of range");
//   U1= sTopoDSBinf[2*id];
//   V1= sTopoDSBinf[2*id+1];
// #ifdef CADMBTB_LOAD_CONTACT
//   printf("CADMBTB_getUBounds UBOUNDS idContact1=%d,U1=%e,U2=%e,\n",id,U1,U2);
// #endif
// }

void CADMBTB_buildLineArtefactLine(unsigned int id,  double* X1, double* Y1, double* Z1,
                                   double* X2, double* Y2, double* Z2)
{
  if(!pAIS_InteractiveContext)
    return;
  assert(id < sNumberOfArtefacts && "CADMBTB_buildArtefactLine id out of range");
  if(sAISArtefacts[id])
  {
    pAIS_InteractiveContext->Erase(sAISArtefacts[id]);
    sAISArtefacts[id]=nullptr;
  }
  if(!X1)
    return;
  gp_Pnt P1;
  P1.SetCoord(*X1, *Y1, *Z1);
  gp_Pnt P2;
  P2.SetCoord(*X2, *Y2, *Z2);
  BRepBuilderAPI_MakeEdge MakeEdge(P1, P2);
  TopoDS_Vertex aVert1 = BRepBuilderAPI_MakeVertex(P1);
  TopoDS_Vertex aVert2 = BRepBuilderAPI_MakeVertex(P2);
  /** Color are listed in Quantity_NameOfColor.hxx header file */
  Quantity_NameOfColor         col =  Quantity_NOC_HOTPINK;
  if(MakeEdge.IsDone() && P1.Distance(P2) > sMinLineLength)
  {
    TopoDS_Compound compound;
    BRep_Builder B;
    B.MakeCompound(compound);
    B.Add(compound,MakeEdge.Edge());
    B.Add(compound,aVert1);
    B.Add(compound,aVert2);
    sAISArtefacts[id] = new AIS_Shape(compound);
    //printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e, P2=%e,%e,%e\n",X1,Y1,Z1,X2,Y2,Z2);
  }
  else
  {
    sAISArtefacts[id] = new AIS_Shape(aVert1);
    //printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e\n",X1,Y1,Z1);
  }
  sAISArtefacts[id]->SetColor(col);
  sAISArtefacts[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  pAIS_InteractiveContext->Display(sAISArtefacts[id],false);
}


void CADMBTB_buildCylinderArtefactLine(unsigned int id,  double* X1, double* Y1, double* Z1,
                                       double* X2, double* Y2, double* Z2, double *radius)
{
  if(!pAIS_InteractiveContext)
    return;
  assert(id < sNumberOfArtefacts && "CADMBTB_buildArtefactLine id out of range");
  if(sAISArtefacts[id])
  {
    pAIS_InteractiveContext->Erase(sAISArtefacts[id]);
    sAISArtefacts[id]=nullptr;
  }
  if(!X1)
    return;
  gp_Pnt P1;
  P1.SetCoord(*X1, *Y1, *Z1);
  gp_Pnt P2;
  P2.SetCoord(*X2, *Y2, *Z2);
  gp_Dir V;
  V.SetCoord(-(*X1-*X2),-(*Y1-*Y2),-(*Z1-*Z2));
  gp_Ax2 gpA2(P1,V);
  double l=sqrt((*X1-*X2)*(*X1-*X2)+(*Y1-*Y2)*(*Y1-*Y2)+(*Z1-*Z2)*(*Z1-*Z2));
  BRepPrim_Cylinder makeCyl(gpA2,*radius,l);
  gp_Ax2 gpA2b(P2,V);
  BRepPrim_Cone makeCone(gpA2b,5*(*radius),0,0.1*l);
  TopoDS_Compound compound;
  BRep_Builder B;
  B.MakeCompound(compound);
  B.Add(compound,makeCyl.Shell());
  B.Add(compound,makeCone.Shell());
  TopoDS_Vertex aVert1 = BRepBuilderAPI_MakeVertex(P1);
  if(P1.Distance(P2) > sMinLineLength)
  {
    sAISArtefacts[id] = new AIS_Shape(compound);
    sAISArtefacts[id]->SetColor(Quantity_NOC_WHITE);
    sAISArtefacts[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
    //Quantity_Color Qc(Quantity_NOC_BLUE1);
    //sAISArtefacts[id]->SetColor(Qc);
    //printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e, P2=%e,%e,%e\n",X1,Y1,Z1,X2,Y2,Z2);
  }
  else
  {
    sAISArtefacts[id] = new AIS_Shape(aVert1);
    //printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e\n",X1,Y1,Z1);
  }
  pAIS_InteractiveContext->Display(sAISArtefacts[id],false);
}

void CADMBTB_buildOrientedLineArtefactLine(unsigned int id,  double* X1, double* Y1, double* Z1,
    double* X2, double* Y2, double* Z2)
{
  if(!pAIS_InteractiveContext)
    return;
  assert(id < sNumberOfArtefacts && "CADMBTB_buildArtefactLine id out of range");
  if(sAISArtefacts[id])
  {
    pAIS_InteractiveContext->Erase(sAISArtefacts[id]);
    sAISArtefacts[id]=nullptr;
  }
  if(!X1)
    return;
  gp_Pnt P1;
  P1.SetCoord(*X1, *Y1, *Z1);
  gp_Pnt P2;
  P2.SetCoord(*X2, *Y2, *Z2);
  BRepBuilderAPI_MakeEdge MakeEdge(P1, P2);
  TopoDS_Vertex aVert1 = BRepBuilderAPI_MakeVertex(P1);
  if(MakeEdge.IsDone() && P1.Distance(P2) > sMinLineLength)
  {
    TopoDS_Compound compound;
    BRep_Builder B;
    gp_Pnt P3;
    B.MakeCompound(compound);
    B.Add(compound,MakeEdge.Edge());
    double dis = 0.1*P1.Distance(P2);
    double P1P2x= *X2-*X1;
    double P1P2y= *Y2-*Y1;
    double P1P2z= *Z2-*Z1;
    if(fabs(P1P2x)>fabs(P1P2y))
    {
      if(fabs(P1P2x)>fabs(P1P2z))
      {
        P3.SetCoord(*X2,*Y2+dis,*Z2+dis);
      }
      else
      {
        P3.SetCoord(*X2+dis,*Y2+dis,*Z2+0);
      }
    }
    else
    {
      if(fabs(P1P2y)>fabs(P1P2z))
      {
        P3.SetCoord(*X2+dis,*Y2+0,*Z2+dis);
      }
      else
      {
        P3.SetCoord(*X2+dis,*Y2+dis,*Z2+0);
      }
    }
    BRepBuilderAPI_MakeEdge MakeEdge2(P2, P3);
    B.MakeCompound(compound);
    B.Add(compound,MakeEdge.Edge());
    if(P2.Distance(P3) > sMinLineLength)
      B.Add(compound,MakeEdge2.Edge());

    sAISArtefacts[id] = new AIS_Shape(compound);
    //printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e, P2=%e,%e,%e\n",X1,Y1,Z1,X2,Y2,Z2);



    // /** V.A. Attempt yo draw arrow with OpenCascade. It Remains to do :
    //  *  + Fix the solid display of the arrow and the length
    //  *  + Manage a list a Presentation  objects, the display and the delete at each time--step.
    //  */
    //  Handle(Prs3d_Presentation) aPresentation = new Prs3d_Presentation (pAIS_InteractiveContext->CurrentViewer()->Viewer());

    // // // void Prs3d_Arrow::Draw	(	const Handle(Prs3d_Presentation)& 	aPresentation,
    // // //                           const gp_Pnt & 	aLocation,
    // // //                           const gp_Dir & 	aDirection,
    // // //                           const Quantity_PlaneAngle 	anAngle,
    // // //                           const Quantity_Length 	aLength
    // // //   )

    // gp_Vec aVec(P1,P2);
    // gp_Dir aDirection(aVec);
    // Quantity_PlaneAngle 	anAngle = 0.2;
    // Quantity_Length 	aLength	= aVec.Magnitude() ;
    // Prs3d_Arrow::Draw	(aPresentation,
    //                    P1,
    //                    aDirection,
    //                    anAngle,
    //                    aLength ) ;

    // aPresentation->Display();

  }
  else
  {
    sAISArtefacts[id] = new AIS_Shape(aVert1);
    //printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e\n",X1,Y1,Z1);
  }




  pAIS_InteractiveContext->Display(sAISArtefacts[id],true);
}

void CADMBTB_setNbOfArtefacts(unsigned int nb)
{
  for(int ii=0; ii<NB_OBJ; ii++)
    sAISArtefacts[ii]=nullptr;

  assert(nb < NB_OBJ);
  sNumberOfArtefacts = nb;
}
void CADMBTB_setContactAISdParam(unsigned int IdParam,unsigned int idContact,unsigned int idShape,double & v)
{
  assert((int)idContact<sNumberOfContacts && "CADMBTB_setContactAISdParam contactId out of range");

  unsigned int idShape1=sNumberOfObj+(2*idContact-2*sNumberOfContacts)+idShape;

  switch(IdParam)
  {
  case 0:
    spAISTrans[idShape1]=v;
    break;
  default:
    printf("Error:  CADMBTB_setAISdParam IdParam out of range.");
  }
}
void CADMBTB_loadArtefactCADFile(const char * fileName,double trans)
{
  if(!pAIS_InteractiveContext)
    return;

  printf("CADMBTB_loadArtefactCADFile using file %s.\n",fileName);
  STEPControl_Reader aReader;
  IFSelect_ReturnStatus status = aReader.ReadFile(fileName);
  if(status == IFSelect_RetDone)
  {
    //Interface_TraceFile::SetDefault();
    bool failsonly = false;
    aReader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);
    int nbr = aReader.NbRootsForTransfer();
    aReader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);
    for(Standard_Integer n = 1; n <= nbr; n++)
    {
      bool ok;
      ok = aReader.TransferRoot(n);
      int nbs = aReader.NbShapes();
      printf("importSTEP Solid, nb shapes: %d",nbs);
      if(nbs > 0)
      {
        for(int i = 1; i <= nbs; i++)
        {
          //	TopoDS_Shape shape = aReaderManette.Shape( i );
          AIS_Shape * pAis=new AIS_Shape(aReader.Shape(i));
          pAis->SetTransparency(trans);
          pAIS_InteractiveContext->Display(pAis, true);

        }
      }
    }
  }
}


TopoDS_Shape CADMBTB_TopoDS(unsigned int numDS)
{
  return sTopoDS[numDS];
}

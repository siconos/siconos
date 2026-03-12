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

#include "CADMBTB_API.hpp"

#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepPrim_Cone.hxx>
#include <BRepPrim_Cylinder.hxx>
#include <BRepTools.hxx>
#include <Geom_Transformation.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <V3d_View.hxx>
#include <boost/math/quaternion.hpp>
#include <filesystem>

#include "CADMBTB_DATA.hpp"
#include "CADMBTB_internalTools.hpp"
#include "Tools.hpp"  // For enum_to_string

unsigned int siconos::mechanisms::data::sNumberOfObj = 0;
TopoDS_Shape siconos::mechanisms::data::sTopoDS[NB_OBJ];
siconos::mechanisms::data::CADMBTB_Type siconos::mechanisms::data::sTopoDSType[NB_OBJ];
AIS_Shape* siconos::mechanisms::data::spAISToposDS[NB_OBJ];
double siconos::mechanisms::data::spAISTrans[NB_OBJ];
gp_Ax3 siconos::mechanisms::data::sStartTopoDS[NB_OBJ];
gp_Trsf siconos::mechanisms::data::sTrsfTopoDS[NB_OBJ];
Geom_Transformation siconos::mechanisms::data::sGeomTrsf[NB_OBJ];
double siconos::mechanisms::data::sTopoDSBinf[2 * NB_OBJ];
double siconos::mechanisms::data::sTopoDSBsup[2 * NB_OBJ];
double siconos::mechanisms::data::sTopoDSBinf2[2 * NB_OBJ];
double siconos::mechanisms::data::sTopoDSBsup2[2 * NB_OBJ];
unsigned int siconos::mechanisms::data::sNumberOfArtefacts = 0;
AIS_Shape* siconos::mechanisms::data::sAISArtefacts[NB_OBJ];
double siconos::mechanisms::data::sMinLineLength = 0.00001;

double* siconos::mechanisms::data::sWorkD = nullptr;
int* siconos::mechanisms::data::sWorkInt = nullptr;
int siconos::mechanisms::data::sNumberOfContacts = 0;
// TopoDS_Face siconos::mechanisms::data::sFaces[2 * NB_OBJ];
// TopExp_Explorer siconos::mechanisms::data::Ex[2 * NB_OBJ];
unsigned int siconos::mechanisms::data::sDumpGraphic = 0;
// TopoDS_Shape sTopoDS1;
// gp_Ax3  data::sStartTopoDS1;
// gp_Trsf  data::sTrsfTopoDS1;
// Geom_Transformation  data::sGeomTrsf1;
// AIS_Shape* data::spAISToposDS1=0;

void siconos::mechanisms::CADMBTB_init(unsigned int NumberOfObj,
                                       unsigned int NumberOfContacts) {
  assert(NumberOfObj < NB_OBJ &&
         "CADMBTB_init NumberOfObj to large, set NB_OBJ in CADMBTB lib");
  data::sNumberOfObj = NumberOfObj;
  data::sNumberOfContacts = NumberOfContacts;
  for (int ii = 0; ii < NB_OBJ; ii++) {
    data::spAISToposDS[ii] = nullptr;
    data::spAISTrans[ii] = 2.0;
    data::sTopoDSType[ii] = data::CADMBTB_Type::None;
  }
  if (data::sNumberOfObj > NB_OBJ)
    printf("***********************CADMBTB_init error too much objects\n");
  int n = 4;
  int sizeD =
      /*x*/ n +
      /*dxmin*/ n +
      /*df1*/ n +
      /*binf*/ n +
      /*bsup*/ n +
      /*rz*/ n * (n + 9) / 2 +
      /*because fortran ?*/ 1;
  int sizeI = 2 * n + 1 + 1;
  data::sWorkD = (double*)calloc(sizeD * data::sNumberOfContacts, sizeof(double));
  data::sWorkInt = (int*)malloc(sizeI * data::sNumberOfContacts * sizeof(int));

  std::filesystem::create_directories("graphic_dump");

  return;
  // data::sTopoDS =
  // (TopoDS_Shape*)malloc(data::sNumberOfObj*sizeof(TopoDS_Shape));
  // data::spAISToposDS=(AIS_Shape**)malloc(data::sNumberOfObj*sizeof(AIS_Shape*));
  // data::sStartTopoDS=(gp_Ax3 *)malloc(data::sNumberOfObj*sizeof(gp_Ax3));
  // data::sTrsfTopoDS=(gp_Trsf *)malloc(data::sNumberOfObj*sizeof(gp_Trsf));
  // data::sGeomTrsf=(Geom_Transformation
  // *)malloc(data::sNumberOfObj*sizeof(Geom_Transformation));

  // /*At the begining: no orientation no translation*/
  // gp_Dir DirAuxz(0,0,1);
  // gp_Dir DirAuxx(1,0,0);
  // gp_Pnt PtAux(0,0,0);
  // gp_Ax3 Ax3Aux(PtAux,DirAuxz,DirAuxx);
  // for (int id = 0;id <NumberOfObj; id++)
  //   data::sStartTopoDS[id]=Ax3Aux;
  // return;
}

void siconos::mechanisms::CADMBTB_initContact(unsigned int idContact) {
  assert((int)idContact < data::sNumberOfContacts &&
         "CADMBTB_initContact contactId out of range");
  unsigned int idFace1 = data::sNumberOfObj + (2 * idContact - 2 * data::sNumberOfContacts);
  unsigned int idFace2 =
      data::sNumberOfObj + (2 * idContact + 1 - 2 * data::sNumberOfContacts);
  if (data::sTopoDSType[idFace1] == data::CADMBTB_Type::Edge) {
    unsigned int aux = idFace1;
    idFace1 = idFace2;
    idFace2 = aux;
  }

  std::cout << "CADMBTB_initContact id=" << idContact
            << ", type1=" << siconos::tools::enum_to_string(data::sTopoDSType[idFace1])
            << ", type2=" << siconos::tools::enum_to_string(data::sTopoDSType[idFace2])
            << "\n";

  int n = 4;
  int sizeD =
      /*x*/ n +
      /*dxmin*/ n +
      /*df1*/ n +
      /*binf*/ n +
      /*bsup*/ n +
      /*rz*/ n * (n + 9) / 2 +
      /*because fortran ?*/ 1;

  double* sauceD = data::sWorkD + idContact * sizeD;
  double* pInSauceD = sauceD + 1;

  double* x = pInSauceD;
  pInSauceD += n;
  // double  f =0;
  // double * g =pInSauceD;
  pInSauceD += n;
  double* dxim = pInSauceD;
  pInSauceD += n;
  // double  df1 =0;
  double* binf = pInSauceD;
  pInSauceD += n;
  double* bsup = pInSauceD;
  pInSauceD += n;
  // double * rz = pInSauceD;

  /*dxim*/

  /*parameter domaine*/
  CADMBTB_getUVBounds(idFace1, binf[0], bsup[0], binf[1], bsup[1]);
  CADMBTB_getUVBounds(idFace2, binf[2], bsup[2], binf[3], bsup[3]);
  /*dxim*/
  dxim[0] = 1e-6 * (bsup[0] - binf[0]);
  dxim[1] = 1e-6 * (bsup[1] - binf[1]);
  dxim[2] = 1e-6 * (bsup[2] - binf[2]);
  dxim[3] = 1e-6 * (bsup[3] - binf[3]);

  x[0] = (binf[0] + bsup[0]) * 0.5;
  x[1] = (binf[1] + bsup[1]) * 0.5;
  x[2] = (binf[2] + bsup[2]) * 0.5;
  x[3] = (binf[3] + bsup[3]) * 0.5;
}
void siconos::mechanisms::CADMBTB_reset() {
  free(data::sWorkD);
  free(data::sWorkInt);
  // free(data::sTopoDS);data::sTopoDS=0;
  // free(data::spAISToposDS);data::spAISToposDS=0;
  // free(data::sStartTopoDS);data::sStartTopoDS=0;
  // free(data::sTrsfTopoDS);data::sTrsfTopoDS=0;
  // free(data::sGeomTrsf);data::sGeomTrsf=0;
}

void siconos::mechanisms::CADMBTB_moveModelFromModel(unsigned int idModel1,
                                                     int unsigned idModel2) {
  assert(data::sNumberOfObj > idModel1 && "CADMBTB_moveModelFromModel idModel1 out of range");
  assert(data::sNumberOfObj > idModel2 && "CADMBTB_moveModelFromModel idModel2 out of range");
  data::sTopoDS[idModel1].Move(data::sTrsfTopoDS[idModel2]);
  /* Move is sufficient, but the following code merge
     the list contains in the TopLoc in an unique item.
     That's necessary because of performance. */
  const TopLoc_Location& aLoc = data::sTopoDS[idModel1].Location();
  const gp_Trsf& T = aLoc.Transformation();
  TopLoc_Location aLocWithoutList(T);
  data::sTopoDS[idModel1].Location(aLocWithoutList);
}
void siconos::mechanisms::CADMBTB_moveGraphicalModelFromModel(unsigned int idGraphicModel,
                                                              unsigned int idModel) {
  if (!data::pAIS_InteractiveContext || !data::spAISToposDS[idGraphicModel]) return;
  assert(data::sNumberOfObj > idGraphicModel &&
         "CADMBTB_moveGraphicModelFromModel idGraphicModel out of range");
  assert(data::sNumberOfObj > idModel &&
         "CADMBTB_moveGraphicModelFromModel idModel out of range");
  // First tentative :
  // data::spAISToposDS[idGraphicModel]->SetTransformation(&(data::sGeomTrsf[idModel]),true,false);

  // Second tentative
  // Handle_AIS_Shape aAS =
  // Handle(AIS_Shape)::DownCast(data::spAISToposDS[idGraphicModel]);
  // aAS->SetTransformation(&(data::sGeomTrsf[idModel]),true,false);

  // Third attempt

  const TopLoc_Location& aLoc = data::sTopoDS[idModel].Location();
  data::pAIS_InteractiveContext->SetLocation(data::spAISToposDS[idGraphicModel], aLoc);

  // data::spAISToposDS1->SetTransformation(&(data::sGeomTrsf1),false,false);
}

void siconos::mechanisms::CADMBTB_moveObjectFromQ(unsigned int id, double& x, double& y,
                                                  double& z, double& q1, double& q2,
                                                  double& q3, double& q4) {
  assert(data::sNumberOfObj > id && "CADMBTB_moveGraphicModelFromModel id out of range");
  ::boost::math::quaternion<double> quattrf(q1, q2, q3, q4);

  // gp_Trsf aTrsf;

  ::boost::math::quaternion<double> quatZ(0, 0, 0, 1);
  ::boost::math::quaternion<double> quatX(0, 1, 0, 0);
  ::boost::math::quaternion<double> quatBuff(0, 0, 0, 0);
  quatBuff = quattrf * quatZ / quattrf;
  //    std::cout<<"Z axis"<<quatBuff<<"\n";
  gp_Dir axeZ(quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());

  quatBuff = quattrf * quatX / quattrf;
  //    std::cout<<"X axis"<<quatBuff<<"\n";
  gp_Dir axeX(quatBuff.R_component_2(), quatBuff.R_component_3(), quatBuff.R_component_4());
  gp_Ax3 aDestAx3(gp_Pnt(x, y, z), axeZ, axeX);
  // Set transformation
  data::sTrsfTopoDS[id].SetDisplacement(data::sStartTopoDS[id], aDestAx3);

  // Perform transformation
  data::sTopoDS[id].Move(data::sTrsfTopoDS[id]);
  /*Move is suffisient, but the following code merge the list contains in the
    TopLoc in an unique item. That's necessary because of performance.*/
  const TopLoc_Location& aLoc = data::sTopoDS[id].Location();
  const gp_Trsf& T = aLoc.Transformation();
  TopLoc_Location aLocWithoutList(T);
  data::sTopoDS[id].Location(aLocWithoutList);

  data::sGeomTrsf[id].SetTrsf(data::sTopoDS[id].Location());
  data::sStartTopoDS[id] = aDestAx3;
}

void siconos::mechanisms::CADMBTB_loadCADFile(unsigned int id, const char* fileName) {
  assert(id < data::sNumberOfObj && "CADMBTB_loadCADFile id out of range");
  bool affected = false;
  printf("CADMBTB_loadCADFile id = %d using file %s.\n", id, fileName);
  STEPControl_Reader aReader;
  IFSelect_ReturnStatus status = aReader.ReadFile(fileName);
  affected = false;
  if (status == IFSelect_RetDone) {
    // Interface_TraceFile::SetDefault();
    bool failsonly = false;
    aReader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity);
    int nbr = aReader.NbRootsForTransfer();
    aReader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity);
    for (Standard_Integer n = 1; n <= nbr; n++) {
      // bool ok = 
      aReader.TransferRoot(n);
      int nbs = aReader.NbShapes();
      printf("importSTEP Solid, nb shapes: %d", nbs);
      if (nbs > 0) {
        for (int i = 1; i <= nbs; i++) {
          //	TopoDS_Shape shape = aReaderManette.Shape( i );
          data::sTopoDS[id] = aReader.Shape(i);

          /*type Face or edge ?*/
          TopoDS_Shape& sh1 = data::sTopoDS[id];
          TopExp_Explorer Ex1;
          int nbFaces = 0;
          Ex1.Init(sh1, TopAbs_FACE);
          while (Ex1.More()) {
            Ex1.Current();
            Ex1.Next();
            nbFaces++;
          }
          printf("CADMBTB_loadCADFile number of face=%i.\n", nbFaces);
          Ex1.Init(sh1, TopAbs_FACE);
          if (Ex1.More()) {
            data::sTopoDSType[id] = data::CADMBTB_Type::Face;
          } else {
            TopExp_Explorer Ex2;
            Ex2.Init(sh1, TopAbs_EDGE);
            if (Ex2.More()) {
              data::sTopoDSType[id] = data::CADMBTB_Type::Edge;
            } else {
              printf("CADMBTB_loadCADFile TopoDS without faces or edges.\n");
            }
          }

          // data::sTopoDS1=aReader.Shape( i );
          affected = true;
          //		sSManette=shape;
          std::cout << "File " << id << " loaded, type="
                    << siconos::tools::enum_to_string(data::sTopoDSType[id]) << "\n";
        }
      }
    }
  }
  if (!affected) {
    printf("ERRRROR CADMBTB_loadCADFile(%s)  id=%d failed.\n", fileName, id);
    exit(1);
  }
}
static int ca = 0;
void siconos::mechanisms::CADMBTB_buildGraphicalModel(unsigned int id) {
  if (!data::pAIS_InteractiveContext) return;
  assert(id < data::sNumberOfObj &&
         "siconos::mechanisms::CADMBTB_buildGraphicModel id out of range");
  data::spAISToposDS[id] = new AIS_Shape(data::sTopoDS[id]);
  //  data::spAISToposDS[id]->SetColor(Quantity_NOC_PINK);
  // data::spAISToposDS[id]->SetColor(Quantity_NOC_BLUE1);
  // data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_BRONZE);
  // data::spAISToposDS[id]->UnsetMaterial ();
  if (ca % 6 == 0) {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_DARKVIOLET);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  } else if (ca % 6 == 1) {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_BLUE1);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  } else if (ca % 6 == 2) {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_GREEN);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  } else if (ca % 6 == 3) {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_RED);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  } else if (ca % 6 == 4) {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_ORANGE);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  } else if (ca % 6 == 5) {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_SALMON);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  } else {
    data::spAISToposDS[id]->SetColor(Quantity_NOC_YELLOW);
    data::spAISToposDS[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  }
  ca++;
  if (data::spAISTrans[id] > 0 && data::spAISTrans[id] < 1)
    data::spAISToposDS[id]->SetTransparency(data::spAISTrans[id]);
  // data::spAISToposDS1= new AIS_Shape( sTopoDS1 );
  data::pAIS_InteractiveContext->Display(data::spAISToposDS[id], false);
}
void siconos::mechanisms::CADMBTB_updateGraphic() {
  if (!data::pAIS_InteractiveContext) return;
  data::pAIS_InteractiveContext->UpdateCurrentViewer();

  if (data::pV3d_View && data::sDumpGraphic) {
    char file[128];
    snprintf(file, 128, "./graphic_dump/snapshot_%i.jpg", data::sCmpDump);
    data::pV3d_View->Dump(file);
    data::sCmpDump++;
  }
}
void siconos::mechanisms::CADMBTB_setLocation(unsigned int id, double& x, double& y,
                                              double& z) {
  assert(id < data::sNumberOfObj && "CADMBTB_setLocation id out of range");
  gp_Pnt PtAux2(x, y, z);
  // data::sStartTopoDS[id].
  data::sStartTopoDS[id].SetLocation(PtAux2);
}
// #define CADMBTB_PRINT_DIST
void siconos::mechanisms::CADMBTB_getMinDistance(unsigned int idContact, unsigned int id1,
                                                 unsigned int id2, double& X1, double& Y1,
                                                 double& Z1, double& X2, double& Y2,
                                                 double& Z2, double& nX, double& nY,
                                                 double& nZ, unsigned int normalFromFace1,
                                                 double& MinDist) {
  assert(id1 < data::sNumberOfObj && "CADMBTB_getMinDistance id1 out of range");
  assert(id2 < data::sNumberOfObj && "CADMBTB_getMinDistance id2 out of range");
  assert((int)idContact < data::sNumberOfContacts &&
         "CADMBTB_getMinDistance idContact out of range");
  MinDist = 1.e9;
  if (data::sTopoDSType[id2] == data::CADMBTB_Type::Edge ||
      data::sTopoDSType[id1] == data::CADMBTB_Type::Edge) {
    _CADMBTB_getMinDistanceFaceEdge_using_n2qn1(idContact, id1, id2, X1, Y1, Z1, X2, Y2, Z2,
                                                nX, nY, nZ, normalFromFace1, MinDist);
  } else {
    _CADMBTB_getMinDistanceFaceFace_using_n2qn1(idContact, id1, id2, X1, Y1, Z1, X2, Y2, Z2,
                                                nX, nY, nZ, normalFromFace1, MinDist);
  }

#ifdef CADMBTB_PRINT_DIST
  printf("  CADMBTB_getMinDistance, P1(%e,%e,%e) P2(%e,%e,%e) n(%e,%e,%e): %e  \n", X1, Y1, Z1,
         X2, Y2, Z2, nX, nY, nZ, MinDist);
#endif
}
void siconos::mechanisms::CADMBTB_computeUVBounds(unsigned int id) {
  assert(id < data::sNumberOfObj && "CADMBTB_computeUVBounds id out of range");
  TopExp_Explorer Ex1;
  TopoDS_Shape& aShape1 = data::sTopoDS[id];
  if (data::sTopoDSType[id] == data::CADMBTB_Type::Face) {
    Ex1.Init(aShape1, TopAbs_FACE);
    BRepTools::UVBounds(TopoDS::Face(Ex1.Current()), data::sTopoDSBinf[2 * id],
                        data::sTopoDSBsup[2 * id], data::sTopoDSBinf[2 * id + 1],
                        data::sTopoDSBsup[2 * id + 1]);
    Ex1.Next();
    if (Ex1.More())
      BRepTools::UVBounds(TopoDS::Face(Ex1.Current()), data::sTopoDSBinf2[2 * id],
                          data::sTopoDSBsup2[2 * id], data::sTopoDSBinf2[2 * id + 1],
                          data::sTopoDSBsup2[2 * id + 1]);
  } else if (data::sTopoDSType[id] == data::CADMBTB_Type::Edge) {
    Ex1.Init(aShape1, TopAbs_EDGE);
    const TopoDS_Edge& edge = TopoDS::Edge(Ex1.Current());
    BRepAdaptor_Curve SC(edge);
    data::sTopoDSBinf[2 * id] = SC.FirstParameter();
    data::sTopoDSBsup[2 * id] = SC.LastParameter();
    data::sTopoDSBinf[2 * id + 1] = 0;
    data::sTopoDSBsup[2 * id + 1] = 0;
  } else {
    printf("CADMBTB_computeUVBounds type unknownn\n");
  }
}
/*Could be call even if the case of an edge.*/
void siconos::mechanisms::CADMBTB_getUVBounds(unsigned int id, double& U1, double& U2,
                                              double& V1, double& V2) {
  assert(id < data::sNumberOfObj && "CADMBTB_getUVBounds id out of range");
  U1 = data::sTopoDSBinf[2 * id];
  V1 = data::sTopoDSBinf[2 * id + 1];
  U2 = data::sTopoDSBsup[2 * id];
  V2 = data::sTopoDSBsup[2 * id + 1];
#ifdef CADMBTB_LOAD_CONTACT
  printf("CADMBTB_getUVBounds UVBOUNDS idContact1=%d,U1=%e,U2=%e,V1=%e,V2=%e\n", id, U1, U2,
         V1, V2);
#endif
}
/*Could be call even if the case of an edge.*/
void siconos::mechanisms::CADMBTB_getUVBounds2(unsigned int id, double& U1, double& U2,
                                               double& V1, double& V2) {
  assert(id < data::sNumberOfObj && "CADMBTB_getUVBounds id out of range");
  U1 = data::sTopoDSBinf2[2 * id];
  V1 = data::sTopoDSBinf2[2 * id + 1];
  U2 = data::sTopoDSBsup2[2 * id];
  V2 = data::sTopoDSBsup2[2 * id + 1];
#ifdef CADMBTB_LOAD_CONTACT
  printf("CADMBTB_getUVBounds2 UVBOUNDS idContact1=%d,U1=%e,U2=%e,V1=%e,V2=%e\n", id, U1, U2,
         V1, V2);
#endif
}

// void CADMBTB_getUBounds(unsigned int id, double& U1, double& U2){
//   assert( id < data::sNumberOfObj && "CADMBTB_getUBounds id out of range");
//   U1= data::sTopoDSBinf[2*id];
//   V1= data::sTopoDSBinf[2*id+1];
// #ifdef CADMBTB_LOAD_CONTACT
//   printf("CADMBTB_getUBounds UBOUNDS idContact1=%d,U1=%e,U2=%e,\n",id,U1,U2);
// #endif
// }

void siconos::mechanisms::CADMBTB_buildLineArtefactLine(unsigned int id, double* X1,
                                                        double* Y1, double* Z1, double* X2,
                                                        double* Y2, double* Z2) {
  if (!data::pAIS_InteractiveContext) return;
  assert(id < data::sNumberOfArtefacts && "CADMBTB_buildArtefactLine id out of range");
  if (data::sAISArtefacts[id]) {
    data::pAIS_InteractiveContext->Erase(data::sAISArtefacts[id], true);
    data::sAISArtefacts[id] = nullptr;
  }
  if (!X1) return;
  gp_Pnt P1;
  P1.SetCoord(*X1, *Y1, *Z1);
  gp_Pnt P2;
  P2.SetCoord(*X2, *Y2, *Z2);
  BRepBuilderAPI_MakeEdge MakeEdge(P1, P2);
  TopoDS_Vertex aVert1 = BRepBuilderAPI_MakeVertex(P1);
  TopoDS_Vertex aVert2 = BRepBuilderAPI_MakeVertex(P2);
  /** Color are listed in Quantity_NameOfColor.hxx header file */
  Quantity_NameOfColor col = Quantity_NOC_HOTPINK;
  if (MakeEdge.IsDone() && P1.Distance(P2) > data::sMinLineLength) {
    TopoDS_Compound compound;
    BRep_Builder B;
    B.MakeCompound(compound);
    B.Add(compound, MakeEdge.Edge());
    B.Add(compound, aVert1);
    B.Add(compound, aVert2);
    data::sAISArtefacts[id] = new AIS_Shape(compound);
    // printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e,
    // P2=%e,%e,%e\n",X1,Y1,Z1,X2,Y2,Z2);
  } else {
    data::sAISArtefacts[id] = new AIS_Shape(aVert1);
    // printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e\n",X1,Y1,Z1);
  }
  data::sAISArtefacts[id]->SetColor(col);
  data::sAISArtefacts[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
  data::pAIS_InteractiveContext->Display(data::sAISArtefacts[id], false);
}

void siconos::mechanisms::CADMBTB_buildCylinderArtefactLine(unsigned int id, double* X1,
                                                            double* Y1, double* Z1, double* X2,
                                                            double* Y2, double* Z2,
                                                            double* radius) {
  if (!data::pAIS_InteractiveContext) return;
  assert(id < data::sNumberOfArtefacts && "CADMBTB_buildArtefactLine id out of range");
  if (data::sAISArtefacts[id]) {
    data::pAIS_InteractiveContext->Erase(data::sAISArtefacts[id], true);
    data::sAISArtefacts[id] = nullptr;
  }
  if (!X1) return;
  gp_Pnt P1;
  P1.SetCoord(*X1, *Y1, *Z1);
  gp_Pnt P2;
  P2.SetCoord(*X2, *Y2, *Z2);
  gp_Dir V;
  V.SetCoord(-(*X1 - *X2), -(*Y1 - *Y2), -(*Z1 - *Z2));
  gp_Ax2 gpA2(P1, V);
  double l =
      sqrt((*X1 - *X2) * (*X1 - *X2) + (*Y1 - *Y2) * (*Y1 - *Y2) + (*Z1 - *Z2) * (*Z1 - *Z2));
  BRepPrim_Cylinder makeCyl(gpA2, *radius, l);
  gp_Ax2 gpA2b(P2, V);
  BRepPrim_Cone makeCone(gpA2b, 5 * (*radius), 0, 0.1 * l);
  TopoDS_Compound compound;
  BRep_Builder B;
  B.MakeCompound(compound);
  B.Add(compound, makeCyl.Shell());
  B.Add(compound, makeCone.Shell());
  TopoDS_Vertex aVert1 = BRepBuilderAPI_MakeVertex(P1);
  if (P1.Distance(P2) > data::sMinLineLength) {
    data::sAISArtefacts[id] = new AIS_Shape(compound);
    data::sAISArtefacts[id]->SetColor(Quantity_NOC_WHITE);
    data::sAISArtefacts[id]->SetMaterial(Graphic3d_NOM_PLASTIC);
    // Quantity_Color Qc(Quantity_NOC_BLUE1);
    // data::sAISArtefacts[id]->SetColor(Qc);
    // printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e,
    // P2=%e,%e,%e\n",X1,Y1,Z1,X2,Y2,Z2);
  } else {
    data::sAISArtefacts[id] = new AIS_Shape(aVert1);
    // printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e\n",X1,Y1,Z1);
  }
  data::pAIS_InteractiveContext->Display(data::sAISArtefacts[id], false);
}

void siconos::mechanisms::CADMBTB_buildOrientedLineArtefactLine(unsigned int id, double* X1,
                                                                double* Y1, double* Z1,
                                                                double* X2, double* Y2,
                                                                double* Z2) {
  if (!data::pAIS_InteractiveContext) return;
  assert(id < data::sNumberOfArtefacts && "CADMBTB_buildArtefactLine id out of range");
  if (data::sAISArtefacts[id]) {
    data::pAIS_InteractiveContext->Erase(data::sAISArtefacts[id], true);
    data::sAISArtefacts[id] = nullptr;
  }
  if (!X1) return;
  gp_Pnt P1;
  P1.SetCoord(*X1, *Y1, *Z1);
  gp_Pnt P2;
  P2.SetCoord(*X2, *Y2, *Z2);
  BRepBuilderAPI_MakeEdge MakeEdge(P1, P2);
  TopoDS_Vertex aVert1 = BRepBuilderAPI_MakeVertex(P1);
  if (MakeEdge.IsDone() && P1.Distance(P2) > data::sMinLineLength) {
    TopoDS_Compound compound;
    BRep_Builder B;
    gp_Pnt P3;
    B.MakeCompound(compound);
    B.Add(compound, MakeEdge.Edge());
    double dis = 0.1 * P1.Distance(P2);
    double P1P2x = *X2 - *X1;
    double P1P2y = *Y2 - *Y1;
    double P1P2z = *Z2 - *Z1;
    if (fabs(P1P2x) > fabs(P1P2y)) {
      if (fabs(P1P2x) > fabs(P1P2z)) {
        P3.SetCoord(*X2, *Y2 + dis, *Z2 + dis);
      } else {
        P3.SetCoord(*X2 + dis, *Y2 + dis, *Z2 + 0);
      }
    } else {
      if (fabs(P1P2y) > fabs(P1P2z)) {
        P3.SetCoord(*X2 + dis, *Y2 + 0, *Z2 + dis);
      } else {
        P3.SetCoord(*X2 + dis, *Y2 + dis, *Z2 + 0);
      }
    }
    BRepBuilderAPI_MakeEdge MakeEdge2(P2, P3);
    B.MakeCompound(compound);
    B.Add(compound, MakeEdge.Edge());
    if (P2.Distance(P3) > data::sMinLineLength) B.Add(compound, MakeEdge2.Edge());

    data::sAISArtefacts[id] = new AIS_Shape(compound);
    // printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e,
    // P2=%e,%e,%e\n",X1,Y1,Z1,X2,Y2,Z2);

    // /** V.A. Attempt yo draw arrow with OpenCascade. It Remains to do :
    //  *  + Fix the solid display of the arrow and the length
    //  *  + Manage a list a Presentation  objects, the display and the delete
    //  at each time--step.
    //  */
    //  Handle(Prs3d_Presentation) aPresentation = new Prs3d_Presentation
    //  (data::pAIS_InteractiveContext->CurrentViewer()->Viewer());

    // // // void Prs3d_Arrow::Draw	(	const
    // Handle(Prs3d_Presentation)& 	aPresentation,
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

  } else {
    data::sAISArtefacts[id] = new AIS_Shape(aVert1);
    // printf("CADMBTB_buildArtefactLine P1 = %e, %e, %e\n",X1,Y1,Z1);
  }

  data::pAIS_InteractiveContext->Display(data::sAISArtefacts[id], true);
}

void siconos::mechanisms::CADMBTB_setNbOfArtefacts(unsigned int nb) {
  for (int ii = 0; ii < NB_OBJ; ii++) data::sAISArtefacts[ii] = nullptr;

  assert(nb < NB_OBJ);
  data::sNumberOfArtefacts = nb;
}

// Note FP: the following function is not declared in any header file
// and not used anywhere --> comment. Remove later ?
// void siconos::mechanisms::CADMBTB_setContactAISdParam(unsigned int IdParam,
//                                                       unsigned int idContact,
//                                                       unsigned int idShape,
//                                                       double& v) {
//   assert((int)idContact < data::sNumberOfContacts &&
//          "CADMBTB_setContactAISdParam contactId out of range");

//   unsigned int idShape1 = data::sNumberOfObj +
//                           (2 * idContact - 2 * data::sNumberOfContacts) +
//                           idShape;

//   switch (IdParam) {
//     case 0:
//       data::spAISTrans[idShape1] = v;
//       break;
//     default:
//       printf("Error:  CADMBTB_setAISdParam IdParam out of range.");
//   }
// }

TopoDS_Shape siconos::mechanisms::CADMBTB_TopoDS(unsigned int numDS) {
  return data::sTopoDS[numDS];
}

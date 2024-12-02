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

#include "CADMBTB_PYTHON_API.hpp"

#include <assert.h>
#include <stdio.h>

#include <AIS_InteractiveContext.hxx>
#include <AIS_Shape.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <V3d_View.hxx>

#include "CADMBTB_DATA.hpp"

AIS_InteractiveContext* siconos::mechanisms::data::pAIS_InteractiveContext = nullptr;
V3d_View* siconos::mechanisms::data::pV3d_View = nullptr;
int siconos::mechanisms::data::sCmpDump = 0;
unsigned int siconos::mechanisms::data::sCmpDumpMan = 0;

void siconos::mechanisms::CADMBTB_loadArtefactCADFile(const char* fileName, double trans) {
  if (!data::pAIS_InteractiveContext) return;

  printf("CADMBTB_loadArtefactCADFile using file %s.\n", fileName);
  STEPControl_Reader aReader;
  IFSelect_ReturnStatus status = aReader.ReadFile(fileName);
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
          AIS_Shape* pAis = new AIS_Shape(aReader.Shape(i));
          pAis->SetTransparency(trans);
          data::pAIS_InteractiveContext->Display(pAis, true);
        }
      }
    }
  }
}

void siconos::mechanisms::CADMBTB_setGraphicContext(AIS_InteractiveContext& aisContext) {
  data::pAIS_InteractiveContext = &aisContext;
}

void siconos::mechanisms::CADMBTB_setGraphicView(V3d_View& aView) { data::pV3d_View = &aView; }
void siconos::mechanisms::CADMBTB_disableGraphic() {
  data::pAIS_InteractiveContext = nullptr;
  data::pV3d_View = nullptr;
}

void siconos::mechanisms::CADMBTB_setContactDParam(unsigned int IdParam,
                                                   unsigned int idContact,
                                                   unsigned int idShape, double v) {
  assert((int)idContact < data::sNumberOfContacts &&
         "CADMBTB_setContactAISdParam contactId out of range");
  unsigned int idShape1 =
      data::sNumberOfObj + (2 * idContact - 2 * data::sNumberOfContacts) + idShape;

  CADMBTB_setShapeDParam(IdParam, idShape1, v);
}

void siconos::mechanisms::CADMBTB_setShapeDParam(unsigned int IdParam, unsigned int idShape,
                                                 double v) {
  assert(idShape < data::sNumberOfObj && "CADMBTB_setShapeDParam idShape out of range");
  switch (IdParam) {
    case 0:
      data::spAISTrans[idShape] = v;
      break;
    default:
      printf("Error:  CADMBTB_setShapeDParam IdParam out of range.");
  }
}
void siconos::mechanisms::CADMBTB_setIParam(unsigned int IdParam, int v) {
  if (v)
    data::sDumpGraphic = 1;
  else
    data::sDumpGraphic = 0;
}
void siconos::mechanisms::CADMBTB_DumpGraphic() {
  if (!data::pAIS_InteractiveContext) return;
  //  data::pAIS_InteractiveContext->UpdateCurrentViewer() ;
  if (data::pV3d_View) {
    char file[16];
    snprintf(file, 16, "manual%d.jpg", data::sCmpDumpMan);
    data::pV3d_View->Dump(file);
    data::sCmpDumpMan++;
  }
}
void siconos::mechanisms::CADMBTB_print_dist(unsigned int v) { data::sCADPrintDist = v; }

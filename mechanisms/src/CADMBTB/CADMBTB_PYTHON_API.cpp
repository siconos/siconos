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

#include "CADMBTB_DATA.hpp"
AIS_InteractiveContext* pAIS_InteractiveContext = nullptr;
V3d_View* pV3d_View = nullptr;
int sCmpDump = 0;
unsigned int sCmpDumpMan = 0;
void CADMBTB_setGraphicContext(AIS_InteractiveContext& aisContext) {
  pAIS_InteractiveContext = &aisContext;
}

void CADMBTB_setGraphicView(V3d_View& aView) { pV3d_View = &aView; }
void CADMBTB_disableGraphic() {
  pAIS_InteractiveContext = nullptr;
  pV3d_View = nullptr;
}

void CADMBTB_setContactDParam(unsigned int IdParam, unsigned int idContact,
                              unsigned int idShape, double v) {
  assert((int)idContact < sNumberOfContacts &&
         "CADMBTB_setContactAISdParam contactId out of range");
  unsigned int idShape1 = sNumberOfObj + (2 * idContact - 2 * sNumberOfContacts) + idShape;

  CADMBTB_setShapeDParam(IdParam, idShape1, v);
}

void CADMBTB_setShapeDParam(unsigned int IdParam, unsigned int idShape, double v) {
  assert(idShape < sNumberOfObj && "CADMBTB_setShapeDParam idShape out of range");
  switch (IdParam) {
    case 0:
      spAISTrans[idShape] = v;
      break;
    default:
      printf("Error:  CADMBTB_setShapeDParam IdParam out of range.");
  }
}
void CADMBTB_setIParam(unsigned int IdParam, int v) {
  if (v)
    sDumpGraphic = 1;
  else
    sDumpGraphic = 0;
}
void CADMBTB_DumpGraphic() {
  if (!pAIS_InteractiveContext) return;
  //  pAIS_InteractiveContext->UpdateCurrentViewer() ;
  if (pV3d_View) {
    char file[16];
    snprintf(file, 16, "manual%d.jpg", sCmpDumpMan);
    pV3d_View->Dump(file);
    sCmpDumpMan++;
  }
}
void CADMBTB_print_dist(unsigned int v) { sCADPrintDist = v; }

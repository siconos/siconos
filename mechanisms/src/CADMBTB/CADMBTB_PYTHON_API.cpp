#include "CADMBTB_DATA.hpp"
#include "CADMBTB_PYTHON_API.hpp"
#include <assert.h>

AIS_InteractiveContext * pAIS_InteractiveContext=0;
V3d_View * pV3d_View=0;
int sCmpDump=0;
unsigned int sCmpDumpMan=0;
void CADMBTB_setGraphicContext(AIS_InteractiveContext & aisContext)
{
  pAIS_InteractiveContext = & aisContext;
}

void CADMBTB_setGraphicView(V3d_View & aView)
{
  pV3d_View= & aView;
}
void CADMBTB_disableGraphic()
{
  pAIS_InteractiveContext =0;
  pV3d_View=0;
}

void CADMBTB_setContactDParam(unsigned int IdParam,unsigned int idContact,unsigned int idShape,double  v)
{
  assert(idContact<sNumberOfContacts && "CADMBTB_setContactAISdParam contactId out of range");
  unsigned int idShape1=sNumberOfObj+(2*idContact-2*sNumberOfContacts)+idShape;

  CADMBTB_setShapeDParam(IdParam,idShape1,v);
}

void CADMBTB_setShapeDParam(unsigned int IdParam,unsigned int idShape,double  v)
{
  assert(idShape<sNumberOfObj && "CADMBTB_setShapeDParam idShape out of range");
  switch(IdParam)
  {
  case 0:
    spAISTrans[idShape]=v;
    break;
  default:
    printf("Error:  CADMBTB_setShapeDParam IdParam out of range.");
  }
}
void CADMBTB_setIParam(unsigned int IdParam, int v)
{
  if(v)
    sDumpGraphic=1;
  else
    sDumpGraphic=0;
}
void CADMBTB_DumpGraphic()
{
  if(!pAIS_InteractiveContext)
    return;
  //  pAIS_InteractiveContext->UpdateCurrentViewer() ;
  if(pV3d_View)
  {
    char file[16];
    sprintf(file,"manual%d.gif",sCmpDumpMan);
    pV3d_View->Dump(file);
    sCmpDumpMan++;
  }

}
void CADMBTB_print_dist(unsigned int v)
{
  sCADPrintDist=v;
}

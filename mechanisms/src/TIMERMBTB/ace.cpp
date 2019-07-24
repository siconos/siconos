#include "ace.h"
#include <math.h>
#include <cstdlib>

aceTime ACE_times[ACE_TIMER_LAST];


using namespace std;




void ACE_INIT_TIME()
{
  ACE_times[ACE_TIMER_MAIN].setName("main ");
  ACE_times[ACE_TIMER_GRAPHIC].setName("OCC Graphic Update ");
  ACE_times[ACE_TIMER_DIST].setName("OCC Distance computation ");
  ACE_times[ACE_TIMER_SICONOS].setName("Siconos  ");
  ACE_times[ACE_TIMER_UPDATE_POS].setName("OCC update position ");
  ACE_times[ACE_TIMER_UV_TO_3D].setName("OCC UV to 3D ");
  ACE_times[ACE_TIMER_UV_POLYNOME].setName("UV polynome evaluation ");
  ACE_times[ACE_TIMER_UV_CLASS].setName("OCC Classifier ");
  ACE_times[ACE_TIMER_UV_GRAD].setName("OCC UV to Grad ");
  ACE_times[ACE_TIMER_3D_PROJ].setName("OCC 3D proj ");
  ACE_times[ACE_TIMER_CAD_VALUE].setName("CAD value ");
  ACE_times[ACE_TIMER_CAD_VALUES].setName("CAD values ");
  ACE_times[ACE_TIMER_CAD_1].setName("CAD 1 ");
  ACE_times[ACE_TIMER_CAD_12].setName("CAD 12 ");
  ACE_times[ACE_TIMER_CAD_13].setName("CAD 13 ");
  ACE_times[ACE_TIMER_CAD_14].setName("CAD 14 ");
  ACE_times[ACE_TIMER_CAD_15].setName("CAD 15 ");
  ACE_times[ACE_TIMER_CAD_16].setName("CAD 16 ");
  ACE_times[ACE_TIMER_CAD_17].setName("CAD 17 ");
  ACE_times[ACE_TIMER_CAD_OK].setName("CAD OK ");
}
void ACE_PRINT_TIME()
{
  for(int i=0; i <ACE_TIMER_LAST; i++)
    ACE_times[i].print();
}

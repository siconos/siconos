/************************************************************************
  			acetime.cpp
**************************************************************************/


#include "acetime.h"
#include <string.h>


#define ACE_WITH_TIMER

aceTime::~aceTime()
{
}
aceTime::aceTime()
{
#ifdef ACE_WITH_TIMER
  mIsRunning=false;
  mCumul=0;
  mCall=0;
  strcpy(mName,"NoName");
#endif
}

void aceTime::start()
{
#ifdef ACE_WITH_TIMER
  gettimeofday(&mStart,NULL);
  mCall++;
  mIsRunning=true;
#endif
}

void aceTime::stop()
{
#ifdef ACE_WITH_TIMER
  if(mIsRunning)
  {
    timeval aux;
    gettimeofday(&aux,NULL);
    mCumul += (aux.tv_sec - mStart.tv_sec)*1000000 +(aux.tv_usec - mStart.tv_usec) ;
    mIsRunning = false;
  }
#endif
}
void aceTime::print(ostream& os)
{
#ifdef ACE_WITH_TIMER
  os << mName<<" : "<<mCumul<<" usec. "<<mCall<<" calls.";
  if(mCall)
    os << "ie "<<mCumul/mCall<<" per call.";
  os <<endl;
#endif
}
void aceTime::setName(char const *Name)
{
#ifdef ACE_WITH_TIMER
  strcpy(mName,Name);
#endif
}



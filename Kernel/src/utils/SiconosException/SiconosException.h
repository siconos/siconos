//$Id: SiconosException.h,v 1.1 2005/01/05 15:59:04 charlety Exp $
/** \class SiconosException
*   \brief This class represent all the exeption in the Siconos platform
*   \author A. Ravoux
*   \version 1.0
*   \date (Creation) 05/25/2004
*
*$Date: 2005/01/05 15:59:04 $
* $Revision: 1.1 $
* $Author: charlety $
* $Source: /CVS/Siconos/SICONOS/src/utils/SiconosException/SiconosException.h,v $
*
* SiconosException should not be throws directly; prefer to use an inherit class
* This exception can be catched by "catch(SiconosException)"
*
*/

#ifndef __SiconosException__
#define __SiconosException__

#include <string>

using namespace std;

// --------------------------------------------------------------------------
class SiconosException
{
public:

  /**
   * \fn SiconosException()
   * \brief constructor
   */
  SiconosException();

  /**
   * \fn SiconosException(string report)
   * \brief constructor with a report
   * \param string report : exception description
   */
  SiconosException(string report);

  /**
   * \fn ~SiconosException()
   * \brief destructor
   */
  ~SiconosException();

  /**
   * \fn string report()
   * \brief return the report of the exception
   * \return string report : exception description
   */
  inline string report()
  {
    return this->reportMsg;
  }


protected:
  /** report message which describe the exception */
  string reportMsg;

};

#endif //__SiconosException__

//$Id: SiconosVectorException.h,v 1.1 2004/09/10 10:00:22 charlety Exp $
/** \class SiconosVectorException
*   \brief This class represent an exeption causing by a SiconosMatrix
*   \author A. Ravoux
*   \version 1.0
*   \date (Creation) 05/25/2004
*
* $Date: 2004/09/10 10:00:22 $
* $Revision: 1.1 $
* $Author: charlety $
* $Source: /CVS/Siconos/SICONOS/src/utils/NewSiconosVector/SiconosVectorException.h,v $
*
*
* SiconosVectorException must be throws when an error is find in a SiconosVector
* This exception can be catched by "catch(SiconosVectorException)" or "catch(SiconosException)"
*
*/

#ifndef __SiconosVectorException__
#define __SiconosVectorException__

#include "SiconosException.h"

// --------------------------------------------------------------------------
class SiconosVectorException : public SiconosException
{
public:

  /**
   * \fn SiconosVectorException()
   * \brief constructor
   */
  SiconosVectorException();

  /**
   * \fn SiconosVectorException(string report)
   * \brief constructor with a report
   * \param string report : exception description
   */
  SiconosVectorException(string report);

  /**
   * \fn ~SiconosVectorException()
   * \brief destructor
   */
  ~SiconosVectorException();

  /**
   * \fn static void selfThrow()
   * \brief static function which throw a SiconosVectorException
   * \exception SiconosVectorException
   */
  static void selfThrow();

  /**
   * \fn static void selfThrow(string report)
   * \brief static function which throw a SiconosVectorException with a report
   * \param string report : exception description
   * \exception SiconosVectorException
   */
  static void selfThrow(string report);

};

#endif //__SiconosVectorException__

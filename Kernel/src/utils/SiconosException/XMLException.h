//$Id: XMLException.h,v 1.1 2005/01/05 15:59:04 charlety Exp $
/** \class XMLException
*   \brief This class represent an exeption causing by an XML class of the platform
*   \author A. Ravoux
*   \version 1.0
*   \date (Creation) 05/25/2004
*
* $Date: 2005/01/05 15:59:04 $
* $Revision: 1.1 $
* $Author: charlety $
* $Source: /CVS/Siconos/SICONOS/src/utils/SiconosException/XMLException.h,v $
*
*
* XMLException must be throws when an error is find in an XML class
* This exception can be catched by "catch(XMLException)" or "catch(SiconosException)"
*
*
*/

#ifndef __XMLException__
#define __XMLException__

#include "SiconosException.h"

// --------------------------------------------------------------------------
class XMLException: public SiconosException
{
public:

  /**
   * \fn XMLException()
   * \brief constructor
   */
  XMLException();

  /**
   * \fn XMLException(string report)
   * \brief constructor with a report
   * \param string report : exception description
   */
  XMLException(string report);

  /**
   * \fn ~XMLException()
   * \brief destructor
   */
  ~XMLException();

  /**
   * \fn static void selfThrow()
   * \brief static function which throw a XMLException
   * \exception XMLException
   */
  static void selfThrow();

  /**
   * \fn static void selfThrow(string report)
   * \brief static function which throw a XMLException with a report
   * \param string report : exception description
   * \exception XMLException
   */
  static void selfThrow(string report);

};

#endif //__XMLException__

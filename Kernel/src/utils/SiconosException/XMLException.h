/** \class XMLException
*   \brief This class represent an exeption causing by an XML class of the platform
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date (Creation) 05/25/2004
*
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
#include <iostream>

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
   * \fn XMLException(const std::string& report)
   * \brief constructor with a report
   * \param std::string report : exception description
   */
  XMLException(const std::string& report);

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
  static void selfThrow() ;

  /**
   * \fn static void selfThrow(const std::string& report)
   * \brief static function which throw a XMLException with a report
   * \param std::string report : exception description
   * \exception XMLException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__XMLException__

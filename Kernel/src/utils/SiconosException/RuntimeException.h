/** \class RuntimeException
*   \brief This class represent a runtime exeption causing by the plateforme
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date (Creation) 05/25/2004
*
*
* RuntimeException can be throws for example when a pointer is used but not allocated
* This exception can be catched by "catch(RuntimeException)" or "catch(SiconosException)"
*
*/

#ifndef __RuntimeException__
#define __RuntimeException__

#include "SiconosException.h"

class RuntimeException: public SiconosException
{
public:

  /**
   * \fn RuntimeException()
   * \brief constructor
   */
  RuntimeException();

  /**
   * \fn RuntimeException(const string& report)
   * \brief constructor with a report
   * \param string report : exception description
   */
  RuntimeException(const std::string& report);

  /**
   * \fn ~RuntimeException()
   * \brief destructor
   */
  ~RuntimeException();

  /**
   * \fn static void selfThrow()
   * \brief static function which throw a RuntimeException
   * \exception RuntimeException
   */
  static void selfThrow() ;

  /**
   * \fn static void selfThrow(string report)
   * \brief static function which throw a RuntimeException with a report
   * \param string report : exception description
   * \exception RuntimeException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__RuntimeException__

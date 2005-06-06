/** \class SiconosMatrixException
*   \brief This class represent an exeption causing by a SiconosMatrix
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date (Creation) 05/25/2004
*
*
*
* SiconosMatrixException must be throws when an error is find in a SiconosMatrix
* This exception can be catched by "catch(SiconosMatrixException)" or "catch(SiconosException)"
*
*/

#ifndef __SiconosMatrixException__
#define __SiconosMatrixException__

#include "SiconosException.h"

// --------------------------------------------------------------------------
class SiconosMatrixException : public SiconosException
{
public:

  /**
   * \fn SiconosMatrixException()
   * \brief constructor
   */
  SiconosMatrixException();

  /**
   * \fn SiconosMatrixException(const std::string& report)
   * \brief constructor with a report
   * \param std::string report : exception description
   */
  SiconosMatrixException(const std::string& report);

  /**
   * \fn ~SiconosMatrixException()
   * \brief destructor
   */
  ~SiconosMatrixException();

  /**
   * \fn static void selfThrow()
   * \brief static function which throw a SiconosMatrixException
   * \exception SiconosMatrixException
   */
  static void selfThrow()  ;

  /**
   * \fn static void selfThrow(const std::string& report)
   * \brief static function which throw a SiconosMatrixException with a report
   * \param std::string report : exception description
   * \exception SiconosMatrixException
   */
  static void selfThrow(const std::string& report) ;


};

#endif //__SiconosMatrixException__

/** \class SiconosVectorException
*   \brief This class represent an exeption causing by a SiconosMatrix
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date (Creation) 05/25/2004
*
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
   * \fn SiconosVectorException(const std::string& report)
   * \brief constructor with a report
   * \param std::string report : exception description
   */
  SiconosVectorException(const std::string& report);

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
  static void selfThrow()  ;

  /**
   * \fn static void selfThrow(const std::string& report)
   * \brief static function which throw a SiconosVectorException with a report
   * \param std::string report : exception description
   * \exception SiconosVectorException
   */
  static void selfThrow(const std::string& report) ;

};

#endif //__SiconosVectorException__

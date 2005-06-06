/** \class SiconosException
*   \brief This class represent all the exeption in the Siconos platform
*  \author SICONOS Development Team - copyright INRIA
*   \version 1.0
*   \date (Creation) 05/25/2004
*
*
* SiconosException should not be throws directly; prefer to use an inherit class
* This exception can be catched by "catch(SiconosException)"
*
*/

#ifndef __SiconosException__
#define __SiconosException__

#include <string>

class SiconosException
{
public:

  /**
   * \fn SiconosException()
   * \brief constructor
   */
  SiconosException();

  /**
   * \fn SiconosException(const string&)
   * \brief constructor with a report
   * \param string report : exception description
   */
  SiconosException(const std::string&);

  /**
   * \fn ~SiconosException()
   * \brief destructor
   */
  virtual ~SiconosException();

  /**
   * \fn string report()
   * \brief return the report of the exception
   * \return string report : exception description
   */
  inline std::string report() const
  {
    return reportMsg;
  } ;


protected:
  /** report message which describe the exception */
  std::string reportMsg;
};

#endif //__SiconosException__

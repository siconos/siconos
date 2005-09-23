#ifndef SICONOSSHAREDLIBRARYEXCEPTION_H
#define SICONOSSHAREDLIBRARYEXCEPTION_H

#include "SiconosException.h"

// --------------------------------------------------------------------------
class SiconosSharedLibraryException : public SiconosException
{
public:

  /**
   * constructor
   */
  SiconosSharedLibraryException();

  /**
   * constructor
   * @param string which describe the exception
   */
  SiconosSharedLibraryException(const std::string& report);

  /**
   * destructor
   */
  ~SiconosSharedLibraryException();

  static void selfThrow() ;

  static void selfThrow(const std::string& report) ;

};

#endif //SICONOSSHAREDLIBRARYEXCEPTION_H

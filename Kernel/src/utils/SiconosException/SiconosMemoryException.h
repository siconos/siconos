#ifndef SICONOSMEMORYEXCEPTION_H
#define SICONOSMEMORYEXCEPTION_H

#include "SiconosException.h"

// --------------------------------------------------------------------------
class SiconosMemoryException : public SiconosException
{
public:

  /**
   * constructor
   */
  SiconosMemoryException();

  /**
   * constructor
   * @param std::string which describe the exception
   */
  SiconosMemoryException(const std::string& report);

  /**
   * destructor
   */
  ~SiconosMemoryException();

  static void selfThrow() ;

  static void selfThrow(const std::string& report) ;

};

#endif //SICONOSMEMORYEXCEPTION_H

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
   * @param string which describe the exception
   */
  SiconosMemoryException(string report);

  /**
   * destructor
   */
  ~SiconosMemoryException();

  static void selfThrow();

  static void selfThrow(string report);

};

#endif //SICONOSMEMORYEXCEPTION_H

/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2021 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


/** ! \file SiconosException.hpp
    \brief Tools to deal with exceptions in Siconos.


*/

#ifndef __SICONOSEXCEPTION_HPP__
#define __SICONOSEXCEPTION_HPP__

#include <exception>
#include <boost/throw_exception.hpp>
#include <boost/exception/diagnostic_information.hpp>
#include <cstring>
#include <iostream>
#include <errno.h>



namespace Siconos
{


  /** Siconos generic exception

      Usage :

      1/ to throw an exception, use :
      THROW_EXCEPTION("some message")
      // or
      THROW_EXCEPTION()

      2/ to catch it :

      try{
      ... // call to functions that may throw exceptions
      }
      catch(Siconos::exception::Exception& e)
      {
      Siconos::exception::process(e);
      }
      catch(...)
      {
      Siconos::exception::process();
      }

      The outputs are :
      * the file, position in the file (line number) and function name
      that thrown the exception,
      * the errno code (0 if unknown),
      * the corresponding description of errno code (if any),
      * the message added when the exception has been thrown.

      */
  class exception : public virtual std::exception, public virtual boost::exception
  {

  public:
    // Boost typedef to handle exception extra outputs using << operator.
    // See https://www.boost.org/doc/libs/1_74_0/libs/exception/doc/error_info.html
    // About boost exceptions, see https://www.boost.org/doc/libs/1_74_0/libs/exception/doc/error_info.html
    // About errno handling, see https://en.cppreference.com/w/cpp/error/errno
    typedef boost::error_info<struct tag_errno_code,int> errno_code;
    typedef boost::error_info<struct tag_errno_description,const char *> errno_description;
    typedef boost::error_info<struct message, std::string> extra_message;
    /** Outputs diagnostic information about unhandled exceptions.
        Must be called inside a catch section.
    */
    static inline void process()
    {
      std::cerr << boost::current_exception_diagnostic_information(true) << std::endl;
    }

    /**
       Outputs diagnostic information about exceptions.
       Must be called inside a catch section.
    */
    static inline void process(Siconos::exception& e)
    {
      std::cerr << boost::diagnostic_information(e, true) << std::endl;
    }

  };

}//namespace Siconos

/** Wrap exception throwing inside Siconos. */
#define THROW_EXCEPTION(X) BOOST_THROW_EXCEPTION(Siconos::exception() << Siconos::exception::extra_message(X) << Siconos::exception::errno_code(errno) <<Siconos::exception::errno_description(std::strerror(errno)));

#endif

/* Siconos-Kernel, Copyright INRIA 2005-2012.
* Siconos is a program dedicated to modeling, simulation and control
* of non smooth dynamical systems.
* Siconos is a free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
* Siconos is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with Siconos; if not, write to the Free Software
* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*
* Contact: Vincent ACARY, siconos-team@lists.gforge.inria.fr
*/

/*! \file TypeName.hpp
  \brief get a string name from visitable classes 
*/

#include <string>
namespace Type
{
#undef REGISTER
#define REGISTER(X) case Type:: X : r.reset(new std::string(#X)); break;

#undef REGISTER_STRUCT
#define REGISTER_STRUCT(X) REGISTER(X)
#undef REGISTER_BASE
#undef REGISTER_BASE_EXTERN
#define REGISTER_BASE(X,Y) REGISTER(X)

#define REGISTER_BASE_EXTERN(X,Y) REGISTER_BASE(X,Y)

inline std11::shared_ptr<std::string> str(const Siconos& X)
{
  std11::shared_ptr<std::string> r;

  switch (X)
  {
    SICONOS_VISITABLES()
  default:
    assert(0);
  }

  return(r);
}


template <class C>
std::string name(const C& c)
{
  return *(Type::str(Type::value(c)));
}

}

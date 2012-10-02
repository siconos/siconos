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

/*! \file Cmp.hpp
Classes related to object ordering in SiconosSet.
*/

#ifndef CMP_H
#define CMP_H

#if (__cplusplus >= 201103L) && !defined(USE_BOOST_FOR_CXX11)
#include <memory>
namespace std11 = std;
#else
#include <boost/shared_ptr.hpp>
namespace std11 = boost;
#endif

/**  Virtual functors class
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) April 25, 2006
 *
 *  Note: this is strongly inspired from tutorial http://www.newty.de/fpt/functor.html
 *
 */
template <class U> class myTFunctor
{
public:

  /** destructor
   */
  virtual ~myTFunctor() {};

  /**
   *  \return int
   */
  virtual U Call() = 0;      // call using function
};

/** derived template class for functors
*  \author SICONOS Development Team - copyright INRIA
*  \version 3.0.0.
*  \date (Creation) April 25, 2006
*
*  Note: this is strongly inspired from tutorial http://www.newty.de/fpt/functor.html
*/
template <class TClass, class U> class myTSpecificFunctor : public myTFunctor<U>
{
private:
  /** pointer to member function */
  U(TClass::*fpt)() const;
  /** pointer to object */
  std11::shared_ptr<TClass> pt2Object;

public:

  /**
   *  constructor - takes pointer to an object and pointer to a member and stores
   * them in two private variables
   * \param pointer to TClass
   * \param pointer to function of type const int(TClass::*_fpt)() const
   */
  myTSpecificFunctor(std11::shared_ptr<TClass> _pt2Object, U(TClass::*_fpt)() const)
  {
    pt2Object = _pt2Object;
    fpt = _fpt;
  };

  /** override function "Call" that executes member function
  *  \return int
  */
  virtual U Call()
  {
    return (*pt2Object.*fpt)();
  };
};

/** Cmp<T,U> objects are used to provide "strict weak ordering" relations,
 * used in SiconosSet to compare objects; See STL doc or
 * SiconosSet.hpp class for example of use.
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 3.0.0.
 *  \date (Creation) April 25, 2006
 *
 *  \code
 *  Cmp<T,U> myCmp;
 *  // a and b objects of type T
 *  //
 *  bool res = myCmp(a,b); // return true if "a<b" else false.
 *
 * \endcode
 *
 * "a<b" is evaluated according to the returning value of the member \n
 * function pointed by fpt.
 *
 * Cmp objects are not supposed to be used directly but as arguments in SiconosSet.hpp.
 * See that class for details.
 *
 *  Note: this is strongly inspired from http://www.josuttis.com/libbook/
 */
template <class T, class U> class Cmp
{
public:

  /** sorting criterion list*/
  enum cmp_mode {normal, reverse};

private:

  /** sorting criterion*/
  cmp_mode mode;

  /**  pointer to member function: return the value that will be used for sorting */
  U(T::*fpt)() const ;

public:

  /** default constructor */
  Cmp() {};

  /**  constructor for sorting criterion
  *  default criterion uses value normal
  * \param a pointer to function of type const int(T::*_fpt)() const
  */
  Cmp(U(T::*_fpt)() const, cmp_mode m = normal): mode(m), fpt(_fpt) {};

  /** comparison of elements
  * \param two pointers to T
  */
  bool operator()(std11::shared_ptr<T> t1, std11::shared_ptr<T> t2) const
  {
    // set functors and call pointers to function to get values to be compared
    myTSpecificFunctor<T, U> specFuncA(t1, fpt);
    myTSpecificFunctor<T, U> specFuncB(t2, fpt);
    myTFunctor<U>* vTable[] = { &specFuncA, &specFuncB};
    U i1 = vTable[0]->Call();
    U i2 = vTable[1]->Call();

    return mode == normal ? i1 < i2 : i2 < i1;
  }
  /** comparison of sorting criteria
  * \param a Cmp
  * \return a bool
  */
  bool operator== (const Cmp& rc)
  {
    return mode == rc.mode;
  }
};
#endif

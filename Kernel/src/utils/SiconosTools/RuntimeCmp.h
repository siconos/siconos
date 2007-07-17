/* Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
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
 * Contact: Vincent ACARY vincent.acary@inrialpes.fr
*/

/*! \file RuntimeCmp.h

*/

#ifndef RUNTIMECMP_H
#define RUNTIMECMP_H

/**  Virtual functors class
 *
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 2.1.1.
 *  \date (Creation) April 25, 2006
 *
 *  Note: this is strongly inspired from tutorial http://www.newty.de/fpt/functor.html
 */
class TFunctor
{
public:

  /** destructor
   */
  virtual ~TFunctor() {};

  /**
   *  \return unsigned long int
   */
  virtual const unsigned long int Call() = 0;      // call using function
};

/** derived template class for functors
*  \author SICONOS Development Team - copyright INRIA
*  \version 2.1.1.
*  \date (Creation) April 25, 2006
*
*  Note: this is strongly inspired from tutorial http://www.newty.de/fpt/functor.html
*/
template <class TClass> class TSpecificFunctor : public TFunctor
{
private:
  /** pointer to member function */
  const unsigned long int (TClass::*fpt)() const;
  /** pointer to object */
  TClass* pt2Object;

public:

  /**
   *  constructor - takes pointer to an object and pointer to a member and stores
   * them in two private variables
   * \param pointer to TClass
   * \param pointer to function of type const unsigned long int(TClass::*_fpt)() const
   */
  TSpecificFunctor(TClass* _pt2Object, const unsigned long int(TClass::*_fpt)() const)
  {
    pt2Object = _pt2Object;
    fpt = _fpt;
  };

  /** override function "Call" that executes member function
  *  \return unsigned long int
  */
  virtual const unsigned long int Call()
  {
    return (*pt2Object.*fpt)();
  };
};

/** Template to provide comparison operator in stl set or map
*   see examples of using in EventsManager.h
*  \author SICONOS Development Team - copyright INRIA
*  \version 2.1.1.
*  \date (Creation) April 25, 2006
*
* Convention: objects to be compared must have a function of type:
*  const unsigned long int name() const (whatever name is)
* The comparison will be based on the return value of this function.
*
*  Note: this is strongly inspired from http://www.josuttis.com/libbook/
*/
template <class T> class RuntimeCmp
{
public:

  /** sorting criterion list*/
  enum cmp_mode {normal, reverse};

private:

  /** sorting criterion*/
  cmp_mode mode;

  /**  pointer to member function: return the value that will be used for sorting */
  const unsigned long int (T::*fpt)() const ;

public:

  /**  constructor for sorting criterion
  *  default criterion uses value normal
  * \param a pointer to function of type const unsigned long int(T::*_fpt)() const
  */
  RuntimeCmp(const unsigned long int(T::*_fpt)() const, cmp_mode m = normal): mode(m), fpt(_fpt) {};

  /** comparison of elements
  * \param two pointers to T
  */
  bool operator()(T* t1, T* t2) const
  {
    // set functors and call pointers to function to get values to be compared
    TSpecificFunctor<T> specFuncA(t1, fpt);
    TSpecificFunctor<T> specFuncB(t2, fpt);
    TFunctor* vTable[] = { &specFuncA, &specFuncB};
    unsigned long int i1 = vTable[0]->Call();
    unsigned long int i2 = vTable[1]->Call();

    return mode == normal ? i1 < i2 : i2 < i1;
  }
  /** comparison of sorting criteria
  * \param a RuntimeCmp
  * \return a bool
  */
  bool operator== (const RuntimeCmp& rc)
  {
    return mode == rc.mode;
  }
};
#endif

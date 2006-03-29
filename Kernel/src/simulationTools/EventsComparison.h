/* Siconos-Kernel version 1.1.4, Copyright INRIA 2005-2006.
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

/** \class EventsComparison
 *  \brief class that provides to compare Events and to sort them in the list of EventsManager
 *  \author SICONOS Development Team - copyright INRIA
 *  \version 1.1.4.
 *  \date (Creation) February 23, 2006
 *
 *   - First criteria is time of event.
 *   - Second criteria, type of event (used if two times are equal).
 *     At the time, we set TimeDiscretisationEvent < NonSmoothEvent.
 *
 */

#ifndef EVENTSCOMPARISON_H
#define EVENTSCOMPARISON_H

#include "Event.h"
#include<set>
#include<string>
#include<iostream>

// type for sorting criterion => with template (to be tested with python)
/* template <class T> */
/* class RuntimeCmp { */
/*   public: */
/*     enum cmp_mode {normal, reverse}; */
/*   private: */
/*     cmp_mode mode; */
/*   public:   */
/*     // constructor for sorting criterion */
/*     // - default criterion uses value normal */
/*     RuntimeCmp(cmp_mode m=normal) : mode(m) {} */
/*     // comparison of elements */
/*     bool operator() (const T& t1, const T& t2) const { */
/*         return mode == normal ? t1 < t2 : t2 < t1; */
/*     } */
/*     // comparison of sorting criteria */
/*     bool operator== (const RuntimeCmp& rc) { */
/*         return mode == rc.mode; */
/*     } */
/* }; */

class Event;
// The same without template.

class EventsComparison
{
public:
  enum cmp_mode {normal, reverse};
private:
  cmp_mode mode;
public:
  // constructor for sorting criterion
  // - default criterion uses value normal
  EventsComparison(cmp_mode m = normal) : mode(m) {}
  // comparison of elements
  bool operator()(Event*, Event*) const ;
  // comparison of sorting criteria
  inline bool operator== (const EventsComparison& rc)
  {
    return mode == rc.mode;
  }
};

#endif

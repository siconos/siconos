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
#include "Event.hpp"
#include "TimeDiscretisation.hpp"
#include "RuntimeException.hpp"
#include <cmath>


double Event::tick = DEFAULT_TICK;

Event::Event(double time, int newType): type(newType), dTime(time)
{
  // Initialize and set timeOfEvent.
  mpz_init_set_d(timeOfEvent, rint(time / tick)) ;
}

Event::~Event()
{
  mpz_clear(timeOfEvent);
}

void Event::display() const
{
  std::cout << "===== Event data display =====" <<std::endl;
  std::cout << " - Type: " << type <<std::endl;
  std::cout << " - time (mpz_t format, double format): (";
  mpz_out_str(stdout, 10, timeOfEvent);
  std::cout << ", " << dTime << ")" <<std::endl;
  std::cout << "===== End of Event display =====" <<std::endl;
}

void Event::update()
{
  if (_td)
  {
    _td->increment();
    setTime(_td->currentTime());
  }
}

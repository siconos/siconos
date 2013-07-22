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


double Event::_tick = DEFAULT_TICK;

Event::Event(double time, int newType, bool reschedule):
  _type(newType), _dTime(time), _k(0), _reschedule(reschedule)
{
  // Initialize and set timeOfEvent.
  mpz_init_set_d(_timeOfEvent, rint(time / _tick));
  mpz_init_set_d(_tickIncrement, 0);
}

Event::~Event()
{
  mpz_clear(_timeOfEvent);
  mpz_clear(_tickIncrement);
}

void Event::update(unsigned int)
{
  if (_td) // if no TimeDiscretisation then do nothing
  {
    _k++;
    if (_td->hGmp())
      incrementTime();
    else
      setTime(_td->getTk(_k));
  }
}
void Event::setTimeDiscretisation(SP::TimeDiscretisation td)
{
  _td = td;
  if (_td->hGmp())
  {
    mpf_t tmp;
    mpf_init_set_d(tmp, _tick);
    mpf_div(tmp, *_td->currentTimeStep(), tmp);
    mpz_init_set_ui(_tickIncrement, mpf_get_ui(tmp));
    mpf_clear(tmp);
  }
}

void Event::display() const
{
  std::cout << "===== Event data display =====" <<std::endl;
  std::cout << " - Type: " << _type <<std::endl;
  std::cout << " - time (mpz_t format, double format): (";
  mpz_out_str(stdout, 10, _timeOfEvent);
  std::cout << ", " << _dTime << ")" <<std::endl;
  std::cout << "===== End of Event display =====" <<std::endl;
}

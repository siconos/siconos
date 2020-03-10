/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2018 INRIA.
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
#include "Event.hpp"
#include "TimeDiscretisation.hpp"
#include "RuntimeException.hpp"
#include <cmath>
#include <iostream>


double Event::_tick = DEFAULT_TICK;
bool Event::_eventCreated = false;

Event::Event(double time, int newType, bool reschedule):
  _type(newType), _dTime(time), _k(0), _reschedule(reschedule)
{
  // Initialize and set timeOfEvent.
  mpz_init_set_d(_timeOfEvent, rint(time / _tick));
  mpz_init_set_d(_tickIncrement, 0);
  _eventCreated = true;
}

Event::~Event()
{
  mpz_clear(_timeOfEvent);
  mpz_clear(_tickIncrement);
}

void Event::update(unsigned int k)
{
  if(_td)  // if no TimeDiscretisation then do nothing
  {
    _k++;
    if(_td->hGmp())
      incrementTime();
    else
      setTime(_td->getTk(_k));
  }
}
void Event::setTimeDiscretisation(SP::TimeDiscretisation td)
{
  _td = td;
  if(_td->hGmp())
  {
    mpf_t tmp;
    mpf_init_set_d(tmp, _tick);
    mpf_div(tmp, *_td->currentTimeStep(), tmp);
    mpz_set_ui(_tickIncrement, mpf_get_ui(tmp));
    mpf_clear(tmp);
  }
}

void Event::setTick(double newTick)
{
  if(_eventCreated)
  {
    std::cout << "Warning: you change tick value for EventsManager -> a new initialization of the object is required. " << std::endl;
  }
  _tick = newTick;
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

/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2022 INRIA.
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

#include "SiconosConst.hpp"
#include "TimeDiscretisation.hpp"
#include <iostream>

bool siconos::simulation::Event::_eventCreated = false;

double siconos::simulation::Event::_tick = siconos::internal::DEFAULT_TICK;

siconos::simulation::Event::Event(double time, EventType newType, bool reschedule)
    : _type(newType), _dTime(time), _reschedule(reschedule)
{
  // Initialize and set timeOfEvent.
  mpz_init_set_d(_timeOfEvent, rint(time / _tick));
  mpz_init_set_d(_tickIncrement, 0);
  _eventCreated = true;
}

siconos::simulation::Event::~Event() noexcept
{
  mpz_clear(_timeOfEvent);
  mpz_clear(_tickIncrement);
}

void siconos::simulation::Event::update(unsigned int k)
{
  if (_td)  // if no TimeDiscretisation then do nothing
  {
    _k++;
    if (_td->hGmp())
      incrementTime();
    else
      setTime(_td->getTk(_k));
  }
}
void siconos::simulation::Event::setTimeDiscretisation(std::shared_ptr<TimeDiscretisation> td)
{
  _td = td;
  if (_td->hGmp()) {
    mpf_t tmp;
    mpf_init_set_d(tmp, _tick);
    mpf_div(tmp, *_td->currentTimeStep(), tmp);
    mpz_set_ui(_tickIncrement, mpf_get_ui(tmp));
    mpf_clear(tmp);
  }
}

void siconos::simulation::Event::setTick(double newTick)
{
  if (_eventCreated) {
    std::cout << "Warning: you change tick value for EventsManager -> a new initialization of "
                 "the object is required. "
              << std::endl;
  }
  _tick = newTick;
}

void siconos::simulation::Event::display() const
{
  std::cout << "===== Event data display =====" << std::endl;
  std::cout << " - Type: " << static_cast<std::underlying_type<EventType>::type>(_type)
            << "\n";
  std::cout << " - time (mpz_t format, double format): (";
  mpz_out_str(stdout, 10, _timeOfEvent);
  std::cout << ", " << _dTime << ")\n";
  std::cout << "===== End of Event display =====\n";
}

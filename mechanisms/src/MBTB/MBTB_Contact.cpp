/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
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

#include "MBTB_Contact.hpp"

#include <iostream>

#include "Interaction.hpp"
#include "MBTB_ContactRelation.hpp"
#include "MBTB_FC3DContactRelation.hpp"

siconos::mechanisms::MBTB_Contact::MBTB_Contact(int id,
                                                const std::string& ContactName,
                                                int indexBody1, int indexBody2,
                                                int indexCAD1, int indexCAD2,
                                                bool withFriction)
    : _ContactName{ContactName},
      _id{id},
      _indexBody1{indexBody1},
      _indexBody2{indexBody2},
      _indexCAD1{indexCAD1},
      _indexCAD2{indexCAD2},
      _withFriction{withFriction} {}

void siconos::mechanisms::MBTB_Contact::setInteraction(
    std::shared_ptr<siconos::modeling::Interaction> newInteraction) {
  std::cout << "siconos::mechanisms::MBTB_Contact::setInteraction(...)\n";
  std::cout << "_interaction before" << _interaction << "\n";
  std::cout << "newinteraction " << newInteraction << "\n";
  newInteraction->display();
  _interaction = newInteraction;
  std::cout << "_interaction after" << _interaction << "\n";
  _interaction->display();
}

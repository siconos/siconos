// SPDX-FileCopyrightText: 2026 INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

#include "MBTB_Contact.hpp"

#include <iostream>

#include "Interaction.hpp"
#include "MBTB_ContactRelation.hpp"
#include "MBTB_FC3DContactRelation.hpp"

siconos::mechanisms::MBTB_Contact::MBTB_Contact(int id, const std::string& ContactName,
                                                int indexBody1, int indexBody2, int indexCAD1,
                                                int indexCAD2, bool withFriction)
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

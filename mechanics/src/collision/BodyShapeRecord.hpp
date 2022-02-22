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

/*! \file BodyShapeRecord.hpp
  \brief Definition of an abstract Body shape record
  The objective of this class is to keep associate ds and static body woth the
shape in a contactor.
*/

#ifndef BodyShapeRecord_h
#define BodyShapeRecord_h

// We need to maintain a record associating each body with a shape,
// contactor, and collision object for each shape type.  We also need
// to access generic shape stuff (group, margin) by a pointer from the
// collision callback, so we need a record base class.
class BodyShapeRecord {
public:
  BodyShapeRecord(SP::SiconosVector b, SP::SecondOrderDS d, SP::SiconosShape sh,
                  SP::SiconosContactor con, SP::StaticBody staticCSR)
      : base(b), ds(d), sshape(sh), contactor(con),
        shape_version(sh->version()), staticBody(staticCSR)
  {
  }
  virtual ~BodyShapeRecord() {}

  SP::SiconosVector base;
  SP::SecondOrderDS ds;
  SP::SiconosShape sshape;
  SP::SiconosContactor contactor;
  unsigned int shape_version;
  SP::StaticBody staticBody;

  VIRTUAL_ACCEPT_VISITORS();
};

#endif /* BodyShapeRecord_h */

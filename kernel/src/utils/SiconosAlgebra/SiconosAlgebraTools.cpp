/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2020 INRIA.
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


#include "SiconosMatrix.hpp"
#include "BlockVector.hpp"
#include "SiconosVector.hpp"
#include "SiconosAlgebraTools.hpp"
#include <iostream>

namespace Siconos {
  namespace Algebra {

    bool isComparableTo(const BlockVector& v1, const BlockVector& v2)
    {
      // return:
      //  - true if both are block but with blocks which are facing each other of the same size.
      //  - false in other cases
      //
      const Index& I1 = *v1.tabIndex();
      const Index& I2 = *v2.tabIndex();

      return (I1 == I2);

    }

    bool isComparableTo(const  SiconosMatrix& m1, const  SiconosMatrix& m2)
    {
      // return:
      // - true if one of the matrices is a Simple and if they have the same dimensions.
      // - true if both are block but with blocks which are facing each other of the same size.
      // - false in other cases

      if((!m1.isBlock() || !m2.isBlock()) && (m1.size(0) == m2.size(0)) && (m1.size(1) == m2.size(1)))
        return true;

      const SP::Index I1R = m1.tabRow();
      const SP::Index I2R = m2.tabRow();
      const SP::Index I1C = m1.tabCol();
      const SP::Index I2C = m2.tabCol();

      return ((*I1R == *I2R) && (*I1C == *I2C));
    }
  }
}

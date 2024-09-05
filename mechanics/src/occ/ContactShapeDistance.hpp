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
#ifndef ContactShapeDistance_hpp
#define ContactShapeDistance_hpp

#include <gp_Pnt.hxx>
#include <gp_Dir.hxx>
namespace siconos::mechanics::occ {
  
struct ContactShapeDistance {
  double value{0.};

  // double x1{0.};
  // double y1{0.};
  // double z1{0.};

  // double x2{0.};
  // double y2{0.};
  // double z2{0.};

  // double nx{0.};
  // double ny{0.};
  // double nz{0.};

  gp_Pnt point1;
  gp_Pnt point2;
  gp_Dir normal;
  
  bool orientates(){
  if(gp_Vec{point1.Coord() - point2.Coord()}.Dot(normal) < 0.){
    normal.Reverse();
    return true;
  }
  return false;
  }
  
};
}  // namespace siconos::mechanics::occ

#endif

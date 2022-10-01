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

/*! \file Mesh.hpp

 */
#ifndef MATERIAL_H
#define MATERIAL_H

#include <vector>
#include <iostream>


/** A simple class for material
 */



enum ANALYSIS_TYPE_2D
{
  PLANE_STRAIN,
  PLANE_STRESS,
  AXYSYMMETRIC
};




class Material
{

protected:
  /* serialization hooks */
  //ACCEPT_SERIALIZATION(Mesh);

  /** mass density */
  double _massDensity;

  /** Young Modulus */
  double _elasticYoungModulus;

  /** Poison coefficient */
  double _poissonCoefficient;

  /** Analysis type in 2D */
  ANALYSIS_TYPE_2D _analysisType2D;

  /** thickness in 2D plane stress analysis */
  double _thickness;

  /** default constructor */
  Material() {};

public:

  /** constructor
   */
  Material(double massDensity, double ElasticYoungModulus,  double poissonCoefficient):
    _massDensity(massDensity),
    _elasticYoungModulus(ElasticYoungModulus),
    _poissonCoefficient(poissonCoefficient),
    _analysisType2D(PLANE_STRAIN),
    _thickness(1.0) {};


  /** destructor */
  ~Material() {};

  double massDensity()
  {
    return _massDensity;
  }

  double elasticYoungModulus()
  {
    return _elasticYoungModulus;
  }

  double poissonCoefficient()
  {
    return _poissonCoefficient;
  }

  ANALYSIS_TYPE_2D analysisType2D()
  {
    return _analysisType2D;
  }

  double thickness()
  {
    return _thickness;
  }

  /** print the data of the Material
   */
  void display(bool brief = true) const;

  //
  //ACCEPT_STD_VISITORS();

};

#endif // MESH_H

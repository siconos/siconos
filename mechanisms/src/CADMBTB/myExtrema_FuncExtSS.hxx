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

#ifndef _myExtrema_FuncExtSS_HeaderFile
#define _myExtrema_FuncExtSS_HeaderFile

#include <Adaptor3d_Surface.hxx>
#ifndef _gp_Pnt_HeaderFile
#include <gp_Pnt.hxx>
#endif
#ifndef _Standard_Real_HeaderFile
#include <Standard_Real.hxx>
#endif
#ifndef _TColStd_SequenceOfReal_HeaderFile
#include <TColStd_SequenceOfReal.hxx>
#endif
#ifndef _Extrema_SequenceOfPOnSurf_HeaderFile
#include <Extrema_SequenceOfPOnSurf.hxx>
#endif
#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
#ifndef _math_FunctionSetWithDerivatives_HeaderFile
#include <math_FunctionSetWithDerivatives.hxx>
#endif
#ifndef _Standard_Integer_HeaderFile
#include <Standard_Integer.hxx>
#endif
class Standard_OutOfRange;
class Adaptor3d_Surface;
class math_Vector;
class math_Matrix;
class Extrema_POnSurf;

#ifndef _Standard_HeaderFile
#include <Standard.hxx>
#endif
#ifndef _Standard_Macro_HeaderFile
#include <Standard_Macro.hxx>
#endif
/**
 * \brief This class has been built from OCC in view of overloding the distance
 * computation between CAD objects.
 */
//! Fonction permettant de rechercher les extrema de la <br>
//!          distance entre deux surfaces. <br>
class myExtrema_FuncExtSS : public math_FunctionSetWithDerivatives {
 public:
  void* operator new(size_t, void* anAddress) { return anAddress; }
  void* operator new(size_t size) { return Standard::Allocate(size); }
  void operator delete(void* anAddress) {
    if (anAddress) Standard::Free((Standard_Address&)anAddress);
  }
  // Methods PUBLIC
  //

  Standard_EXPORT myExtrema_FuncExtSS();

  Standard_EXPORT myExtrema_FuncExtSS(const Adaptor3d_Surface& S1,
                                      const Adaptor3d_Surface& S2);

  //! sets the field mysurf of the function. <br>
  Standard_EXPORT void Initialize(const Adaptor3d_Surface& S1, const Adaptor3d_Surface& S2);

  Standard_EXPORT Standard_Integer NbVariables() const;

  Standard_EXPORT Standard_Integer NbEquations() const;

  //! Calcul de Fi(U,V). <br>
  Standard_EXPORT Standard_Boolean Value(const math_Vector& UV, math_Vector& F);

  //! Calcul de Fi'(U,V). <br>
  Standard_EXPORT Standard_Boolean Derivatives(const math_Vector& UV, math_Matrix& DF);

  //! Calcul de Fi(U,V) et Fi'(U,V). <br>
  Standard_EXPORT Standard_Boolean Values(const math_Vector& UV, math_Vector& F,
                                          math_Matrix& DF);

  //! Memorise l'extremum trouve. <br>
  Standard_EXPORT virtual Standard_Integer GetStateNumber();

  //! Renvoie le nombre d'extrema trouves. <br>
  Standard_EXPORT Standard_Integer NbExt() const;

  //! Renvoie la valeur de la Nieme distance. <br>
  Standard_EXPORT Standard_Real Value(const Standard_Integer N) const;

  //! Renvoie le Nieme extremum sur S1. <br>
  Standard_EXPORT Extrema_POnSurf PointOnS1(const Standard_Integer N) const;

  //! Renvoie le Nieme extremum sur S2. <br>
  Standard_EXPORT Extrema_POnSurf PointOnS2(const Standard_Integer N) const;

 protected:
  // Methods PROTECTED
  //

  // Fields PROTECTED
  //

 private:
  // Methods PRIVATE
  //

  Standard_EXPORT Adaptor3d_Surface* Bidon() const;

  // Fields PRIVATE
  //
  Adaptor3d_Surface* myS1;
  Adaptor3d_Surface* myS2;
  gp_Pnt myP1;
  gp_Pnt myP2;
  Standard_Real myU1;
  Standard_Real myV1;
  Standard_Real myU2;
  Standard_Real myV2;
  TColStd_SequenceOfReal myValue;
  Extrema_SequenceOfPOnSurf myPoint1;
  Extrema_SequenceOfPOnSurf myPoint2;
  Standard_Boolean myS1init;
  Standard_Boolean myS2init;
};

// other Inline functions and methods (like "C++: function call" methods)
//

#endif

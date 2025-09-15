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

#ifndef _myExtrema_GenExtSS_HeaderFile
#define _myExtrema_GenExtSS_HeaderFile

#ifndef _Standard_Boolean_HeaderFile
#include <Standard_Boolean.hxx>
#endif
#ifndef _Standard_Real_HeaderFile
#include <Standard_Real.hxx>
#endif
#ifndef _Standard_Integer_HeaderFile
#include <Standard_Integer.hxx>
#endif
// #ifndef _Handle_TColgp_HArray2OfPnt_HeaderFile
// #include <Handle_TColgp_HArray2OfPnt.hxx>v
// #endif
#include <TColgp_HArray2OfPnt.hxx>
// #include "Standard_Handle.hxx"
// typedef Handle(TColgp_HArray2OfPnt) TColgp_HArray2OfPnt;

#ifndef _myExtrema_FuncExtSS_HeaderFile
#include "myExtrema_FuncExtSS.hxx"
#endif
#include <Adaptor3d_Surface.hxx>

class TColgp_HArray2OfPnt;
// class Handle_TColgp_HArray2OfPnt;
class StdFail_NotDone;
class Standard_OutOfRange;
class Standard_TypeMismatch;
class Adaptor3d_Surface;
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
//! It calculates all the extremum distances <br>
//!          between two surfaces. <br>
//!          These distances can be minimum or maximum. <br>
class myExtrema_GenExtSS {
 public:
  void* operator new(size_t, void* anAddress) { return anAddress; }
  void* operator new(size_t size) { return Standard::Allocate(size); }
  void operator delete(void* anAddress) {
    if (anAddress) Standard::Free((Standard_Address&)anAddress);
  }
  // Methods PUBLIC
  //

  Standard_EXPORT myExtrema_GenExtSS();

  //! It calculates all the distances. <br>
  //!          The function F(u,v)=distance(S1(u1,v1),S2(u2,v2)) has an <br>
  //!          extremum when gradient(F)=0. The algorithm searchs <br>
  //!          all the zeros inside the definition ranges of the <br>
  //!          surfaces. <br>
  //!          NbU and NbV are used to locate the close points <br>
  //!          to find the zeros. <br>
  Standard_EXPORT myExtrema_GenExtSS(const Adaptor3d_Surface& S1, const Adaptor3d_Surface& S2,
                                     const Standard_Integer NbU, const Standard_Integer NbV,
                                     const Standard_Real Tol1, const Standard_Real Tol2);

  //! It calculates all the distances. <br>
  //!          The function F(u,v)=distance(P,S(u,v)) has an <br>
  //!          extremum when gradient(F)=0. The algorithm searchs <br>
  //!          all the zeros inside the definition ranges of the <br>
  //!          surface. <br>
  //!          NbU and NbV are used to locate the close points <br>
  //!          to find the zeros. <br>
  Standard_EXPORT myExtrema_GenExtSS(const Adaptor3d_Surface& S1, const Adaptor3d_Surface& S2,
                                     const Standard_Integer NbU, const Standard_Integer NbV,
                                     const Standard_Real U1min, const Standard_Real U1sup,
                                     const Standard_Real V1min, const Standard_Real V1sup,
                                     const Standard_Real U2min, const Standard_Real U2sup,
                                     const Standard_Real V2min, const Standard_Real V2sup,
                                     const Standard_Real Tol1, const Standard_Real Tol2);

  Standard_EXPORT void Initialize(const Adaptor3d_Surface& S2, const Standard_Integer NbU,
                                  const Standard_Integer NbV, const Standard_Real Tol2);

  Standard_EXPORT void Initialize(const Adaptor3d_Surface& S2, const Standard_Integer NbU,
                                  const Standard_Integer NbV, const Standard_Real U2min,
                                  const Standard_Real U2sup, const Standard_Real V2min,
                                  const Standard_Real V2sup, const Standard_Real Tol2);

  //! the algorithm is done with S1 <br>
  //!          An exception is raised if the fields have not <br>
  //!          been initialized. <br>
  Standard_EXPORT void Perform(const Adaptor3d_Surface& S1, const Standard_Real Tol1);

  //! the algorithm is done withS1 <br>
  //!          An exception is raised if the fields have not <br>
  //!          been initialized. <br>
  Standard_EXPORT void Perform(const Adaptor3d_Surface& S1, const Standard_Real U1min,
                               const Standard_Real U1sup, const Standard_Real V1min,
                               const Standard_Real V1sup, const Standard_Real Tol1);

  //! Returns True if the distances are found. <br>
  Standard_EXPORT Standard_Boolean IsDone() const;

  //! Returns the number of extremum distances. <br>
  Standard_EXPORT Standard_Integer NbExt() const;

  //! Returns the value of the Nth resulting distance. <br>
  Standard_EXPORT Standard_Real Value(const Standard_Integer N) const;

  //! Returns the point of the Nth resulting distance. <br>
  Standard_EXPORT Extrema_POnSurf PointOnS1(const Standard_Integer N) const;

  //! Returns the point of the Nth resulting distance. <br>
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
  Standard_Boolean myDone;
  Standard_Boolean myInit;
  Standard_Real myu1min;
  Standard_Real myu1sup;
  Standard_Real myv1min;
  Standard_Real myv1sup;
  Standard_Real myu2min;
  Standard_Real myu2sup;
  Standard_Real myv2min;
  Standard_Real myv2sup;
  Standard_Integer myusample;
  Standard_Integer myvsample;
  Handle_TColgp_HArray2OfPnt mypoints1;
  Handle_TColgp_HArray2OfPnt mypoints2;
  Standard_Real mytol1;
  Standard_Real mytol2;
  myExtrema_FuncExtSS myF;
  // Adaptor3d_Surface* myS1;
  Adaptor3d_Surface* myS2;
  //  Standard_Real myGuessU1;
  //  Standard_Real myGuessV1;
  //  Standard_Real myGuessU2;
  //  Standard_Real myGuessV2;
};

// other Inline functions and methods (like "C++: function call" methods)
//

#endif

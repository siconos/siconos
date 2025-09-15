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

#ifndef _myBRepExtrema_ExtFF_HeaderFile
#define _myBRepExtrema_ExtFF_HeaderFile

#ifndef _myExtrema_ExtSS_HeaderFile
#include "myExtrema_ExtSS.hxx"
#endif
#ifndef _Standard_Integer_HeaderFile
#include <Standard_Integer.hxx>
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
#ifndef _Standard_Real_HeaderFile
#include <Standard_Real.hxx>
#endif

class BRepAdaptor_Surface;
class StdFail_NotDone;
class Standard_OutOfRange;
class Standard_TypeMismatch;
class TopoDS_Face;
class gp_Pnt;
// #include "Standard_Handle.hxx"

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

class myBRepExtrema_ExtFF {
 public:
  using Handle_BRepAdaptor_Surface = Handle(BRepAdaptor_Surface);

  void* operator new(size_t, void* anAddress) { return anAddress; }
  void* operator new(size_t size) { return Standard::Allocate(size); }
  void operator delete(void* anAddress) {
    if (anAddress) Standard::Free((Standard_Address&)anAddress);
  }
  // Methods PUBLIC
  //

  Standard_EXPORT myBRepExtrema_ExtFF();

  //! It calculates all the distances. <br>
  Standard_EXPORT myBRepExtrema_ExtFF(const TopoDS_Face& F1, const TopoDS_Face& F2);

  Standard_EXPORT myBRepExtrema_ExtFF(const TopoDS_Face& F1, const TopoDS_Face& F2, int id1,
                                      int id2);
  Standard_EXPORT void Initialize(const TopoDS_Face& F2);

  //! An exception is raised if the fields have not been <br>
  //!          initialized. <br>
  //!          Be careful: this method uses the Face F2 only for <br>
  //!          classify, not for the fields. <br>
  Standard_EXPORT void Perform(const TopoDS_Face& F1, const TopoDS_Face& F2);

  //! True if the distances are found. <br>
  Standard_EXPORT Standard_Boolean IsDone() const;

  //! Returns True if the surfaces are parallel. <br>
  Standard_EXPORT Standard_Boolean IsParallel() const;

  //! Returns the number of extremum distances. <br>
  Standard_EXPORT Standard_Integer NbExt() const;

  //! Returns the value of the <N>th extremum distance. <br>
  Standard_EXPORT Standard_Real Value(const Standard_Integer N) const;

  //! Returns the parameters on the  Face F1 of the  <N>th <br>
  //!          extremum distance. <br>
  Standard_EXPORT void ParameterOnFace1(const Standard_Integer N, Standard_Real& U,
                                        Standard_Real& V) const;

  //! Returns the parameters on the  Face F2 of the  <N>th <br>
  //!          extremum distance. <br>
  Standard_EXPORT void ParameterOnFace2(const Standard_Integer N, Standard_Real& U,
                                        Standard_Real& V) const;

  //! Returns the Point of the <N>th extremum distance. <br>
  Standard_EXPORT gp_Pnt PointOnFace1(const Standard_Integer N) const;

  //! Returns the Point of the <N>th extremum distance. <br>
  Standard_EXPORT gp_Pnt PointOnFace2(const Standard_Integer N) const;

  int _id1;
  int _id2;

 protected:
  // Methods PROTECTED
  //

  // Fields PROTECTED
  //

 private:
  // Methods PRIVATE
  //

  // Fields PRIVATE
  //
  myExtrema_ExtSS myExtrem;
  Standard_Integer mynbext;
  TColStd_SequenceOfReal mydist;
  Extrema_SequenceOfPOnSurf myPointsOnS1;
  Extrema_SequenceOfPOnSurf myPointsOnS2;
  Handle_BRepAdaptor_Surface myHS;
};

// other Inline functions and methods (like "C++: function call" methods)
//

#endif

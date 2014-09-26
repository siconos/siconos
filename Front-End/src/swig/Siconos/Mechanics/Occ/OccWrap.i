// -*- c++ -*-
// SWIG interface for Siconos Mechanics/Occ
%module(directors="1", allprotected="1") OccWrap


%include start.i

#undef WITH_IO
#undef WITH_SERIALIZATION

#ifdef WITH_IO
%{
#include <SiconosFull.hpp>
%}
#endif

%include picklable.i

%include path.i

%{
#include <SiconosKernel.hpp>
%}

%include handleException.i
%include sharedPointers.i
%include KernelTypes.i

%import Kernel.i

%{
#include <ExternalBody.hpp>
%}
%import ../ContactDetection/Base.i


%include pyRegister.i


%{
#include <MechanicsFwd.hpp>
%}
%include <MechanicsFwd.hpp>

// do not wrap visitor visit : this leads to a huge amount of wrapper
// code generation and this fails at compile time on shared_ptr freearg
%ignore visit;

%{
#include <TopoDS_Shape.hxx>
%}
%shared_ptr(TopoDS_Shape)
%import <Standard_Macro.hxx>
#define DEFINE_STANDARD_ALLOC
%include <TopoDS_Shape.hxx>

// force the definition of SWIGTYPE_p_Interaction...
typedef Interaction Interaction;

// due to undefined private copy constructors
%feature("notabstract") OccTimeStepping;

%typecheck(SWIG_TYPECHECK_INTEGER) (const OccContactShape & reference_shape) ()
%{
  // director mess
  int res;
  res = SWIG_ConvertPtr(argv[1], 0, SWIGTYPE_p_std11__shared_ptrT_OccContactShape_t, 0);
  _v = SWIG_CheckState(res);
  if(!_v)
  {
    res = SWIG_ConvertPtr(argv[1], 0, SWIGTYPE_p_OccContactShape, 0);
  }
%}


PY_FULL_REGISTER(ContactShapeDistance);
PY_FULL_REGISTER(OccContactShape);
PY_FULL_REGISTER(OccContactFace);
PY_FULL_REGISTER(OccContactEdge);
PY_FULL_REGISTER(ContactPoint);
PY_FULL_REGISTER(OccBody);
PY_FULL_REGISTER(OccR);
PY_FULL_REGISTER(OccTimeStepping);
PY_FULL_REGISTER(OccSpaceFilter);


%inline
%{
  #include <BRepTools.hxx>
  #include <BRep_Builder.hxx>
  #include <sstream>

  SP::OccContactShape importContactShape(PyObject * o)
  {

  if (PyObject_HasAttrString(o, "exportBrepToString"))
  {
    PyObject * ostr =
      PyObject_CallMethod(o,
                          const_cast<char *>("exportBrepToString"), NULL);

    std::stringstream brep_stream;
    brep_stream << PyString_AsString(ostr);

    SP::OccContactShape shape(new OccContactShape());
    BRep_Builder brep_builder;
    BRepTools::Read(shape->data(), brep_stream, brep_builder);

    shape->computeUVBounds();

    return shape;
  }
  else
  {
    /* void pointer should return None on Python side */
    return SP::OccContactShape();
  };

  }

  /* fix: use generated dynamic casting instead! */
  SP::OccBody cast_OccBody(SP::DynamicalSystem ds)
  {
    return std11::dynamic_pointer_cast<OccBody>(ds);
  };




%}


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

PY_FULL_REGISTER(OccShape);
PY_FULL_REGISTER(OccBody);
PY_FULL_REGISTER(OccR);

%inline
%{
  #include <BRepTools.hxx>
  #include <BRep_Builder.hxx>
  #include <sstream>

  SP::OccShape OccShapeFromFreeCAD(PyObject * o)
  {

  if (PyObject_HasAttrString(o, "exportBrepToString"))
  {
    PyObject * ostr =
      PyObject_CallMethod(o,
                          "exportBrepToString", NULL);

    std::stringstream brep_stream;
    brep_stream << PyString_AsString(ostr);

    SP::OccShape shape(new OccShape());
    BRep_Builder brep_builder;
    BRepTools::Read(*shape, brep_stream, brep_builder);

    return shape;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "cannot import brep");
  };

  }

%}

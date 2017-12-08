/* cadmbtb.i this file contains exported API of the CADMBTB library.*/
%module(package="siconos.mechanics.mechanisms") cadmbtb
%{
#include "CADMBTB_API.hpp"
#include "CADMBTB_PYTHON_API.hpp"
#include <string>
#include <sstream>
#include <BRepTools.hxx>
%}

// handle stl data types
%include stl.i

%include "CADMBTB_PYTHON_API.hpp"

%inline %{

std::string CADMBTB_TopoDSAsString(unsigned int numDS)
{
  std::stringstream out;
  const TopoDS_Shape& shape = CADMBTB_TopoDS(numDS);
  BRepTools::Write(shape, out);
  return out.str();
}

%}

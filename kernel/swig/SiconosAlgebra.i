// The SiconosAlgebra classes are handled seperately from the other Kernel classes
// This is because we define typemaps for them. It look like swig forgets some typemaps
// after being told more about a class (more specifically after applying PY_REGISTER_WITHOUT_DIRECTOR)
// Hence, we declare them fully here, and just after we define the typemaps (note the %include KernelTypes.i at the end)


%include KernelTypes.i

%warnfilter(509) SiconosMatrix::PLUForwardBackwardInPlace;
%warnfilter(509) SimpleMatrix::PLUForwardBackwardInPlace;
%warnfilter(509) SimpleMatrix::SolveByLeastSquares;
%warnfilter(504) SimpleMatrix::ACCEPT_STD_VISITORS();
%warnfilter(509) SiconosVector;

%ignore operator<<(std::ostream& os, const SiconosMatrix& bv);
%ignore operator<<(std::ostream& os, const SimpleMatrix& bv);
%ignore operator<<(std::ostream& os, const SiconosVector& bv);
%ignore operator<<(std::ostream& os, const BlockVector& bv);

%include SiconosMatrix.hpp
%include SimpleMatrix.hpp
%include SiconosVector.hpp
%include BlockVector.hpp

%extend SiconosMatrix{
  std::string __str__() { return $self->toString(); }
}
%extend SimpleMatrix{
  std::string __str__() { return $self->toString(); }
}
%extend SiconosVector{
  std::string __str__() { return $self->toString(); }
}
%extend BlockVector{
  std::string __str__() { return $self->toString(); }
}

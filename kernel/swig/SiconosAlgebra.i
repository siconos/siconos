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

/* Difficult to teach SIWG about the iterator base class, so disable
 * "nothing known about base class" warning. */
%warnfilter(401) SiconosVectorIterator;
%warnfilter(401) SiconosVectorConstIterator;

%ignore operator<<(std::ostream& os, const SiconosMatrix& bv);
%ignore operator<<(std::ostream& os, const SimpleMatrix& bv);
%ignore operator<<(std::ostream& os, const SiconosVector& bv);
%ignore operator<<(std::ostream& os, const BlockVector& bv);
%ignore operator std::vector<double>();

%ignore SiconosVectorIteratorTypeTpl::operator=;
%ignore SiconosVectorIteratorTypeTpl::operator++;
%ignore SiconosVectorIteratorTypeTpl;
%ignore SiconosVectorIteratorType;
%ignore SiconosVectorConstIteratorType;
%ignore SiconosVector::begin();
%ignore SiconosVector::end();
%ignore SiconosVector::begin() const;
%ignore SiconosVector::end() const;

%include SiconosMatrix.hpp
%include SimpleMatrix.hpp
%include SiconosVector.hpp
%include SiconosVectorIterator.hpp
%include BlockVector.hpp

%extend SiconosMatrix{
  std::string __str__() { return $self->toString(); }
}
%extend SimpleMatrix{
  std::string __str__() { return $self->toString(); }
}
%extend SiconosVector{
  std::string __str__() { return $self->toString(); }
  PyObject *__getitem__(size_t i) {
    if (i < $self->size())
      return PyFloat_FromDouble($self->getValue(i));
    else {
      PyErr_SetNone(PyExc_IndexError);
      return NULL;
    }
  }
  PyObject *__setitem__(size_t i, double value) {
    if (i < $self->size()) {
      $self->setValue(i, value);
      return Py_None;
    } else {
      PyErr_SetNone(PyExc_IndexError);
      return NULL;
    }
  }
  double __len__() { return $self->size(); }
  SiconosVectorIterator __iter__() {
    return SiconosVectorIterator($self->begin());
  }
%insert("python") %{
    def __array__(self):
        import numpy
        return numpy.fromiter(self, dtype=float)
%}
}
%extend BlockVector{
  std::string __str__() { return $self->toString(); }
}

%extend SiconosVectorIterator{
%insert("python") %{
    def __iter__(self):
        return SiconosVectorIterator(self)
%}
  PyObject* next()
  {
    if (*$self != (*$self).v->end())
      return PyFloat_FromDouble(*(*$self)++);
    PyErr_SetNone(PyExc_StopIteration);
    return NULL;
  }
%insert("python") %{
    __next__ = next
%}
}

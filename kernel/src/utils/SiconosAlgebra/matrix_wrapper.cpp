#include "matrix_wrapper.h"
#include "SiconosPointers.hpp"

MatrixWrapper::MatrixWrapper(siconos::algebra::SiconosMatrix &m) : _matrix(siconos::pointers::createSPtr(m)) {}

siconos::algebra::SiconosMatrix& MatrixWrapper::get_matrix() {
    return *_matrix;
}

std::shared_ptr<siconos::algebra::SiconosMatrix> MatrixWrapper::get_shared_ptr() {
    return _matrix;
}



#ifndef MATRIX_WRAPPER_H
#define MATRIX_WRAPPER_H

#include "SiconosMatrix.hpp"

class MatrixWrapper {
public:
    MatrixWrapper(siconos::algebra::SiconosMatrix &m);

    siconos::algebra::SiconosMatrix& get_matrix();

    std::shared_ptr<siconos::algebra::SiconosMatrix> get_shared_ptr();

private:
    std::shared_ptr<siconos::algebra::SiconosMatrix> _matrix;
};


#endif
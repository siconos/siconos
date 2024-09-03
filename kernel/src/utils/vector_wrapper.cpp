#include "vector_wrapper.h"
#include "SiconosPointers.hpp"

VectorWrapper::VectorWrapper(siconos::algebra::SiconosVector &v) : _vector(siconos::pointers::createSPtr(v)) {
    // std::cout << "contructor : " << _vector->contents.size() << "\t" << std::endl;
}

siconos::algebra::SiconosVector& VectorWrapper::get_vector() {
    std::cout << "get_vector : " << _vector->rows() << "\t" << _vector->cols() << std::endl;
    return *_vector;
}

std::shared_ptr<siconos::algebra::SiconosVector> VectorWrapper::get_shared_ptr() {
    std::cout << _vector->rows() << "\t" << _vector->cols() << std::endl;
    return _vector;
}

void VectorWrapper::print_vector() const {
    std::cout << "print_vector : " << _vector->rows() << "\t" << _vector->cols() << std::endl;
}

//************************************* */

MyVectorWrapper::MyVectorWrapper(MyVectorClass &v) : _vector(siconos::pointers::createSPtr(v)) {
    // std::cout << "contructor : " << _vector->contents.size() << "\t" << std::endl;
}


std::shared_ptr<MyVectorClass> MyVectorWrapper::get_shared_ptr() {
    std::cout << _vector->contents.size() << std::endl;
    return _vector;
}

void MyVectorWrapper::print_vector() const {
    std::cout << "MyVectorWrapper_print_vector " << _vector->contents[0] << "\t "  << _vector->contents[1] << "\t "  << _vector->contents[2] << "\t "<< std::endl;
}



//************************************* */


MyEigenWrapper::MyEigenWrapper(MyEigenClass &v) : _vector(siconos::pointers::createSPtr(v)) {
    // std::cout << "contructor : " << _vector->contents.size() << "\t" << std::endl;
}


std::shared_ptr<MyEigenClass> MyEigenWrapper::get_shared_ptr() {
    std::cout << _vector->contents.size() << std::endl;
    return _vector;
}

void MyEigenWrapper::change_vector() {
    _vector->contents(0) = 11.;
}

void MyEigenWrapper::print_vector() const {
    std::cout << "MyEigenWrapper_print_vector " << _vector->contents[0] << "\t "  << _vector->contents[1] << "\t "  << _vector->contents[2] << "\t "<< std::endl;
}

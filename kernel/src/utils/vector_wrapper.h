#ifndef VECTOR_WRAPPER_H
#define VECTOR_WRAPPER_H

#include "SiconosVector.hpp"

struct MyVectorClass {
    MyVectorClass(std::vector<float> &v) : contents(v){};
    std::vector<float> contents;
    void print_v() {
        std::cout << contents[0] << std::endl;
    }
};

struct MyEigenClass {
    siconos::algebra::SiconosVector contents;
};


class MyVectorWrapper {
public:
    MyVectorWrapper(MyVectorClass& v);

    std::shared_ptr<MyVectorClass> get_shared_ptr();
    
    void print_vector() const;

private:
    std::shared_ptr<MyVectorClass> _vector {nullptr};
};


class MyEigenWrapper {
public:
    MyEigenWrapper(MyEigenClass& v);

    std::shared_ptr<MyEigenClass> get_shared_ptr();

    void change_vector();
    
    void print_vector() const;

private:
    std::shared_ptr<MyEigenClass> _vector {nullptr};
};






class VectorWrapper {
public:
    VectorWrapper(siconos::algebra::SiconosVector& v);

    siconos::algebra::SiconosVector& get_vector();

    std::shared_ptr<siconos::algebra::SiconosVector> get_shared_ptr();
    
    void print_vector() const;

private:
    std::shared_ptr<siconos::algebra::SiconosVector> _vector {nullptr};
};

#endif
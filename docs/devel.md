  
# Rules for developers

## 1 - Default behavior
If you don't know what to do, follow : [C++ core guidelines](https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#S-introduction)

Style: respect [Google C++ style guide](https://google.github.io/styleguide/cppguide.html)

## Plugin function types (for DS, relations ...)

- Declare all function prototypes in FunctionTypes.hpp
- Naming convention: FunctionXYZ..._O or RFunctionXYZ_O with: 
    -  XYZ... the type of read-only in parameters, with S for scalar type, V for vector type, M for dense matrix type and Ms for sparse matrix type
    - O(=S, V, M or Ms) the read-write in-out parameter used to collect the result
    - RFunction if the return type is non void (and corresponds to the read-write in-out parameter)

    example:
    ```c++
    using FunctionVS_M =
    std::function<void(const Eigen::Ref<const siconos::algebra::SiconosVector>&,
                       double, Eigen::Ref<siconos::algebra::MapType>)>;

    ```
  
 ## access (method) to stl container type class attributes

- Assuming an attribute of stl type, use std::span to provide access to its content.

Example
```
private:
    std::vector<T> attr_;

public:

    std::span<const T> attr() const {return attr_};
```

## pass eigen objects to functions or methods


- if the function is to be wrapped with pybind11, use Eigen::Ref, if not use const&

Why? See [Eigen doc - Pass by reference](https://github.com/pybind/pybind11/blob/master/docs/advanced/cast/eigen.rst)

Example:

```
somefunc(const Eigen::Ref<const siconos::algebra::SiconosVector>& ); // pb11
// or
somefunc(const siconos::algebra;;SiconosVector&); // not for pb11
```
This is required for a proper/optimized wrapping between numpy like objects and Eigen vectors

Many examples are available in Lagrangian DS classes.

## Read-only access to a std::shared<SiconosVector> or std::unique<SiconosVector> attribute

Rq: if it's a read-only attribute, unique_ptr<T> (compared to std::shared) is probably the best choice for the arg. type
And it's worth checking if "T" only won't be better.

```
std::unique_ptr<siconos::algebra::SiconosVector> name_
siconos::algebra::SiconosVector name2_

inline auto name() const {
    return siconos::algebra::ConstMapVectorType(name_->data(), name_->size());
  }

inline auto name2() const {
    return siconos::algebra::ConstMapVectorType(name2_.data(), name2_.size());
  }

```



## Work vectors in integrators

#### related to DS

Each osi has usually access to work_ds: a container of vectors, saved in the graph of ds.

```c++
// For a given graph node (ds) dsi:
work_ds = *_dynamicalSystemsGraph->properties(*dsi).workVectors;
```

- Access to its elements (SiconosVector) is defined thanks to an enum. 
- Each integrator has its own enum

Example:
```c++
auto vfree = work_ds[tools::enum_to_index(work::ds::vfree)];
```
- enum must be a class enum
- if a derived class has its own enum, any method that uses this work_ds must be re-implemented in derived classes
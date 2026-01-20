  
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
  

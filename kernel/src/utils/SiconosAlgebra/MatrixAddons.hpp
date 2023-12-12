// LU
// template<typename OtherDerived>
// Eigen::FullPivLU<MatrixBase<Derived>> lu = Eigen::FullPivLU<MatrixBase<Derived>>::FullPivLU();
bool isFactorized = false;
FullPivLU<Matrix>* lu_siconos = nullptr;
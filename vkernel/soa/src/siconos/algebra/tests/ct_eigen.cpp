#include "siconos/algebra/eigen.hpp"

namespace siconos::algebra {
static_assert(vector<int, 3>::SizeAtCompileTime != Eigen::Dynamic);
static_assert(unbounded_vector<int>::SizeAtCompileTime == Eigen::Dynamic);
static_assert(vector<int, 3>::ColsAtCompileTime == 1);
static_assert(match::fixed_size_matrix<matrix<int, 1, 1>>);
static_assert(
    !match::fixed_size_matrix<Eigen::Matrix<int, Eigen::Dynamic, 1>>);
static_assert(!match::fixed_size_matrix<unbounded_matrix<int>>);
static_assert(match::fixed_size_vector<vector<int, 3>>);
static_assert(!match::fixed_size_vector<unbounded_vector<int>>);
static_assert(match::variable_size_matrix<unbounded_matrix<int>>);
static_assert(match::variable_size_vector<unbounded_vector<int>>);
static_assert(match::any_matrix<matrix<int, 2, 2>>);
static_assert(match::any_matrix<unbounded_matrix<int>>);
static_assert(match::any_matrix<Eigen::DiagonalMatrix<int, 3>>);
static_assert(!match::matrix<int>);
static_assert(!match::vector<int>);

static_assert(std::is_same_v<matrix<int, 1, 2>, trans_t<matrix<int, 2, 1>>>);
static_assert(std::is_same_v<prod_t<matrix<int, 2, 2>, matrix<int, 2, 1>>,
                             matrix<int, 2, 1>>);
static_assert(std::is_same_v<prod_t<matrix<int, 1, 2>, matrix<int, 2, 1>>,
                             matrix<int, 1, 1>>);
static_assert(
    std::is_same_v<prod_t<matrix<int, 1, 2>, trans_t<matrix<int, 1, 2>>>,
                   matrix<int, 1, 1>>);

}  // namespace siconos::algebra
int main() {};

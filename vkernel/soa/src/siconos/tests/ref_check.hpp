#pragma once
// Reference data for the vkernel tests: read a whitespace separated matrix and
// compare it to a computed one. Eigen only: the vkernel tests link numerics,
// not kernel, so algebra::io::compareRefFile is not available here.

#include <cppunit/TestAssert.h>

#include <Eigen/Dense>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace siconos::test {

namespace detail {

inline bool is_unsigned_integer(const std::string& word) {
  if (word.empty()) return false;
  for (char c : word)
    if (c < '0' || c > '9') return false;
  return true;
}

inline double to_double(const std::string& path, const std::string& word) {
  std::size_t used = 0;
  double value = 0.;
  try {
    value = std::stod(word, &used);
  } catch (const std::exception&) {
    used = 0;
  }
  if (used != word.size()) throw std::runtime_error(path + ": not a number: '" + word + "'");
  return value;
}

// the words of one line, empty lines giving none
inline std::vector<std::string> words_of(const std::string& line) {
  std::istringstream tokens(line);
  std::vector<std::string> words;
  for (std::string word; tokens >> word;) words.push_back(word);
  return words;
}

}  // namespace detail

/// Reads every whitespace separated number of `path` into a matrix. A leading
/// line of exactly two integer tokens ("rows cols", the tutorials .ref format)
/// is skipped; the column count comes from the first data line, which must be
/// 5 columns, 7 columns, ... consistently for the whole file.
inline Eigen::MatrixXd read_ref(const std::string& path) {
  std::ifstream file(path);
  if (!file) throw std::runtime_error(path + ": cannot open reference file");

  std::vector<double> values;
  std::size_t ncols = 0;
  bool first_line = true;
  std::string line;
  while (std::getline(file, line)) {
    const auto words = detail::words_of(line);
    if (words.empty()) continue;
    if (first_line) {
      first_line = false;
      // optional "rows cols" header of the tutorials .ref files
      if (words.size() == 2 && detail::is_unsigned_integer(words[0]) &&
          detail::is_unsigned_integer(words[1]))
        continue;
    }
    if (ncols == 0)
      ncols = words.size();
    else if (words.size() != ncols)
      throw std::runtime_error(path + ": line with " + std::to_string(words.size()) +
                               " values, expected " + std::to_string(ncols));
    for (const auto& word : words) values.push_back(detail::to_double(path, word));
  }
  if (ncols == 0) throw std::runtime_error(path + ": no value in the reference file");

  Eigen::MatrixXd matrix(values.size() / ncols, ncols);
  for (std::size_t i = 0; i < values.size(); ++i) matrix(i / ncols, i % ncols) = values[i];
  return matrix;
}

/// What a comparison found: the relative error, and the worst element with its
/// (0 based) row and column and its absolute difference.
struct RefError {
  double relative;
  double worst;
  Eigen::Index row;
  Eigen::Index col;
};

/// ||dat - ref||_F / max(||ref||_F, 1), the largest element difference and where
/// it is. Both matrices must have the same dimensions.
inline RefError relative_error(const Eigen::MatrixXd& dat, const Eigen::MatrixXd& ref) {
  if (dat.rows() != ref.rows() || dat.cols() != ref.cols())
    throw std::runtime_error("relative_error: dimension mismatch");

  const Eigen::MatrixXd difference = dat - ref;
  RefError error;
  error.worst = difference.cwiseAbs().maxCoeff(&error.row, &error.col);
  error.relative = difference.norm() / std::max(ref.norm(), 1.);
  return error;
}

}  // namespace siconos::test

/// Asserts that the matrix `dat` matches the reference file `ref_file` within a
/// relative Frobenius tolerance `tol`. A dimension mismatch is an explicit
/// failure, never a partial comparison; a failure reports both dimensions, the
/// computed error, the tolerance and the worst element with its 1 based
/// position.
#define CPPUNIT_ASSERT_REF(dat, ref_file, tol)                                               \
  do {                                                                                       \
    const Eigen::MatrixXd& ref_check_dat = (dat);                                            \
    const Eigen::MatrixXd ref_check_ref = ::siconos::test::read_ref(ref_file);               \
    {                                                                                        \
      std::ostringstream ref_check_msg;                                                      \
      ref_check_msg << ref_file << ": computed " << ref_check_dat.rows() << "x"              \
                    << ref_check_dat.cols() << ", reference " << ref_check_ref.rows() << "x" \
                    << ref_check_ref.cols();                                                 \
      CPPUNIT_ASSERT_MESSAGE(ref_check_msg.str(),                                            \
                             ref_check_dat.rows() == ref_check_ref.rows() &&                 \
                                 ref_check_dat.cols() == ref_check_ref.cols());              \
    }                                                                                        \
    const auto ref_check_error =                                                             \
        ::siconos::test::relative_error(ref_check_dat, ref_check_ref);                       \
    if (!(ref_check_error.relative <= (tol))) { /* NaN fails too */                          \
      std::ostringstream ref_check_msg;                                                      \
      ref_check_msg << ref_file << ": relative Frobenius error " << ref_check_error.relative \
                    << " > tolerance " << (tol) << "; worst element at row "                 \
                    << ref_check_error.row + 1 << ", column " << ref_check_error.col + 1     \
                    << ": " << ref_check_error.worst;                                        \
      CPPUNIT_FAIL(ref_check_msg.str());                                                     \
    }                                                                                        \
  } while (false)

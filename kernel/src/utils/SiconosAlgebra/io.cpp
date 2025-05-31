/* Siconos is a program dedicated to modeling, simulation and control
 * of non smooth dynamical systems.
 *
 * Copyright 2024 INRIA.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
// #include "SiconosAlgebra.hpp"
#include "io.hpp"

#include <filesystem>
#include <fstream>
#include <iterator>

#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"
// #include <boost/numeric/ublas/io.hpp>

bool siconos::algebra::io::read(const std::string &fileName, SiconosVector &m,
                                const std::ios_base::openmode &mode, int prec,
                                WriteType inputType, const std::ios::fmtflags &flags) {
  // Read and check the file
  std::ifstream infile(fileName, mode);
  std::filesystem::path p1{fileName};
  if (std::filesystem::is_empty(p1)) {
    THROW_EXCEPTION("The given file is empty!");
  }

  if (!infile.good()) THROW_EXCEPTION("");

  infile.flags(flags);
  infile.precision(prec);

  if (mode == BINARY_IN) {
    double *x = m.data();
    if (inputType == WriteType::python) {
      unsigned int dim;
      infile.read((char *)(&dim), sizeof(m.size()));
    }
    infile.read((char *)(&x[0]), m.size() * sizeof(double));
  } else {
    // Read the dimension of the vector in the first line of the input file
    // Just use to check that sizes are consistents.
    if (inputType == WriteType::python) {
      unsigned int dim;
      infile >> dim;
      if (dim != m.size()) m.resize(dim);
    }
    copy((std::istream_iterator<double>(infile)), std::istream_iterator<double>(),
         (m.begin()));
  }
  infile.close();
  return true;
}

bool siconos::algebra::io::write(const std::string &fileName, const SiconosVector &m,
                                 const std::ios_base::openmode &mode, int prec,
                                 const WriteType outputType, const std::ios::fmtflags &flags) {
  std::ofstream outfile(fileName, mode);
  outfile.flags(flags);

  if (!outfile.good()) THROW_EXCEPTION("");
  outfile.precision(prec);
  if (mode == BINARY_OUT) {
    const double *x = m.data();
    if (outputType == WriteType::python) {
      unsigned int dim = m.size();
      outfile.write((char *)&dim, sizeof(dim));
    }
    outfile.write((char *)(&x[0]), sizeof(double) * m.size());
  } else {
    if (outputType == WriteType::python) outfile << m.size() << std::endl;
    std::copy(m.begin(), m.end(), std::ostream_iterator<double>(outfile, " "));
  }
  outfile.close();
  return true;
}

bool siconos::algebra::io::read(const std::string &filename, SiconosMatrix &m,
                                const std::ios_base::openmode &mode) {
  std::ifstream infile(filename, mode);

  if (!infile.good()) THROW_EXCEPTION("");

  if (infile.peek() == std::ifstream::traits_type::eof()) {
    THROW_EXCEPTION("the given file is empty!");
  }

  infile.precision(15);
  infile.setf(std::ios::scientific);

  // Dim of the matrix are given in the first line.
  // Just use to check that sizes are consistents.

  unsigned int s1, s2;
  infile >> s1;
  infile >> s2;

  if (s1 != m.rows() || s2 != m.cols()) m.resize(s1, s2);

  // Note: using istream stl iterator seems to be 2-times faster than << with a
  // loop over matrix data.
  //  copy((std::istream_iterator<double>(infile)),
  //  std::istream_iterator<double>(), (p->data()).begin());
  // But it fails with column-major saving ... (ok if user write its matrix in a
  // column-major way)

  for (unsigned int i = 0; i < s1; i++) {
    for (unsigned int j = 0; j < s2; j++) {
      infile >> m(i, j);
    }
  }

  infile.close();
  return true;
}

bool siconos::algebra::io::write(const std::string &filename, const SiconosMatrix &m,
                                 const std::ios_base::openmode &mode, WriteType outputType) {
  // Open file and various checks
  std::ofstream outfile(filename, mode);

  if (!outfile.good()) THROW_EXCEPTION("");

  // if (m.isBlock()) THROW_EXCEPTION("not yet implemented for BlockMatrix");

  outfile.precision(15);
  outfile.setf(std::ios::scientific);
  // Writing

  if (outputType == WriteType::python) outfile << m.rows() << " " << m.cols() << std::endl;

  double tmp;
  for (decltype(m.rows()) i = 0; i < m.rows(); i++) {
    for (decltype(m.cols()) j = 0; j < m.cols(); j++) {
      tmp = m(i, j);
      if (fabs(tmp) < std::numeric_limits<double>::min()) tmp = 0.0;
      outfile << tmp << " ";
      assert(outfile.good());
    }
    outfile << std::endl;
  }

  outfile.close();
  return true;
}

double siconos::algebra::io::compareRefFile(const SiconosMatrix &data, std::string filename,
                                            double epsilon, std::vector<int> index,
                                            const std::ios_base::openmode mode, bool verbose) {
  auto ref = std::make_shared<SiconosMatrix>(0, 0);
  bool compare = false;
  // SiconosMatrix ref{0, 0};
  try {
    compare = read(filename, *ref, mode);
  } catch (...) {
    if (verbose)
      std::cout << "Warning: reference file " << filename
                << " not found, no comparison performed." << std::endl;
    siconos::exception::process();
  }
  if (!compare) return -1.0;

  if (verbose) std::cout << "\n ===> Comparison with reference file " << filename << std::endl;

  SiconosVector err(data.cols());
  siconos::algebra::normInfByColumn(data - *ref, err);

  if (verbose) {
    std::cout << "Error vector:\n";
    siconos::algebra::print(err);
  }

  double error = 0.;
  /* Scalar error = max of columns */
  if (index.empty()) {
    error = *std::max_element(err.begin(), err.end());
  } else {
    for (auto &i : index) {
      if (error < err(i)) error = err(i);
    }
  }
  error = std::abs(error);

  if (verbose) std::cout << "\nError max = " << error << "\n\n";
  if (error > epsilon) {
    if (verbose) {
      std::cout << "\nWarning. The results are rather different from the "
                   "reference file.\n\n";
    }
  }

  return error;
}

siconos::algebra::SiconosVector siconos::algebra::io::readVectorFromJson(
    const nlohmann::json &jin) {
  // jin might be a vector of double or a list of points
  // so we must update vec_size by checking the first element in jin
  auto first = jin[0];
  auto vec_size = jin.size();
  if (first.is_array()) vec_size *= first.size();
  siconos::algebra::SiconosVector vec{vec_size};
  size_t element_index = 0;
  for (const auto &element : jin) {
    if (element.is_array()) {
      for (const auto &sub_element : element)
        vec.coeffRef(element_index++) = sub_element.get<double>();
    } else
      vec.coeffRef(element_index++) = element.get<double>();
  }

  return vec;  // RVO
}

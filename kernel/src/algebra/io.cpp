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

#include "SiconosAlgebraAddons.hpp"
#include "SiconosException.hpp"
#include "SiconosMatrix.hpp"
#include "SiconosVector.hpp"

bool siconos::algebra::io::read(const std::string& fileName, SiconosVector& m,
                                const std::ios_base::openmode& mode, int prec,
                                WriteType inputType, const std::ios::fmtflags& flags) {
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
    double* x = m.data();
    if (inputType == WriteType::python) {
      siconos::algebra::Index dim;
      infile.read((char*)(&dim), sizeof(m.size()));
    }
    infile.read((char*)(&x[0]), m.size() * sizeof(double));
  } else {
    // Read the dimension of the vector in the first line of the input file
    // Just use to check that sizes are consistents.
    if (inputType == WriteType::python) {
      siconos::algebra::Index dim;
      infile >> dim;
      if (dim != m.size()) m.resize(dim);
    }
    copy((std::istream_iterator<double>(infile)), std::istream_iterator<double>(),
         (m.begin()));
  }
  infile.close();
  return true;
}

bool siconos::algebra::io::write(const std::string& fileName, const SiconosVector& m,
                                 const std::ios_base::openmode& mode, int prec,
                                 const WriteType outputType, const std::ios::fmtflags& flags) {
  std::ofstream outfile(fileName, mode);
  outfile.flags(flags);

  if (!outfile.good()) THROW_EXCEPTION("");
  outfile.precision(prec);
  if (mode == BINARY_OUT) {
    const double* x = m.data();
    if (outputType == WriteType::python) {
      siconos::algebra::Index dim = m.size();
      outfile.write((char*)&dim, sizeof(dim));
    }
    outfile.write((char*)(&x[0]), sizeof(double) * m.size());
  } else {
    if (outputType == WriteType::python) outfile << m.size() << std::endl;
    std::copy(m.begin(), m.end(), std::ostream_iterator<double>(outfile, " "));
  }
  outfile.close();
  return true;
}

siconos::algebra::SiconosDenseMatrix siconos::algebra::io::readDenseMatrix(
    const std::string& filename, const std::ios_base::openmode& mode) {
  std::ifstream infile(filename, mode);
  if (!infile.is_open()) throw std::runtime_error("Cannot open file: " + filename);

  if (infile.peek() == std::ifstream::traits_type::eof())
    throw std::runtime_error("The given file is empty!");

  // Binary read
  if (mode & std::ios::binary) {
    char magic[4];
    uint16_t version, type;
    siconos::algebra::Index rows, cols;

    infile.read(magic, 4);
    infile.read(reinterpret_cast<char*>(&version), sizeof(version));
    infile.read(reinterpret_cast<char*>(&type), sizeof(type));
    infile.read(reinterpret_cast<char*>(&rows), sizeof(rows));
    infile.read(reinterpret_cast<char*>(&cols), sizeof(cols));

    if (std::strncmp(magic, "SICO", 4) != 0 || type != 1)
      throw std::runtime_error("Invalid file format or not a dense matrix");

    SiconosDenseMatrix m(rows, cols);
    infile.read(reinterpret_cast<char*>(m.data()), sizeof(double) * rows * cols);

    infile.read(reinterpret_cast<char*>(&rows), sizeof(int));
    infile.read(reinterpret_cast<char*>(&cols), sizeof(int));

    SiconosDenseMatrix mat(rows, cols);
    infile.read(reinterpret_cast<char*>(mat.data()), sizeof(double) * rows * cols);
    return mat;
  } else {
    siconos::algebra::Index rows, cols;
    infile >> rows >> cols;

    SiconosDenseMatrix m(rows, cols);

    for (siconos::algebra::Index i = 0; i < rows; ++i)
      for (siconos::algebra::Index j = 0; j < cols; ++j) infile >> m(i, j);

    return m;  // RVO or move
  }
}

siconos::algebra::SiconosSparseMatrix siconos::algebra::io::readSparseMatrix(
    const std::string& filename, const std::ios_base::openmode& mode) {
  std::ifstream infile(filename, mode);
  if (!infile.is_open()) {
    throw std::runtime_error("Cannot open file: " + filename);
  }

  if (infile.peek() == std::ifstream::traits_type::eof()) {
    throw std::runtime_error("The given file is empty!");
  }

  if (mode & std::ios_base::binary) {
    char magic[4];
    infile.read(magic, 4);
    if (std::string(magic, 4) != "SICO") throw std::runtime_error("Invalid magic header");

    int rows, cols, nnz;
    infile.read(reinterpret_cast<char*>(&rows), sizeof(int));
    infile.read(reinterpret_cast<char*>(&cols), sizeof(int));
    infile.read(reinterpret_cast<char*>(&nnz), sizeof(int));

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz);
    for (int k = 0; k < nnz; ++k) {
      int i, j;
      double v;
      infile.read(reinterpret_cast<char*>(&i), sizeof(int));
      infile.read(reinterpret_cast<char*>(&j), sizeof(int));
      infile.read(reinterpret_cast<char*>(&v), sizeof(double));
      triplets.emplace_back(i, j, v);
    }
    SiconosSparseMatrix mat(rows, cols);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
  } else {
    int rows, cols, nnz;
    infile >> rows >> cols >> nnz;
    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(nnz);
    for (int k = 0; k < nnz; ++k) {
      int i, j;
      double v;
      infile >> i >> j >> v;
      triplets.emplace_back(i, j, v);
    }
    SiconosSparseMatrix mat(rows, cols);
    mat.setFromTriplets(triplets.begin(), triplets.end());
    return mat;
  }
}

void siconos::algebra::io::write(const std::string& filename,
                                 const siconos::algebra::SiconosDenseMatrix& mat,
                                 const std::ios_base::openmode& mode,
                                 const WriteType outputType) {
  std::ofstream outfile(filename, mode);
  if (!outfile.is_open()) throw std::runtime_error("Cannot open file " + filename);

  if (mode & std::ios::binary) {
    // We enforce row/col output when binary is on.
    assert(outputType == WriteType::python);
    // Unique id for siconos binary files
    const char magic[4] = {'S', 'I', 'C', 'O'};
    uint16_t version = 1;
    uint16_t type = 1;
    // int32_t rows = static_cast<int32_t>(mat.rows());
    // int32_t cols = static_cast<int32_t>(mat.cols());
    auto rows = mat.rows();
    auto cols = mat.cols();

    outfile.write(magic, 4);
    outfile.write(reinterpret_cast<const char*>(&version), sizeof(version));
    outfile.write(reinterpret_cast<const char*>(&type), sizeof(type));
    outfile.write(reinterpret_cast<const char*>(&rows), sizeof(rows));
    outfile.write(reinterpret_cast<const char*>(&cols), sizeof(cols));
    outfile.write(reinterpret_cast<const char*>(mat.data()), sizeof(double) * rows * cols);
  } else {
    outfile.precision(15);
    outfile.setf(std::ios::scientific);
    if (outputType == WriteType::python) outfile << mat.rows() << " " << mat.cols() << "\n";
    for (siconos::algebra::Index i = 0; i < mat.rows(); ++i) {
      for (siconos::algebra::Index j = 0; j < mat.cols(); ++j) {
        auto tmp = mat(i, j);
        if (std::abs(tmp) < std::numeric_limits<double>::min()) tmp = 0.0;
        outfile << tmp << " ";
      }
      outfile << "\n";
    }
  }
}

// --- Sparse Write ---
void siconos::algebra::io::writeSparseMatrix(const std::string& filename,
                                             const SiconosSparseMatrix& mat,
                                             const std::ios_base::openmode& mode) {
  std::ofstream outfile(filename, mode);
  if (!outfile.is_open()) throw std::runtime_error("Cannot open file " + filename);

  if (mode & std::ios_base::binary) {
    // Check for unique id in siconos binary files
    outfile.write("SICO", 4);
    siconos::algebra::Index rows = mat.rows(), cols = mat.cols(), nnz = mat.nonZeros();
    outfile.write(reinterpret_cast<const char*>(&rows), sizeof(siconos::algebra::Index));
    outfile.write(reinterpret_cast<const char*>(&cols), sizeof(siconos::algebra::Index));
    outfile.write(reinterpret_cast<const char*>(&nnz), sizeof(siconos::algebra::Index));

    for (int k = 0; k < mat.outerSize(); ++k)
      for (SiconosSparseMatrix::InnerIterator it(mat, k); it; ++it) {
        siconos::algebra::Index i = it.row(), j = it.col();
        double v = it.value();
        outfile.write(reinterpret_cast<const char*>(&i), sizeof(siconos::algebra::Index));
        outfile.write(reinterpret_cast<const char*>(&j), sizeof(siconos::algebra::Index));
        outfile.write(reinterpret_cast<const char*>(&v), sizeof(double));
      }
  } else {
    outfile << mat.rows() << " " << mat.cols() << " " << mat.nonZeros() << "\n";
    for (siconos::algebra::Index k = 0; k < mat.outerSize(); ++k)
      for (SiconosSparseMatrix::InnerIterator it(mat, k); it; ++it)
        outfile << it.row() << " " << it.col() << " " << it.value() << "\n";
  }
}

double siconos::algebra::io::compareRefFile(const SiconosDenseMatrix& data,
                                            std::string filename, double epsilon,
                                            std::vector<int> index,
                                            const std::ios_base::openmode mode, bool verbose) {
  // SiconosMatrix ref{0, 0};
  auto ref = readDenseMatrix(filename, mode);

  if (verbose) std::cout << "\n ===> Comparison with reference file " << filename << std::endl;
  auto err = siconos::algebra::normInfByColumn(data - ref);

  if (verbose) {
    std::cout << "Error vector:\n";
    siconos::algebra::print(err);
  }

  double error = 0.;
  /* Scalar error = max of columns */
  if (index.empty()) {
    error = *std::max_element(err.begin(), err.end());
  } else {
    for (auto& i : index) {
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
    const nlohmann::json& jin) {
  // jin might be a vector of double or a list of points
  // so we must update vec_size by checking the first element in jin
  auto first = jin[0];
  auto vec_size = jin.size();
  if (first.is_array()) vec_size *= first.size();
  siconos::algebra::SiconosVector vec{vec_size};
  Index element_index = 0;
  for (const auto& element : jin) {
    if (element.is_array()) {
      for (const auto& sub_element : element)
        vec.coeffRef(element_index++) = sub_element.get<double>();
    } else
      vec.coeffRef(element_index++) = element.get<double>();
  }

  return vec;  // RVO
}

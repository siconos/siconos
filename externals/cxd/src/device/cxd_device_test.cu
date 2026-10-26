/* CXD is a program dedicated to convex collision detection.
 * It is part of the Siconos Software.
 *
 * Copyright 2026 INRIA.
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
// The device backend's parity test: the deterministic corpus, the committed
// golden files, the two device passes, and the phase-4 grading over both.
//
// The rig mirrors grade_batch_against_golden's sweep exactly: directions
// outermost, then the forward/mirrored shape order, then the cases; the
// pool is the mirrored a/b/b/a packing; statuses grade exactly; hits grade
// {depth, normal} within the phase-4 bands (double 1e-6 / 1e-4 rads; float
// 1e-4 / 1e-3); refused rows are expected NaN in place from the device.
// Alignment against the committed rows is name-checked per block, so a
// corpus drift cannot silently pass over the wrong reference line.
//
// Exit code 0 = parity held. Any mismatch prints the row and fails.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "batch_parity_corpus.hpp"
#include "cxd_batch_device.cuh"

#ifndef CXD_TEST_GOLDEN_DIR
#define CXD_TEST_GOLDEN_DIR "."
#endif

namespace {

// The angle between two hit normals, in radians; both widened to double.
template <std::size_t Dim>
double normal_angle(const std::array<double, Dim>& a,
                    const std::array<double, Dim>& b)
{
  double dot = 0.;
  double length_a = 0.;
  double length_b = 0.;
  for (std::size_t k = 0; k < Dim; ++k) {
    dot += a[k] * b[k];
    length_a += a[k] * a[k];
    length_b += b[k] * b[k];
  }
  length_a = std::sqrt(length_a);
  length_b = std::sqrt(length_b);
  if (length_a < 1e-300 || length_b < 1e-300) return 3.14159265358979323846;

  double cosine = dot / (length_a * length_b);
  if (cosine > 1.) cosine = 1.;
  if (cosine < -1.) cosine = -1.;
  return std::acos(cosine);
}

const char* dev_status_name(const unsigned char st)
{
  switch (st) {
    case cxd_dev_separated: return "separated";
    case cxd_dev_enclosed: return "enclosed";
    case cxd_dev_max_iterations: return "max_iterations";
    case cxd_dev_invalid_input: return "invalid_input";
  }
  return "unknown_status";
}

struct golden_row {
  std::string name;
  std::string status;
  bool has_answer = false;
  double depth = 0.;
  std::vector<double> normal;
};

// The committed rows only: comment and blank lines dropped, then one
// "name order dir status [depth normal...]" row each.
std::vector<golden_row> read_golden_rows(const std::string& path)
{
  std::ifstream in(path);
  if (!in) {
    std::cerr << "golden file missing: " << path << '\n';
    std::exit(1);
  }
  std::vector<golden_row> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream parse(line);
    golden_row row;
    int order = 0;
    std::size_t direction = 0;
    std::string depth_token;
    parse >> row.name >> order >> direction >> row.status;
    if (parse >> depth_token) {
      if (depth_token == "refused") {
        row.has_answer = false;
      }
      else {
        row.has_answer = true;
        row.depth = std::stod(depth_token);
        double component = 0.;
        while (parse >> component) row.normal.push_back(component);
      }
    }
    rows.push_back(row);
  }
  return rows;
}

unsigned char status_value(const std::string& name)
{
  if (name == "separated") return cxd_dev_separated;
  if (name == "enclosed") return cxd_dev_enclosed;
  if (name == "max_iterations") return cxd_dev_max_iterations;
  if (name == "invalid_input") return cxd_dev_invalid_input;
  return 255;
}

// One corpus kind (double or float) at one dimension, graded against its
// committed golden file.
template <std::size_t Dim, typename Scalar, bool FloatCorpus>
bool run_corpus(const double depth_tol, const double angle_tol)
{
  constexpr bool float_corpus = FloatCorpus;
  using point = std::array<Scalar, Dim>;
  using case_type = parity_case<point>;

  const std::string dim_tag = float_corpus ? "float" : "double";
  const std::string title = dim_tag + "_" + std::to_string(Dim);

  std::vector<case_type> cases;
  if constexpr (float_corpus)
    cases = make_parity_cases_float<Dim>();
  else
    cases = make_parity_cases<Dim>();

  const std::string path = std::string(CXD_TEST_GOLDEN_DIR) + "/batch_parity_" +
                           dim_tag + "_" + std::to_string(Dim) + ".golden";
  const auto golden = read_golden_rows(path);
  const auto directions = make_parity_directions<Dim>();
  const std::size_t expected_rows = directions.size() * 2 * cases.size();
  if (golden.size() != expected_rows) {
    std::cerr << title << ": row count " << golden.size() << " != expected "
              << expected_rows << '\n';
    return false;
  }

  // The mirrored pool both orders read from: a-shapes, b-shapes, then the
  // same two groups mirrored.
  std::vector<Scalar> flat_vertices;
  std::vector<cxd_dev_entry> entries;
  const auto push_shape = [&](const std::vector<point>& shape) {
    cxd_dev_entry entry{};
    entry.vertex_offset = static_cast<std::uint32_t>(entries.empty()
        ? 0
        : entries.back().vertex_offset + entries.back().vertex_count);
    entry.vertex_count = static_cast<std::uint32_t>(shape.size());
    for (const auto& vertex : shape)
      for (std::size_t k = 0; k < Dim; ++k) flat_vertices.push_back(vertex[k]);
    entries.push_back(entry);
  };
  // Group order a/b/b/a: the forward block reads a at [0,n) and b at
  // [n,2n); the mirrored block reads the same shapes at [2n,3n) and
  // [3n,4n).
  const auto pack_group = [&](const bool take_a) {
    for (const auto& one : cases) push_shape(take_a ? one.a : one.b);
  };
  pack_group(true);
  pack_group(false);
  pack_group(false);
  pack_group(true);

  const auto forward = static_cast<std::uint32_t>(cases.size());

  cxd_device::device_pool<Scalar, Dim> pool;
  pool.upload(flat_vertices.data(), flat_vertices.size(), entries.data(),
              entries.size());

  const cxd_dev_options options{1e-6, 64};

  std::size_t row = 0;
  std::size_t graded_answered = 0;
  std::size_t graded_refused = 0;
  bool ok = true;

  for (const auto& direction : directions)
    for (int order = 0; order < 2; ++order) {
      const bool mirrored = (order == 1);
      const auto left_offset =
          mirrored ? static_cast<std::uint32_t>(2u * cases.size())
                   : static_cast<std::uint32_t>(0);
      const auto right_offset =
          mirrored ? static_cast<std::uint32_t>(3u * cases.size()) : forward;

      std::vector<cxd_dev_pair<Scalar, Dim>> pairs;
      pairs.reserve(cases.size());
      for (std::size_t i = 0; i < cases.size(); ++i) {
        cxd_dev_pair<Scalar, Dim> one{};
        one.shape_id_A = left_offset + static_cast<std::uint32_t>(i);
        one.shape_id_B = right_offset + static_cast<std::uint32_t>(i);
        for (std::size_t k = 0; k < Dim; ++k)
          one.init_dir[k] = mirrored ? Scalar(-direction[k]) : Scalar(direction[k]);
        pairs.push_back(one);
      }

      const auto results = cxd_device::gjk_batch<Scalar, Dim>(pool, pairs, options);

      for (std::size_t i = 0; i < pairs.size(); ++i) {
        // Alignment first: the block's rows must name the same cases in the
        // same order the generator wrote.
        if (golden[row + i].name != cases[i].name) {
          std::cerr << title << ": row " << (row + i) << " names '"
                    << golden[row + i].name << "' but the corpus has '"
                    << cases[i].name << "'\n";
          return false;
        }
        if (golden[row + i].status != dev_status_name(results[i].st)) {
          std::cerr << title << ": row " << (row + i) << " ("
                    << cases[i].name << " order " << order
                    << ", direction " << (&direction - directions.data())
                    << ") golden status " << golden[row + i].status
                    << " but device status " << dev_status_name(results[i].st)
                    << '\n';
          ok = false;
        }
      }

      // Vacuity: a golden sweep of misses keeps nothing graded, the same
      // guard the host run places.
      bool golden_has_hit = false;
      for (std::size_t i = 0; i < pairs.size(); ++i)
        if (golden[row + i].status == "enclosed") golden_has_hit = true;
      if (!golden_has_hit) {
        std::cerr << title << ": block at row " << row
                  << " has no enclosed golden row\n";
        return false;
      }

      std::vector<std::uint32_t> hit_pairs;
      for (std::size_t i = 0; i < results.size(); ++i)
        if (results[i].st == cxd_dev_enclosed)
          hit_pairs.push_back(static_cast<std::uint32_t>(i));

      const auto answers =
          cxd_device::epa_batch<Scalar, Dim>(pool, pairs, hit_pairs, results,
                                             options);

      for (std::size_t h = 0; h < hit_pairs.size(); ++h) {
        const auto& g = golden[row + hit_pairs[h]];
        std::array<double, Dim> got{};
        for (std::size_t k = 0; k < Dim; ++k) got[k] = double(answers[h].normal[k]);
        const double got_depth = double(answers[h].depth);

        if (!g.has_answer) {
          ++graded_refused;
          if (!std::isnan(got_depth)) {
            std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                      << g.name << ") golden refuses the hit but the device"
                      " answers depth " << got_depth << '\n';
            ok = false;
            continue;
          }
          bool all_nan = true;
          for (const double component : got)
            if (!std::isnan(component)) all_nan = false;
          if (!all_nan) {
            std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                      << g.name << ") golden refuses but device normal is"
                      " not all-NaN\n";
            ok = false;
          }
          continue;
        }

        ++graded_answered;
        if (std::isnan(got_depth)) {
          std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                    << g.name << ") golden depth " << g.depth
                    << " but the device answers NaN\n";
          ok = false;
          continue;
        }
        if (got_depth < g.depth - depth_tol || got_depth > g.depth + depth_tol) {
          std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                    << g.name << ") golden depth " << g.depth
                    << " but device depth " << got_depth << '\n';
          ok = false;
        }
        if (g.normal.size() != Dim) {
          std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                    << g.name << ") golden normal has " << g.normal.size()
                    << " components\n";
          return false;
        }
        // A zero normal carries no direction: both sides answer zero and
        // the angle between two zero normals is not defined. Mirror the
        // host's assert_penetration_agrees zero_normal special case.
        bool zero_normal = true;
        for (std::size_t k = 0; k < Dim; ++k)
          if (g.normal[k] != 0.) zero_normal = false;
        if (zero_normal) {
          for (std::size_t k = 0; k < Dim; ++k)
            if (got[k] < -depth_tol || got[k] > depth_tol) {
              std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                        << g.name << ") golden zero normal but device normal "
                        << got[k] << '\n';
              ok = false;
            }
          continue;
        }

        std::array<double, Dim> expected{};
        for (std::size_t k = 0; k < Dim; ++k) expected[k] = g.normal[k];
        const double angle = normal_angle<Dim>(got, expected);
        if (angle <= angle_tol) continue;
        std::cerr << title << ": row " << (row + hit_pairs[h]) << " ("
                  << g.name << ") normal angle " << angle
                  << " rad exceeds the " << angle_tol << " band\n";
        ok = false;
      }

      row += pairs.size();
    }

  // Every hit must be graded: at least one answered hit across the corpus.
  // The host no longer refuses any 2D hit (the degenerate-witness bootstrap
  // answers every case), so the device port matches it point-for-point.
  if (graded_answered == 0) {
    std::cerr << title << ": no answered hit graded\n";
    return false;
  }

  if (ok)
    std::cout << title << ": " << row << " rows graded, "
              << graded_answered << " answered, " << graded_refused
              << " refused\n";
  return ok;
}

}  // namespace

int main()
{
  int device_count = 0;
  if (cudaGetDeviceCount(&device_count) != cudaSuccess || device_count < 1) {
    std::cerr << "no CUDA device visible\n";
    return 1;
  }
  cudaDeviceProp prop{};
  cudaGetDeviceProperties(&prop, 0);
  std::cout << "device: " << prop.name << '\n';

  bool ok = true;
  ok = run_corpus<2, double, false>(1e-6, 1e-4) && ok;
  ok = run_corpus<3, double, false>(1e-6, 1e-4) && ok;
  ok = run_corpus<2, float, true>(1e-4, 1e-3) && ok;
  ok = run_corpus<3, float, true>(1e-4, 1e-3) && ok;

  if (!ok) {
    std::cerr << "device parity FAILED\n";
    return 1;
  }
  std::cout << "device parity ok\n";
  return 0;
}

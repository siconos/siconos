// SPDX-FileCopyrightText: INRIA
//
// SPDX-License-Identifier: Apache-2.0
//
// Siconos is a program dedicated to modeling, simulation and control
// of non smooth dynamical systems.

/*! \file TransportCableManager.h

*/

#pragma once
#include <nlohmann/json.hpp>
#include <string>

#include "TransportCableModel.h"
#include "TransportCableResult.h"

/*! \file TransportCableManager.h

*/

namespace siconos::fem::cable {

/** A class to
 *
 */
class TransportCableManager {
 private:
  std::unique_ptr<TransportCableModel> model_{nullptr};
  TransportCableResult m_results{};

  void computeDS(double a_tolContact = 1e-3, double a_mus = 0.8, double a_mup = 1.1);

  /** Update mass matrix (attribute of m_results)
      \param elem_length elements length
      \param elem_rho linear density
   */
  void compute_mass(double elem_length, double elem_rho);

  /** Update external forces vector (attribute of m_results)
      \param elem_length elements length
      \param elem_rho linear density
  */
  void compute_external_load(double elem_length, double elem_rho);

  // Rule of five
  TransportCableManager(const TransportCableManager&) = delete;
  TransportCableManager(TransportCableManager&&) = delete;
  TransportCableManager& operator=(const TransportCableManager&) = delete;
  TransportCableManager& operator=(TransportCableManager&&) = delete;

 public:
  TransportCableManager() = default;
  TransportCableManager(const std::string& a_filename);

  ~TransportCableManager() noexcept = default;

  /** \return true if the model is valid and the manager ready to work */
  bool isReady() const { return model_->isLoaded(); }

  /** Read a TransportCableModel from a json input
   *  \param a_input a valid json object
   */
  void importModel(const nlohmann::json& a_input);

  /** Read a TransportCableModel from a json file
   *  \param a_filename input file name
   *  \return code error
   */
  void importModel(const std::string& a_filename);

  void computeFEM(const nlohmann::json& a_args, const std::string& a_outfile,
                  nlohmann::ordered_json& output);
  void exportTC(const nlohmann::json& a_args, const std::string& a_outfile,
                nlohmann::ordered_json& output);

  void simulation(const nlohmann::json& a_model, const nlohmann::json& a_args,
                  const std::string& a_filename, const std::string& a_outfile,
                  nlohmann::ordered_json& output);
};
}  // namespace siconos::fem::cable

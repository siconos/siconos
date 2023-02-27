#pragma once
#include "TransportCableModel.h"
#include "TransportCableResult.h"

template <typename T> T getParam(const json &a_arg, const std::string &a_name, T a_default)
{
  T vRet = a_default;
  if (!a_arg.is_null()) {
    if (a_arg.contains(a_name)) {
      a_arg.at(a_name).get_to(vRet);
    }
  }
  return vRet;
}

class TransportCableManager {

private:
  // Rule of five
  TransportCableManager(const TransportCableManager &) = delete;
  TransportCableManager(TransportCableManager &&) = delete;
  TransportCableManager &operator=(const TransportCableManager &) = delete;
  TransportCableManager &operator=(TransportCableManager &&) = delete;

public:
  TransportCableManager() = default;
  TransportCableManager(const std::string &a_filename);

  ~TransportCableManager() noexcept = default;

  int importModel(const json &a_input, const std::string &a_filename = "");
  int computeFEM(const json &a_args, const std::string &a_outfile, ojson &output);
  int exportTC(const json &a_args, const std::string &a_outfile, ojson &output);

  int simulation(const json &a_model, const json &a_args, const std::string &a_filename,
                 const std::string &a_outfile, ojson &output);

private:
  TransportCableModel m_model;
  TransportCableResult m_results;

#ifndef NSICONOS
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
#endif
};

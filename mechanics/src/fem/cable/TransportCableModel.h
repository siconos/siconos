#pragma once
#include "Cable.h"
#include "Carriers.h"
#include "Pile.h"
#include "TransportCable.h"

class TransportCableModel {
public:
  TransportCableModel() = default;

  /** Read data from a json file
      \param a_filename input file name (json)
  */
  TransportCableModel(const std::string &a_filename);

  virtual ~TransportCableModel() noexcept = default;

  int from_file(const std::string &a_fileName);
  int from_json(const json &j);
  bool isLoaded();

  const Cable &get_cable() const;
  const Carriers &get_carriers() const;
  const std::vector<Pile> &get_piles1() const;
  const std::vector<Pile> &get_piles2() const;

private:
  Cable m_cable;
  Carriers m_carriers;
  Pile m_stationUp; // drive station
  Pile m_stationDown;
  std::vector<Pile> m_piles; // Roller batteries
  // Rule of five
  TransportCableModel(const TransportCableModel &) = delete;
  TransportCableModel(TransportCableModel &&) = delete;
  TransportCableModel &operator=(const TransportCableModel &) = delete;
  TransportCableModel &operator=(TransportCableModel &&) = delete;

  void clear();
  int validate();
  std::vector<Pile> m_pilesUp;
  std::vector<Pile> m_pilesDown;
};

#include "TransportCableModel.h"

#include <algorithm>
#include <iostream>

TransportCableModel::TransportCableModel(const std::string &a_filename)
{
  from_file(a_filename);
}

int TransportCableModel::from_file(const std::string &a_fileName)
{
  int res = EXIT_FAILURE;
  json input;
  std::ifstream file(a_fileName);
  if (file.is_open()) {
    file >> input;
    res = from_json(input);
  }
  return res;
}

int TransportCableModel::from_json(const json &j)
{
  int res = EXIT_SUCCESS;
  std::map<std::string, BaseModel &> vElems = {{"cable", m_cable},
                                               {"carriers", m_carriers},
                                               {"stationUp", m_stationUp},
                                               {"stationDown", m_stationDown}};
  clear();
  try {
    for (auto &vElem : vElems) {
      vElem.second.from_json(j, vElem.first);
    }
    if (j.contains("piles")) {
      const json &jpiles = j["piles"];
      if (jpiles.is_array()) {
        m_piles.reserve(jpiles.size());
        for (const auto &jp : jpiles) {
          m_piles.push_back(Pile());
          m_piles.back().from_json(jp);
        }
      }
    }
    res = validate();
  }
  catch (const json::exception &) {
    res = EXIT_FAILURE;
  }

  return res;
}

bool TransportCableModel::isLoaded()
{
  return (m_pilesUp.size() != 0 && m_pilesDown.size() != 0);
}

const Cable &TransportCableModel::get_cable() const { return m_cable; }

const Carriers &TransportCableModel::get_carriers() const { return m_carriers; }

const std::vector<Pile> &TransportCableModel::get_piles1() const { return m_pilesUp; }

const std::vector<Pile> &TransportCableModel::get_piles2() const { return m_pilesDown; }

void TransportCableModel::clear()
{
  // raz du modèle
  m_piles.clear();
  m_pilesUp.clear();
  m_pilesDown.clear();
}

int TransportCableModel::validate()
{
  // Validation du modèle
  // Création des 2 lignes à partir de la définition des poteaux
  if (m_stationUp > m_stationDown) {
    bool vOk = true;
    if (m_piles.size()) {
      // les x des poteaux doivent être croissant
      std::sort(m_piles.begin(), m_piles.end());

      vOk = (m_piles.front() > m_stationDown && m_piles.back() < m_stationUp);
    }
    if (vOk) {
      m_pilesUp.push_back(Pile(m_stationDown, true));
      for (auto &p : m_piles) {
        m_pilesUp.push_back(p); // Copy
      }
      m_pilesUp.push_back(Pile(m_stationUp, true));

      for (auto &p : m_pilesUp) {
        m_pilesDown.push_back(p); // Copy
        p.transform(true);        // useless ???
        m_pilesDown.back().transform(false);
      }
      return EXIT_SUCCESS;
    }
  }

  return EXIT_FAILURE;
}

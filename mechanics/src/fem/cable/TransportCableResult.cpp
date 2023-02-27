#include "TransportCableResult.h"

#ifndef NSICONOS
#include "SiconosVector.hpp"
#include "SimpleMatrix.hpp"
#endif
#include <fstream>

TransportCableResult::TransportCableResult() { rope2.set_Down(true); }


void TransportCableResult::prepareSupport()
{
  supports.clear();
  puller12idx = -1;
  puller21idx = -1;
  rope1.prepareSupport(supports, puller12idx);
  rope2.prepareSupport(supports, puller21idx);

  // Bottom pylon/pulley
  supports[puller12idx]->prepare(rope1.get_LastPile(), rope2.get_LastPile(),
                                 rope1.get_LastT());
  // Top pylon/pulley
  supports[puller21idx]->prepare(rope2.get_FirstPile(), rope1.get_FirstPile(), rope2.get_T0());
}

void TransportCableResult::prepareIneqConstraint(int nb_node)
{
  contacts.clear();
  contacts.resize(nb_node, 0);
  g.clear();
  g.resize(nb_node, 1);
  G.clear();
  G.resize(nb_node);
  for (auto &gg : G) {
    gg.clear();
    gg.resize(nb_node);
  }
  T.clear();
  T.resize(nb_node);
  for (auto &tt : T) {
    tt.clear();
    tt.resize(nb_node);
  }
}

int TransportCableResult::exportTC(const std::string &a_fileName, ojson &a_output,
                                   const std::string &a_option)
{
  int res = EXIT_SUCCESS;
  try {
    res = to_json(a_output, a_option);
  }
  catch (const json::exception &ex) {
    // "Error exporting model " << ex.what());
    throw ex;
  }
  if (a_fileName != "") {
    res = EXIT_FAILURE;
    std::ofstream outputFile(a_fileName);
    if (outputFile.is_open()) {
      outputFile << a_output.dump(2);
      res = EXIT_SUCCESS;
    }
  }
  return res;
}

int TransportCableResult::to_json(ojson &j, const std::string &a_option)
{
  if (a_option == "fem") {
    j["g"] = g;
  }
  else if (a_option == "ropeway") {
    rope1.to_json(j["rope1"]);
    rope2.to_json(j["rope2"]);
  }
  else {
    j["supports"] = ojson::array();
    j["pulleys"] = ojson::array();
    int ns = supports.size();
    for (int i = 0; i < ns; i++) {
      if (i == puller12idx || i == puller21idx) {
        std::shared_ptr<Pulley> pulley = std::static_pointer_cast<Pulley>(supports[i]);
        j["pulleys"].push_back(*pulley);
      }
      else {
        j["supports"].push_back(*supports[i]);
      }
    }
    j["q"] = q;
    j["R"] = R;
    j["tension"] = TS;
    j["punct"] = punct;
    j["cable_length"] = length;
    j["contacts"] = contacts;
  }
  return EXIT_SUCCESS;
}

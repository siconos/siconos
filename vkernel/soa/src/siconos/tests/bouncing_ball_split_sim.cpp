#include "bouncing_ball_split.hpp"

namespace {
auto simulation_handle(siconos::config::data_t& data, std::size_t simulation) {
  return siconos::storage::make_handle(
      data, siconos::storage::index<siconos::config::simulation, std::size_t>(simulation));
}
}  // namespace

void simulation_initialize(siconos::config::data_t& data, std::size_t simulation) {
  simulation_handle(data, simulation).initialize();
}

unsigned simulation_compute_one_step(siconos::config::data_t& data, std::size_t simulation) {
  return simulation_handle(data, simulation).compute_one_step();
}

bool simulation_has_next_event(siconos::config::data_t& data, std::size_t simulation) {
  return simulation_handle(data, simulation).has_next_event();
}

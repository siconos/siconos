#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/moreau_jean_schedule.hpp"
#include "siconos/simul/simul_head.hpp"
#include "siconos/storage/schedule.hpp"

namespace siconos::simul {
template <match::item... Items>
struct time_stepping : item {
  using items = gather<Items...>;
  using time_discretization_t = nth_t<0, items>;
  using one_step_integrator_t = nth_t<1, items>;
  using one_step_nonsmooth_problem_t = nth_t<2, items>;
  using topology_t = typename one_step_integrator_t::topology;

  using formulation_t = typename one_step_nonsmooth_problem_t::problem_t::formulation_t;

  struct attributes {
    some::item_ref<time_discretization_t> time_discretization;
    some::item_ref<one_step_integrator_t> one_step_integrator;
    some::item_ref<one_step_nonsmooth_problem_t> one_step_nonsmooth_problem;
    some::item_ref<topology_t> topology;
  };

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    void __init__() {
      auto& data = self()->data();
      time_discretization() = storage::add<time_discretization_t>(data);
      one_step_integrator() = storage::add<one_step_integrator_t>(data);
      one_step_nonsmooth_problem() = storage::add<one_step_nonsmooth_problem_t>(data);
      topology() = storage::add<topology_t>(data);
    }

    void __del__() {
      // not implemented
      throw std::runtime_error("Stable pointers need to be implemented");
      storage::remove(time_discretization());
      storage::remove(one_step_integrator());
      storage::remove(one_step_nonsmooth_problem());
      storage::remove(topology());
    }

    decltype(auto) time_discretization() {
      return storage::make_ref_handle(self()->data(),
                                      storage::attr<"time_discretization">(*self()));
    }
    decltype(auto) one_step_integrator() {
      return make_ref_handle(self()->data(), storage::attr<"one_step_integrator">(*self()));
    }
    decltype(auto) one_step_nonsmooth_problem() {
      return make_ref_handle(self()->data(),
                             storage::attr<"one_step_nonsmooth_problem">(*self()));
    }
    decltype(auto) topology() {
      return make_ref_handle(self()->data(), storage::attr<"topology">(*self()));
    }

    decltype(auto) current_step() { return (*self()).time_discretization().step(); }

    decltype(auto) time_step() { return self()->time_discretization().h(); }

    // the schedule of a time step, given by the integrator
    using time_step_schedule = typename one_step_integrator_t::time_step_schedule;

    auto compute_one_step() {
      using env_t = decltype(self()->env());
      using indice_t = typename env_t::indice;

      indice_t step = current_step();

      // v, q at step + 1: the free state if there is no contact
      auto ctx = step_context{*self(), one_step_integrator(), step, time_step()};
      storage::schedule::run<time_step_schedule>(ctx);

      current_step() += 1;

      std::print("step {}\n", current_step());
      return ctx.involved;
    }

    // the time step as text, stage by stage
    std::string describe_time_step() {
      return storage::schedule::describe<time_step_schedule>();
    }

    template <typename Formulation>
    void solve_nonsmooth_problem(auto step) {  // for a LCP:
                                               // M z = w + q
                                               // 0 <= z _|_ w >= 0
      auto osi = self()->one_step_integrator();
      //      resize(osi.lambda_vector_assembled(), ninter);
      //      resize(osi.ydot_vector_assembled(), ninter);

      self()->one_step_nonsmooth_problem().template solve<Formulation>(
          osi.w_matrix_assembled(),       // M
          osi.q_nsp_vector_assembled(),   // q
          osi.lambda_vector_assembled(),  // z
          osi.ydot_vector_assembled(),    // w
          osi.mu_vector_assembled());     // mu
    }

    bool has_next_event() {
      return time_discretization().t0() + current_step() * time_discretization().h() <
             time_discretization().tmax();
    }

    void initialize() { one_step_integrator().initialize(current_step()); }

    auto methods() {
      // using env_t = decltype(self()->env());
      // using indice = typename env_t::indice;
      // using scalar = typename env_t::scalar;

      return collect(method("initialize", &interface<Handle>::initialize),
                     method("current_step", &interface<Handle>::current_step),
                     method("time_step", &interface<Handle>::time_step),
                     method("compute_one_step", &interface<Handle>::compute_one_step),
                     method("describe_time_step", &interface<Handle>::describe_time_step),
                     method("has_next_event", &interface<Handle>::has_next_event));
    }
  };
};
}  // namespace siconos::simul

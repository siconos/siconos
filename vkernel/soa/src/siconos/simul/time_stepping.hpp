#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {
template <match::item... Items>
struct time_stepping : item {
  using items = gather<Items...>;
  using time_discretization_t = nth_t<0, items>;
  using one_step_integrator_t = nth_t<1, items>;
  using one_step_nonsmooth_problem_t = nth_t<2, items>;
  using topology_t = typename one_step_integrator_t::topology;

  using formulation_t =
      typename one_step_nonsmooth_problem_t::problem_t::formulation_t;
  using attributes = gather<
      attribute<"time_discretization", some::item_ref<time_discretization_t>>,
      attribute<"one_step_integrator", some::item_ref<one_step_integrator_t>>,
      attribute<"one_step_nonsmooth_problem",
                some::item_ref<one_step_nonsmooth_problem_t>>,
      attribute<"topology", some::item_ref<topology_t>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    void __init__()
    {
      auto& data = self()->data();
      time_discretization() = storage::add<time_discretization_t>(data);
      one_step_integrator() = storage::add<one_step_integrator_t>(data);
      one_step_nonsmooth_problem() =
          storage::add<one_step_nonsmooth_problem_t>(data);
      topology() = storage::add<topology_t>(data);
    }

    void __del__()
    {
      // not implemented
        throw std::runtime_error("Stable pointers need to be implemented");
        storage::remove(time_discretization());
        storage::remove(one_step_integrator());
        storage::remove(one_step_nonsmooth_problem());
        storage::remove(topology());
    }

    decltype(auto) time_discretization()
    {
      return storage::make_ref_handle(
          self()->data(), storage::attr<"time_discretization">(*self()));
    }
    decltype(auto) one_step_integrator()
    {
      return make_ref_handle(self()->data(),
                             storage::attr<"one_step_integrator">(*self()));
    }
    decltype(auto) one_step_nonsmooth_problem()
    {
      return make_ref_handle(
          self()->data(),
          storage::attr<"one_step_nonsmooth_problem">(*self()));
    }
    decltype(auto) topology()
    {
      return make_ref_handle(self()->data(),
                             storage::attr<"topology">(*self()));
    }

    decltype(auto) current_step()
    {
      return (*self()).time_discretization().step();
    }

    decltype(auto) time_step() { return self()->time_discretization().h(); }

    auto compute_one_step()
    {
      using env_t = decltype(self()->env());
      using indice_t = typename env_t::indice;

      auto osi = one_step_integrator();

      indice_t step = current_step();

      indice_t total_number_of_involved_ds = 0;

      // loop over different kind of dynamic systems
      mp::for_each(osi.elements(), [&](auto elem) {
        // udpate M^-1 * fext
        // do nothing if fext is time_invariant
        elem.update_iteration_matrix(step);

        // update jacobians
        // do nothing if lagrangian_r is time_invariant
        elem.update_h_matrices(step);

        // vfree stored in v(step+1)
        elem.compute_free_state(step, time_step());

        // xfree stored in positions(step+1)
        elem.update_positions(step, time_step());

        // -> y & ydot (step & step+1)
        elem.compute_output(step);
        elem.compute_output(step + 1);

        // compute active interactions
        auto [ninter, nds] =
            elem.compute_active_interactions(step, time_step());

        total_number_of_involved_ds += nds;
      });

      if (total_number_of_involved_ds > 0) {
        // resize assembled matrices and vectors
        osi.assemble_setup();

        mp::for_each(osi.elements(), [&](auto elem) {
          // a least one activated interaction
          elem.compute_h_matrices(step + 1);

          elem.assemble_h_matrix_for_involved_ds(step);
          elem.assemble_mass_matrix_for_involved_ds(step);
          elem.assemble_vectors(step);

          elem.nsl_effect_on_free_output(step);
          elem.compute_q_nsp_vector_assembled(step);
        });

        // H M^-1 H^t
        osi.compute_w_matrix(step);
        self()->template solve_nonsmooth_problem<formulation_t>(step);

        // velocity_vector_assembled <- mass_matrix^-1 * (h_matrix^t * lambda)
        osi.compute_input();

        mp::for_each(osi.elements(), [&](auto elem) {
          // lambda_vector_assembled -> lambdas
          elem.keep_lambdas(step);

          // velocity_vector_assembled -> velocities
          elem.update_velocities(step, time_step());
          elem.update_positions(step, time_step());
        });
      }

      else {
        print(".");
      }

      current_step() += 1;

      print("step {}\n", current_step());
      return total_number_of_involved_ds;
    }

    template <typename Formulation>
    void solve_nonsmooth_problem(auto step)
    {  // for a LCP:
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

    bool has_next_event()
    {
      return time_discretization().t0() +
                 current_step() * time_discretization().h() <
             time_discretization().tmax();
    }

    void initialize() { one_step_integrator().initialize(current_step()); }

    auto methods()
    {
      // using env_t = decltype(self()->env());
      // using indice = typename env_t::indice;
      // using scalar = typename env_t::scalar;

      return collect(
          method("initialize", &interface<Handle>::initialize),
          method("current_step", &interface<Handle>::current_step),
          method("time_step", &interface<Handle>::time_step),
          method("compute_one_step", &interface<Handle>::compute_one_step),
          method("has_next_event", &interface<Handle>::has_next_event));
    }
  };
};
}  // namespace siconos::simul

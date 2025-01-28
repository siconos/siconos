#pragma once

#include "siconos/algebra/eigen.hpp"
#include "siconos/model/nslaws.hpp"
#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {
template <match::item... Items>
struct time_stepping : item<> {
  using items = gather<Items...>;
  using time_discretization_t = nth_t<0, items>;
  using one_step_integrator_t = nth_t<1, items>;
  using one_step_nonsmooth_problem_t = nth_t<2, items>;
  using topology_t = typename one_step_integrator_t::topology;

  using formulation_t =
      typename one_step_nonsmooth_problem_t::problem_t::formulation_t;
  using attributes = attributes_of_items<Items...>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) time_discretization()
    {
      return make_handle<time_discretization_t>(*self());
    }
    decltype(auto) one_step_integrator()
    {
      return make_handle<one_step_integrator_t>(*self());
    }
    decltype(auto) one_step_nonsmooth_problem()
    {
      return make_handle<one_step_nonsmooth_problem_t>(*self());
    }
    decltype(auto) topology() { return make_handle<topology_t>(*self()); }

    decltype(auto) current_step()
    {
      return (*self()).time_discretization().step();
    }

    decltype(auto) time_step() { return self()->time_discretization().h(); }

    auto compute_one_step()
    {
      auto osi = one_step_integrator();
      auto step = current_step();

      // update jacobians
      // do nothing if lagrangian_r is time_invariant
      osi.update_h_matrices(step);

      // do nothing if fext is time_invariant
      osi.update_iteration_matrix(step);

      // vfree stored in v(step+1)
      osi.compute_free_state(step, time_step());

      // xfree stored in positions(step+1)
      osi.update_positions(step, time_step());

      // -> y & ydot (step & step+1)
      osi.compute_output(step);
      osi.compute_output(step + 1);

      // compute active interactions
      auto [ninter, nds] = osi.compute_active_interactions(step, time_step());

      if (nds > 0) {
        // a least one activated interaction

        //        print("ninter, nds = {},{}\n", ninter, nds);
        osi.compute_h_matrices(step + 1);
        osi.assemble_h_matrix_for_involved_ds(step, ninter, nds);
        osi.assemble_mass_matrix_for_involved_ds(step, nds);

        osi.resize_assembled_vectors(step, ninter, nds);

        // H M^-1 H^t
        osi.compute_w_matrix(step);
        osi.nsl_effect_on_free_output(step);
        osi.compute_q_nsp_vector_assembled(step, ninter);

        self()->template solve_nonsmooth_problem<formulation_t>(step, ninter);

        // lambda_vector_assembled -> lambdas
        osi.keep_lambdas(step);

        // velocity_vector_assembled <- mass_matrix^-1 * (h_matrix^t * lambda)
        osi.compute_input();

        // velocity_vector_assembled -> velocities
        osi.update_velocities(step, time_step());
        osi.update_positions(step, time_step());
      }
      else {
        print(".");
      }

      current_step() += 1;

      print("step {}\n", current_step());
      return nds;
    }

    template <typename Formulation>
    void solve_nonsmooth_problem(auto step, auto ninter)
    {  // for a LCP:
       // M z = w + q
       //  z _|_ w
      auto osi = self()->one_step_integrator();

      //      resize(osi.lambda_vector_assembled(), ninter);
      //      resize(osi.ydot_vector_assembled(), ninter);

      self()->one_step_nonsmooth_problem().template solve<Formulation>(
          osi.w_matrix(),                 // M
          osi.q_nsp_vector_assembled(),   // q
          osi.lambda_vector_assembled(),  // z
          osi.ydot_vector_assembled());   // w
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
          method("time_discretization",
                 &interface<Handle>::time_discretization),
          method("one_step_integrator",
                 &interface<Handle>::one_step_integrator),
          method("one_step_nonsmooth_problem",
                 &interface<Handle>::one_step_nonsmooth_problem),
          method("topology", &interface<Handle>::topology),
          method("current_step", &interface<Handle>::current_step),
          method("time_step", &interface<Handle>::time_step),
          method("compute_one_step", &interface<Handle>::compute_one_step),
          method("has_next_event", &interface<Handle>::has_next_event));
    }
  };
};
}  // namespace siconos::simul

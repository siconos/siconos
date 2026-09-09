#pragma once

#include <FrictionContactProblem.h>
#include <FrictionContact_options.h>
#include <LinearComplementarityProblem.h>
#include <NonSmoothDrivers.h>
#include <NumericsVerbose.h>
#include <SolverOptions.h>
#include <fclib_interface.h>
#include <lcp_cst.h>

#include <chrono>
#include <format>

#include "siconos/algebra/algebra.hpp"
#include "siconos/simul/simul_head.hpp"

namespace siconos {

namespace simul {

struct solver_options : storage::data_holder<SolverOptions> {
  using attributes = typename storage::data_holder<SolverOptions>::attributes;

  template <typename Handle>
  struct interface : storage::data_holder<SolverOptions>::template interface<
                         Handle> {
    using storage::data_holder<SolverOptions>::template interface<
        Handle>::self;

    void create(int solver_id = SICONOS_FRICTION_2D_LEMKE)
    {
      // need to be fixed : solver_options_delete not called
      self()->instance().reset(solver_options_create(solver_id),
                               [](SolverOptions* so) {
                                 solver_options_delete(so);
                               });
    }

    decltype(auto) iparam(auto index)
    {
      return self()->instance()->iparam[index];
    }
    void set_iparam(auto index, auto value)
    {
      self()->instance()->iparam[index] = value;
    }

    decltype(auto) dparam(auto index)
    {
      return self()->instance()->dparam[index];
    }
    void set_dparam(auto index, auto value)
    {
      self()->instance()->dparam[index] = value;
    }

    decltype(auto) solver_id() { return self()->instance()->solverId; }

    void set_solver_id(auto value) { self()->instance()->solverId = value; }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;
      using scalar = typename env_t::scalar;

      return collect(
          method("create", &interface<Handle>::create),
          method("iparam", &interface<Handle>::iparam<indice>),
          method("dparam", &interface<Handle>::dparam<indice>),
          method("solver_id", &interface<Handle>::solver_id),
          method("set_iparam",
                 &interface<Handle>::set_iparam<indice, indice>),
          method("set_dparam",
                 &interface<Handle>::set_dparam<indice, scalar>),
          method("set_solver_id", &interface<Handle>::set_solver_id<indice>));
    }
  };
};

template <typename Formulation>
struct nonsmooth_problem : storage::data_holder<Formulation> {
  using attributes = typename storage::data_holder<Formulation>::attributes;

  using formulation_t = Formulation;
  template <typename Handle>
  struct interface : storage::data_holder<Formulation>::template interface<
                         Handle> {
    using storage::data_holder<Formulation>::template interface<Handle>::self;

    void create(int solver_id = SICONOS_FRICTION_2D_LEMKE)
    {
      self()->instance().reset(new Formulation);
    };

    //    ~interface() { solver_options_delete(&*_options); };
    auto methods()
    {
      return collect(method("create", &interface<Handle>::create));
    }
  };
};

struct trace_params : item {
  using attributes =
      gather<attribute<"maxiter", some::indice>,
             attribute<"counter", some::indice>,
             attribute<"title", some::specific<std::string>>,
             attribute<"description", some::specific<std::string>>,
             attribute<"math_info", some::specific<std::string>>,
             attribute<"filename", some::specific<std::string>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    void __init__() { counter() = 1; }
    auto maxiter() { return attr<"maxiter">(*self()); }
    decltype(auto) counter() { return attr<"counter">(*self()); }
    auto title() { return attr<"title">(*self()); }
    auto description() { return attr<"description">(*self()); }
    auto math_info() { return attr<"math_info">(*self()); }
    auto filename() { return attr<"filename">(*self()); }
  };
};

template <typename NonsmoothProblem>
struct one_step_nonsmooth_problem : item {
  using problem_t = NonsmoothProblem;
  using attributes =
      gather<attribute<"solver_duration_seconds", some::scalar>,
             attribute<"number_of_contacts", some::indice>,
	     attribute<"number_of_iterations", some::indice>,
             attribute<"level", some::indice>,
             attribute<"trace", some::boolean>,
             attribute<"trace_params", some::item_ref<trace_params>>,
             attribute<"verbose", some::boolean>,
             attribute<"options", some::item_ref<solver_options>>,
             attribute<"problem", some::item_ref<NonsmoothProblem>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) trace() { return attr<"trace">(*self()); }
    decltype(auto) trace_params()
    {
      return storage::make_ref_handle(self()->data(),
                                      attr<"trace_params">(*self()));
    }

    decltype(auto) options()
    {
      return storage::make_ref_handle(self()->data(),
                                      attr<"options">(*self()));
    }
    decltype(auto) problem()
    {
      return storage::make_ref_handle(self()->data(),
                                      attr<"problem">(*self()));
    }
    decltype(auto) level() { return attr<"level">(*self()); };

    template <typename Formulation, match::matrix W, match::vector V>
    void solve(algebra::mat<W>& w_mat, algebra::vec<V>& q_vec,
               algebra::vec<V>& z_vec, algebra::vec<V>& w_vec,
               algebra::vec<V>& mu_vec)  // mu_vec can be empty for LCP
    {
      using env_t = decltype(self()->env());
      using scalar = typename env_t::scalar;

      numerics_set_verbose(storage::attr<"verbose">(*self()));
      if constexpr (std::derived_from<Formulation,
                                      LinearComplementarityProblem>) {
        // w_mat cannot be sparse
        auto w_mat_dense = NM_create(NM_DENSE, size0(w_mat), size1(w_mat));
        NM_to_dense(w_mat._m, w_mat_dense);

        self()->problem().instance()->size = size0(w_mat);
        self()->problem().instance()->M = w_mat_dense;
        self()->problem().instance()->q = q_vec._v->matrix0;

        linearComplementarity_driver(&*self()->problem().instance(),
                                     z_vec._v->matrix0, w_vec._v->matrix0,
                                     &*options().instance());

        NM_free(w_mat_dense);
      }
      else if constexpr (std::derived_from<Formulation,
                                           FrictionContactProblem>) {
        self()->problem().instance()->dimension = 2;
        self()->problem().instance()->numberOfContacts =
            size0(w_mat) / self()->problem().instance()->dimension;
        self()->problem().instance()->M = w_mat._m;
        self()->problem().instance()->q = q_vec._v->matrix0;
        self()->problem().instance()->mu = mu_vec._v->matrix0;

        attr<"number_of_contacts">(*self()) =
            self()->problem().instance()->numberOfContacts;

        if (!trace()) {
          const auto start_time{std::chrono::steady_clock::now()};
          fc2d_driver(&*self()->problem().instance(), z_vec._v->matrix0,
                      w_vec._v->matrix0, &*options().instance());
          const auto end_time{std::chrono::steady_clock::now()};

          const std::chrono::duration<scalar> elapsed_seconds{end_time -
                                                              start_time};

          attr<"solver_duration_seconds">(*self()) = elapsed_seconds.count();

        }
        else {
          auto z_bck = algebra::copy(z_vec);
          auto w_bck = algebra::copy(w_vec);

          const auto start_time{std::chrono::steady_clock::now()};
          fc2d_driver(&*self()->problem().instance(), z_vec._v->matrix0,
                      w_vec._v->matrix0, &*options().instance());
          const auto end_time{std::chrono::steady_clock::now()};

          const std::chrono::duration<scalar> elapsed_seconds{end_time -
                                                              start_time};
          attr<"solver_duration_seconds">(*self()) = elapsed_seconds.count();

	  attr<"number_of_iterations">(*self()) = options().iparam(SICONOS_IPARAM_ITER_DONE);

          // trace_params must be set!
          if (options().iparam(SICONOS_IPARAM_ITER_DONE) >
              trace_params().maxiter()) {
            auto solver_maxiter = options().iparam(SICONOS_IPARAM_MAX_ITER);
            auto n_format_string = std::to_string(solver_maxiter).length();

            // format alone ambiguous for clang-19
            auto iter_filename = std::format(
                "{}-i{:0{}d}-{}-{}.hdf5", trace_params().filename(),
                solver_maxiter, n_format_string, size0(w_mat),
                trace_params().counter());

            frictionContact_fclib_write(&*self()->problem().instance(),
                                        trace_params().title().c_str(),
                                        trace_params().description().c_str(),
                                        trace_params().math_info().c_str(),
                                        iter_filename.c_str(),
                                        3  // dof
            );

            frictionContact_fclib_write_guess(
                z_vec._v->matrix0, w_vec._v->matrix0,
                trace_params().filename().c_str());

            trace_params().counter()++;
          }
        }
      }
    }
    auto methods()
    {
      return collect(method("options", &interface<Handle>::options),
                     method("problem", &interface<Handle>::problem),
                     method("level", &interface<Handle>::level));
    }
  };
};
}  // namespace simul
}  // namespace siconos

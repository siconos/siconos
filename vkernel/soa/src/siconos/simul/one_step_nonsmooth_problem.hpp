#pragma once

#include <FrictionContactProblem.h>
#include <Friction_cst.h>
#include <LinearComplementarityProblem.h>
#include <NonSmoothDrivers.h>
#include <NumericsVerbose.h>
#include <SolverOptions.h>
#include <lcp_cst.h>

#include "siconos/algebra/algebra.hpp"
#include "siconos/simul/simul_head.hpp"

namespace siconos {

namespace simul {

struct solver_options : storage::data_holder<SolverOptions> {
  using attributes = typename storage::data_holder<SolverOptions>::attributes;

  template <typename Handle>
  struct interface : storage::data_holder<SolverOptions>::template interface<
                         Handle> {
    using default_interface<Handle>::self;
    void create(int solver_id = SICONOS_FRICTION_2D_LEMKE)
    {
      self()->instance().reset(solver_options_create(solver_id),
                               [](SolverOptions* so) {
                                 solver_options_delete(so);
                                 delete so;
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
    using default_interface<Handle>::self;

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

template <typename NonsmoothProblem>
struct one_step_nonsmooth_problem : item<> {
  using problem_t = NonsmoothProblem;
  using attributes =
      gather<attribute<"level", some::indice>,
             attribute<"verbose", some::boolean>,
             attribute<"options", some::item_ref<solver_options>>,
             attribute<"problem", some::item_ref<NonsmoothProblem>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) options()
    {
      return storage::handle(self()->data(), attr<"options">(*self()));
    };
    decltype(auto) problem()
    {
      return storage::handle(self()->data(), attr<"problem">(*self()));
    };
    decltype(auto) level() { return attr<"level">(*self()); };

    template <typename Formulation, match::matrix W, match::vector V>
    void solve(algebra::mat<W>& w_mat, algebra::vec<V>& q_vec,
               algebra::vec<V>& z_vec, algebra::vec<V>& w_vec,
               algebra::vec<V>& mu_vec) // mu_vec can be empty for LCP
    {
      using fmt::print;

      numerics_set_verbose(storage::attr<"verbose">(*self()));
      if constexpr (std::derived_from<Formulation,
                                      LinearComplementarityProblem>) {
        // w_mat cannot be sparse
        auto w_mat_dense = NM_create(NM_DENSE, size0(w_mat), size1(w_mat));
        NM_to_dense(w_mat._m, w_mat_dense);

        /*        print("w_mat_dense:\n");
                  NM_display(w_mat_dense);*/
        self()->problem().instance()->size = size0(w_mat);
        self()->problem().instance()->M = w_mat_dense;
        self()->problem().instance()->q = q_vec._v->matrix0;

        // print("LCP [\n");

        // print("W:\n");
        // algebra::display(w_mat);
        // print("----\n");
        // print("----\n");

        // print("q:\n");
        // algebra::display(q_vec);
        // print("----\n");

        // print("z:\n");
        // algebra::display(z_vec);
        // print("----\n");

        // print("w:\n");
        // algebra::display(w_vec);
        // print("----\n");

        linearComplementarity_driver(&*self()->problem().instance(),
                                     z_vec._v->matrix0, w_vec._v->matrix0,
                                     &*options().instance());

        // print("q:\n");
        // algebra::display(q_vec);
        // print("----\n");

        // print("z:\n");
        // algebra::display(z_vec);
        // print("----\n");

        // print("w:\n");
        // algebra::display(w_vec);
        // print("----\n");

        // print("]\n\n");
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

        fc2d_driver(&*self()->problem().instance(), z_vec._v->matrix0,
                    w_vec._v->matrix0, &*options().instance());

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

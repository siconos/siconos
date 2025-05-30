#pragma once

#include "siconos/simul/simul_head.hpp"
#include "siconos/algebra/algebra.hpp"

namespace siconos::simul {

struct moreau_jean_assembled : item<> {
  using size_1_t = some::indice_value<1>;

  using attributes = gather<
      attribute<"theta", some::scalar>, attribute<"gamma", some::scalar>,
      attribute<"constraint_activation_threshold", some::scalar>,
      attribute<"mass_matrix_assembled",
                some::assembled_matrix<
                    some::matrix<some::scalar, size_1_t, size_1_t>>>,
      attribute<"h_matrix_assembled", some::assembled_matrix<some::matrix<
                                          some::scalar, size_1_t, size_1_t>>>,
      attribute<"w_matrix_assembled", some::assembled_matrix<some::matrix<
                                          some::scalar, size_1_t, size_1_t>>>,
      attribute<"lambda_vector_assembled",
                some::assembled_vector<some::vector<some::scalar, size_1_t>>>,
      attribute<"y_vector_assembled",
                some::assembled_vector<some::vector<some::scalar, size_1_t>>>,
      attribute<"ydot_vector_assembled",
                some::assembled_vector<some::vector<some::scalar, size_1_t>>>,
      attribute<"velocity_vector_assembled",
                some::assembled_vector<some::vector<some::scalar, size_1_t>>>,
      attribute<"p0_vector_assembled",
                some::assembled_vector<some::vector<some::scalar, size_1_t>>>,
      attribute<"q_nsp_vector_assembled", some::assembled_vector<some::vector<
                                              some::scalar, size_1_t>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mass_matrix_assembled()
    {
      return storage::attr<"mass_matrix_assembled">(*self());
    }

    decltype(auto) h_matrix_assembled()
    {
      return storage::attr<"h_matrix_assembled">(*self());
    }

    decltype(auto) w_matrix_assembled()
    {
      return storage::attr<"w_matrix_assembled">(*self());
    }

    decltype(auto) minv_fext_vector_assembled()
    {
      return storage::attr<"minv_fext_vector_assembled">(*self());
    }
    decltype(auto) lambda_vector_assembled()
    {
      return storage::attr<"lambda_vector_assembled">(*self());
    }

    decltype(auto) velocity_vector_assembled()
    {
      return storage::attr<"velocity_vector_assembled">(*self());
    }

    decltype(auto) p0_vector_assembled()
    {
      return storage::attr<"p0_vector_assembled">(*self());
    }
    decltype(auto) q_nsp_vector_assembled()
    {
      return attr<"q_nsp_vector_assembled">(*self());
    }

    decltype(auto) y_vector_assembled()
    {
      return attr<"y_vector_assembled">(*self());
    }
    decltype(auto) ydot_vector_assembled()
    {
      return attr<"ydot_vector_assembled">(*self());
    }

    decltype(auto) theta() { return attr<"theta">(*self()); }
    decltype(auto) gamma() { return attr<"gamma">(*self()); }
    decltype(auto) constraint_activation_threshold()
    {
      return attr<"constraint_activation_threshold">(*self());
    }

    void assemble_setup(auto raw_inter_size, auto raw_ds_size)
    {
      algebra::resize(mass_matrix_assembled(), raw_ds_size, raw_ds_size);
      algebra::resize(h_matrix_assembled(), raw_inter_size, raw_ds_size);
      algebra::resize(lambda_vector_assembled(), raw_inter_size);
      algebra::resize(velocity_vector_assembled(), raw_ds_size);
      algebra::resize(y_vector_assembled(), raw_inter_size);
      algebra::resize(ydot_vector_assembled(), raw_inter_size);
      algebra::resize(q_nsp_vector_assembled(), raw_inter_size);
    }

    void compute_w_matrix(auto step)
    {
      algebra::compute_kkt_matrix(h_matrix_assembled(),
                                  mass_matrix_assembled(),
                                  w_matrix_assembled());
    }

    void compute_input()
    {
      auto &h_matrix = h_matrix_assembled();
      auto &lambda = lambda_vector_assembled();
      auto &p0 = p0_vector_assembled();
      auto &velo = velocity_vector_assembled();
      auto &mass_matrix = mass_matrix_assembled();

      resize(p0, size1(h_matrix));
      resize(velo, size1(h_matrix));

      transpose(h_matrix);
      prodt1(h_matrix, lambda, p0);
      solve(mass_matrix, p0, velo);
    }

    auto methods()
    {
      using env_t = decltype(self()->env());
      using indice = typename env_t::indice;

      return collect(
          method("compute_input", &interface<Handle>::compute_input),
          method("compute_w_matrix",
                 &interface<Handle>::compute_w_matrix<indice>));
    }
  };
};
}  // namespace siconos::simul

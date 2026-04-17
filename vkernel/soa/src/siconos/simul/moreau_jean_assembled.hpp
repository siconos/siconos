#pragma once

#include "siconos/algebra/algebra.hpp"
#include "siconos/simul/simul_head.hpp"

namespace siconos::simul {

struct moreau_jean_assembled : item {
  using size_1_t = some::indice_value<1>;

  using attributes = gather<
      attribute<"theta", some::scalar>, attribute<"gamma", some::scalar>,
      attribute<"constraint_activation_threshold", some::scalar>,
      attribute<"mass_matrix_assembled",
                some::assembled_matrix<
                    some::matrix<some::scalar, size_1_t, size_1_t>>>,
      attribute<"iteration_matrix_assembled",
                some::assembled_matrix<
                    some::matrix<some::scalar, size_1_t, size_1_t>>>,
      attribute<"k_matrix_assembled", some::assembled_matrix<some::matrix<
                                          some::scalar, size_1_t, size_1_t>>>,
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
      attribute<"q_nsp_vector_assembled",
                some::assembled_vector<some::vector<some::scalar, size_1_t>>>,
      attribute<"mu_vector_assembled", some::assembled_vector<some::vector<
                                           some::scalar, size_1_t>>>>;

  template <typename Handle>
  struct interface : default_interface<Handle> {
    using default_interface<Handle>::self;

    decltype(auto) mass_matrix_assembled()
    {
      return storage::attr<"mass_matrix_assembled">(*self());
    }

    decltype(auto) iteration_matrix_assembled()
    {
      return storage::attr<"iteration_matrix_assembled">(*self());
    }

    decltype(auto) k_matrix_assembled()
    {
      return storage::attr<"k_matrix_assembled">(*self());
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

    decltype(auto) mu_vector_assembled()
    {
      return attr<"mu_vector_assembled">(*self());
    }

    decltype(auto) theta() { return attr<"theta">(*self()); }
    decltype(auto) gamma() { return attr<"gamma">(*self()); }
    decltype(auto) constraint_activation_threshold()
    {
      return attr<"constraint_activation_threshold">(*self());
    }
  };
};
}  // namespace siconos::simul

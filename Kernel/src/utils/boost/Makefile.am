# Siconos-Kernel version 2.1.1, Copyright INRIA 2005-2007.
# Siconos is a program dedicated to modeling, simulation and control
# of non smooth dynamical systems.	
# Siconos is a free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# Siconos is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Siconos; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contact: Vincent ACARY vincent.acary@inrialpes.fr 
#	
# Makefile.am --- Boost bindings

include $(top_srcdir)/config/Make-rules

atlasdir = $(boostdir)/numeric/bindings/atlas
traitsdir = $(boostdir)/numeric/bindings/traits
detaildir = $(boostdir)/numeric/bindings/traits/detail

atlas_HEADERS = numeric/bindings/atlas/cblas1.hpp numeric/bindings/atlas/cblas2_overloads.hpp numeric/bindings/atlas/cblas_enum.hpp numeric/bindings/atlas/clapack.hpp numeric/bindings/atlas/cblas1_overloads.hpp \
		  numeric/bindings/atlas/cblas3.hpp numeric/bindings/atlas/cblas.hpp numeric/bindings/atlas/clapack_inc.hpp numeric/bindings/atlas/cblas2.hpp numeric/bindings/atlas/cblas3_overloads.hpp \
		  numeric/bindings/atlas/cblas_inc.hpp numeric/bindings/atlas/clapack_overloads.hpp

traits_HEADERS = numeric/bindings/traits/algorithm.hpp numeric/bindings/traits/sparse_traits.hpp numeric/bindings/traits/traits.hpp numeric/bindings/traits/ublas_banded.hpp numeric/bindings/traits/ublas_vector2.hpp \
		  numeric/bindings/traits/boost_array.hpp numeric/bindings/traits/fortran.h numeric/bindings/traits/std_valarray.hpp numeric/bindings/traits/transpose.hpp numeric/bindings/traits/ublas_hermitian.hpp numeric/bindings/traits/ublas_vector.hpp \
		  numeric/bindings/traits/c_array.hpp numeric/bindings/traits/std_vector.hpp numeric/bindings/traits/type.h numeric/bindings/traits/ublas_matrix.hpp numeric/bindings/traits/vector_raw.hpp numeric/bindings/traits/config.hpp numeric/bindings/traits/matrix_raw.hpp \
		  numeric/bindings/traits/symm_herm_raw.hpp numeric/bindings/traits/type.hpp numeric/bindings/traits/ublas_sparse.hpp numeric/bindings/traits/vector_traits.hpp numeric/bindings/traits/matrix_traits.hpp numeric/bindings/traits/tnt.hpp numeric/bindings/traits/type_traits.hpp numeric/bindings/traits/ublas_symmetric.hpp

detail_HEADERS = numeric/bindings/traits/detail/array.hpp numeric/bindings/traits/detail/symm_herm_traits.hpp numeric/bindings/traits/detail/ublas_uplo.hpp \
		 numeric/bindings/traits/detail/array_impl.hpp numeric/bindings/traits/detail/generate_const.hpp \
		 numeric/bindings/traits/detail/ublas_ordering.hpp numeric/bindings/traits/detail/utils.hpp


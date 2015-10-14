//
//  Copyright Markus Rickert 2008
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/vector_proxy.hpp>
#include <boost/numeric/bindings/blas/level3/gemm.hpp>
#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/bindings/std/valarray.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace bindings = boost::numeric::bindings;

int
main(int argc, char** argv)
{
	{
		// a * b' = C ; a' * b = d
		
		boost::numeric::ublas::vector<double> a(3);
		for (std::size_t i = 0; i < a.size(); ++i) a(i) = i;
		std::cout << "a=" << a << std::endl;
		
		boost::numeric::ublas::vector<double> b(3);
		for (std::size_t i = 0; i < b.size(); ++i) b(i) = i;
		std::cout << "b=" << b << std::endl;
		
		boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> c(3, 3);
		boost::numeric::bindings::blas::gemm(
			1.0, a, bindings::trans(b), 0.0, c
		);
		std::cout << "c=" << c << std::endl;
		
		boost::numeric::ublas::vector<double> d(1);
		boost::numeric::bindings::blas::gemm(
			1.0, bindings::trans(a), b, 0.0, d
		);
		std::cout << "d=" << d << std::endl;
	}
	
	std::cout << std::endl;
	
	{
		// a * b' = C ; a' * b = d
		
		std::vector<double> a(3);
		for (std::size_t i = 0; i < a.size(); ++i) a[i] = i;
		std::cout << "a=[" << a.size() << "](";
		for (std::size_t i = 0; i < a.size(); ++i) std::cout << (i > 0 ? "," : "") << a[i];
		std::cout << ")" << std::endl;
		
		std::valarray<double> b(3);
		for (std::size_t i = 0; i < b.size(); ++i) b[i] = i;
		std::cout << "b=[" << b.size() << "](";
		for (std::size_t i = 0; i < b.size(); ++i) std::cout << (i > 0 ? "," : "") << b[i];
		std::cout << ")" << std::endl;
		
		boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> c(3, 3);
		boost::numeric::bindings::blas::gemm(
			1.0, a, bindings::trans(b), 0.0, c
		);
		std::cout << "c=" << c << std::endl;
		
		std::vector<double> d(1);
		boost::numeric::bindings::blas::gemm(
			1.0, bindings::trans(a), b, 0.0, d
		);
		std::cout << "d=[" << d.size() << "](";
		for (std::size_t i = 0; i < d.size(); ++i) std::cout << (i > 0 ? "," : "") << d[i];
		std::cout << ")" << std::endl;
	}
	
	std::cout << std::endl;
	
	{
		// a * b' = C ; a' * b = d
		
		double a[3];
		for (std::size_t i = 0; i < 3; ++i) a[i] = i;
		std::cout << "a=[" << 3 << "](";
		for (std::size_t i = 0; i < 3; ++i) std::cout << (i > 0 ? "," : "") << a[i];
		std::cout << ")" << std::endl;
		
		std::vector<double> b(3);
		for (std::size_t i = 0; i < b.size(); ++i) b[i] = i;
		std::cout << "b=[" << b.size() << "](";
		for (std::size_t i = 0; i < b.size(); ++i) std::cout << (i > 0 ? "," : "") << b[i];
		std::cout << ")" << std::endl;
		
		boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> c(3, 3);
		boost::numeric::bindings::blas::gemm(
			1.0, a, bindings::trans(b), 0.0, c
		);
		std::cout << "c=" << c << std::endl;
		
		std::valarray<double> d(1);
		boost::numeric::bindings::blas::gemm(
			1.0, bindings::trans(a), b, 0.0, d
		);
		std::cout << "d=[" << d.size() << "](";
		for (std::size_t i = 0; i < d.size(); ++i) std::cout << (i > 0 ? "," : "") << d[i];
		std::cout << ")" << std::endl;
	}
	
	return 0;
}


//
//  Copyright (c) 2007
//  Tsai, Dung-Bang	
//  National Taiwan University, Department of Physics
// 
//  E-Mail : dbtsai (at) gmail.com
//  Begine : 2007/11/20
//  Last modify : 2007/11/22
//  Version : v0.1
//
//  EXPGM_PAD computes the matrix exponential exp(H) for general matrixs,
//  including complex and real matrixs using the irreducible (p,p) degree
//  rational Pade approximation to the exponential 
//  exp(z) = r(z)=(+/-)( I+2*(Q(z)/P(z))).
//
//  Usage : 
//
//    U = expm_pad(H)
//    U = expm_pad(H, p)
//    
//    where p is internally set to 6 (recommended and gererally satisfactory).
//
//  See also MATLAB supplied functions, EXPM and EXPM1.
//
//  Reference :
//  EXPOKIT, Software Package for Computing Matrix Exponentials.
//  ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
//
//  Permission to use, copy, modify, distribute and sell this software
//  and its documentation for any purpose is hereby granted without fee,
//  provided that the above copyright notice appear in all copies and
//  that both that copyright notice and this permission notice appear
//  in supporting documentation.  The authors make no representations
//  about the suitability of this software for any purpose.
//  It is provided "as is" without express or implied warranty.
//

#ifndef _BOOST_UBLAS_EXPM_
#define _BOOST_UBLAS_EXPM_
/// @cond

#include <complex>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace boost { namespace numeric { namespace ublas {

template<typename MATRIX> MATRIX expm_pad(const MATRIX &H, const unsigned int p = 13)
{
	typedef typename MATRIX::value_type value_type;
        typedef typename MATRIX::size_type size_type;
	typedef double real_value_type;	// Correct me. Need to modify.
	assert(H.size1() == H.size2());	
	const size_type n = H.size1();
	const identity_matrix<value_type> I(n);
	matrix<value_type> U(n,n),H2(n,n),P(n,n),Q(n,n);
	real_value_type norm = 0.0;

// Calcuate Pade coefficients  (1-based instead of 0-based as in the c vector)
	vector<real_value_type> c(p+2);
	c(1)=1;  
	for(size_type i = 1; i <= p; ++i)
		c(i+1) = c(i) * ((p + 1.0 - i)/(i * (2.0 * p + 1 - i)));
// Calcuate the infinty norm of H, which is defined as the largest row sum of a matrix
	for(size_type i=0; i<n; ++i)
	{
		real_value_type temp = 0.0;
		for(size_type j=0;j<n;j++)
			temp += std::abs<real_value_type>(H(i,j)); // Correct me, if H is complex, can I use that abs?
		norm = std::max<real_value_type>(norm, temp);
	}
	if (norm == 0.0) 
	{
		std::cerr<<"Error! Null input in the routine EXPM_PAD.\n";
		exit(0);
	}
// Scaling, seek s such that || H*2^(-s) || < 1/2, and set scale = 2^(-s)
 	size_type s = 0;
	real_value_type scale = 1.0;

  
	//if(norm > 0.5)
	{
		s = std::max<int>(0, static_cast<int>((log(norm) / log(2.0) + 2.0)));
		scale /= static_cast<real_value_type>(std::pow(2.0, (double)s));
		U.assign(scale * H); // Here U is used as temp value due to that H is const
	}
// Horner evaluation of the irreducible fraction, see the following ref above.
// Initialise P (numerator) and Q (denominator) 
	H2.assign( prod(U, U) );
	Q.assign( c(p+1)*I );
	P.assign( c(p)*I );
	size_type odd = 1;
	for( size_type k = p - 1; k > 0; --k)
	{
		if( odd == 1)
		{
			Q = ( prod(Q, H2) + c(k) * I ); 
		}
		else
		{
			P = ( prod(P, H2) + c(k) * I );
		}
		odd = 1 - odd;
	}
	if( odd == 1)
	{
		Q = ( prod(Q, U) );	
		Q -= P ;
		//U.assign( -(I + 2*(Q\P)));
	}
	else
	{
		P = (prod(P, U));
		Q -= P;
		//U.assign( I + 2*(Q\P));
	}
// In origine expokit package, they use lapack ZGESV to obtain inverse matrix,
// and in that ZGESV routine, it uses LU decomposition for obtaing inverse matrix.
// Since in ublas, there is no matrix inversion template, I simply use the build-in
// LU decompostion package in ublas, and back substitute by myself.
//
//////////////// Implement Matrix Inversion ///////////////////////
	permutation_matrix<size_type> pm(n); 
	int res = lu_factorize(Q, pm);
	if( res != 0)
	{
		std::cerr << "Error in the matrix inversion in template expm_pad.\n";
		exit(0);
	}
	H2 = I;  // H2 is not needed anymore, so it is temporary used as identity matrix for substituting.
	lu_substitute(Q, pm, H2); 
	if( odd == 1)
		U.assign( -(I + 2.0 * prod(H2, P)));
 	else
		U.assign( I + 2.0 * prod(H2, P));
// Squaring 
	for(size_type i = 0; i < s; ++i)
	{
		U = (prod(U,U));
	}
	return U;
}

}}}

/// @endcond
#endif

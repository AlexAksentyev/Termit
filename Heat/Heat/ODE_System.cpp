#include "stdafx.h"
#include "ODE_System.h"
#include <boost\numeric\odeint.hpp>


bool ODE_System::integrate(boost::numeric::ublas::vector<temperature>& initials)
{
	using namespace boost::numeric::odeint;
	
	integrate_const( runge_kutta4<vector>() , this->rhs , initials , this->parameters.time_span[0] , this->parameters.time_span[1] , this->parameters.time_step, Observer( std::ofstream(this->parameters.output_file) ) ); // for now it outputs into a file only

	return true;
}

bool ODE_System::InvertMatrix(const matrix& input, boost::numeric::ublas::compressed_matrix<mx_elem>& inverse)
{
	typedef boost::numeric::ublas::permutation_matrix<std::size_t> pmatrix;

	// create a working copy of the input
	matrix A(input);

	// create a permutation matrix for the LU-factorization
	pmatrix pm(A.size1());

	// perform LU-factorization
	int res = lu_factorize(A, pm);
	if (res != 0)
		return false;

	// create identity matrix of "inverse"
	inverse.assign(boost::numeric::ublas::identity_matrix<mx_elem> (A.size1()));

	// backsubstitute to get the inverse
	lu_substitute(A, pm, inverse);

	return true;
}